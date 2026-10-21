"""Compare the uniform revised preparation with immutable raw and v1 evidence."""
import csv
import json
import numpy as np
from prepare_controls import ROOT, read_structure
from analyze_controls import pose_metrics
from analyze_generation_comparison import design_pose_metrics, contacts_summary, screen_sample
from direction_check import histidines, release_ratio
from summarize_direction import terminal_distance
from receptor_preparation_v2 import sha
from audit_expansion_his_context import distance


def direction_metrics(free, bound, terminus_distance):
    shifts = [b-f for f,b in zip(free,bound)]
    ratio = release_ratio(free,bound)
    checks = {'all_shifts_le_minus_0_5': all(d <= -.5 for d in shifts),
              'any_free_pka_6_to_7': any(6 <= f <= 7 for f in free),
              'his_only_ratio_ge_2': ratio >= 2, 'target_terminus_ge_10A': terminus_distance >= 10}
    return {'free_pkas': free, 'bound_pkas': bound, 'shifts': shifts,
            'his_only_Kd6_over_Kd7_4': ratio, 'artificial_terminus_distance_A': terminus_distance,
            'checks': checks, 'direction_pass': all(checks.values())}


def main():
    out = ROOT/'reports/receptor_preparation_v2'
    prep = json.loads((out/'results.json').read_text())
    rawpath = ROOT/'reports/expansion_his_confirmation/results.json'
    raw = json.loads(rawpath.read_text())
    assert prep['source_results_sha256'] == sha(rawpath)
    selected = {(r['id'],r['seed'],r['sample']): r for r in raw['rows'] if r['confidence_selected']}
    assert len(prep['rows']) == len(selected) == 28
    assert {(r['id'],r['seed'],r['sample']) for r in prep['rows']} == set(selected)
    manifest = json.loads((ROOT/'data/expansion_his_confirmation/manifest.json').read_text())
    cases = {c['id']: c for c in manifest['cases']}; arms = {a['id']: a for a in manifest['arms']}
    conserved = {int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    v1 = json.loads((ROOT/'reports/expansion_his_direction/results.json').read_text())
    original_direction = {(v['id'],m['seed']): m for v in v1['variants'] for m in v['models']}
    original_both = {v['id']:v['all_selected_seeds_full_screen'] for v in raw['candidate_summaries']}
    rows = []
    for p in sorted(prep['rows'],key=lambda x:(x['id'],x['seed'])):
        previous = selected[(p['id'],p['seed'],p['sample'])]; case = cases[p['id']]
        r = {'id':p['id'],'seed':p['seed'],'sample':p['sample'],'kind':case['kind'],
             'raw_pass':previous.get('screen',{}).get('pass',previous.get('control_pass')),
             'preparation_status':p['status'], 'prepared_eligible':False}
        assert sha(ROOT/p['source']) == p['source_sha256'] == prep['source_hashes'][p['source']]
        if p['status'] != 'ok':
            r['error'] = p['error']; rows.append(r); continue
        path = ROOT/p['dest']/'prepared.pdb'; assert sha(path) == p['prepared_sha256']
        s = read_structure(path); ref = read_structure(ROOT/'data/expansion_his_confirmation'/f"{p['id']}.pdb")
        r.update(displacement=p['displacement'], internal_geometry=p['prepared_internal_geometry'], propka=p['propka'],
                 prepared_source=str(path.relative_to(ROOT)),prepared_sha256=sha(path))
        if case['kind'] == 'novel_design':
            scored = {'metrics':previous['metrics'], 'pose':design_pose_metrics(s,ref), 'geometry':p['prepared_geometry'],
                      'contacts':contacts_summary(s,arms[case['arm']],conserved)}
            r['prepared_screen'] = screen_sample(scored,manifest['screen']); r['prepared_pose'] = scored['pose']
        elif case['kind'] == 'positive':
            pose = pose_metrics(s,ref); r['prepared_pose'] = pose
            checks = {'original_confidence':previous['metrics']['design_to_target_iptm'] >= .5,
                      'fold':pose['binder_CA_fold_rmsd_A'] <= 2,
                      'interface':pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and pose['binder_interface_CA_rmsd_after_target_fit_A'] <= 2,
                      'recall':pose['native_contact_recall'] >= .5}
            r['prepared_screen'] = {'checks':checks,'pass':all(checks.values())}
        else:
            r['prepared_screen'] = {'pass':previous['control_pass'], 'note':'Original presumed-negative confidence check retained; minimization cannot recompute prediction confidence or prove nonbinding.'}
        r['prepared_eligible'] = bool(r['raw_pass'] and p['displacement']['pass'] and p['prepared_internal_geometry']['pass'] and
                                     r['prepared_screen']['pass'] and all(x['status']=='ok' for x in p['propka'].values()))
        binder = list(s[0]['B']); target_atoms = [a for a in s[0]['A'].get_atoms() if a.element!='H']
        terminus = [a for a in binder[-1] if a.name in {'C','O','OXT'}]
        r['tag_context'] = {'binder_length':len(binder), 'terminal_carboxyl_to_target_A':distance(terminus,target_atoms),
                            'introduced_histidines':[], 'interpretation':'Untagged coordinates only. Distance does not establish tag accessibility or binding compatibility.'}
        for mutation in case.get('mutations',[]):
            pos = int(mutation[1:-1]); ring = [s[0]['B'][pos][a] for a in ('CG','ND1','CD2','CE1','NE2')]
            r['tag_context']['introduced_histidines'].append({'position':pos, 'ring_to_C_carboxyl_A':distance(ring,terminus),
                'ring_to_N_terminus_A':distance(ring,[binder[0]['N']])})
        # Same ten models as v1: no newly rescued raw candidate enters the direction cohort.
        if (r['id'],r['seed']) in original_direction:
            old = original_direction[(r['id'],r['seed'])]
            d = {'status':'not_scored', 'direction_pass':False, 'original_v1_status':old['status'],
                 'original_v1_shifts':old.get('shifts'), 'original_v1_direction_pass':old['direction_pass']}
            if r['prepared_eligible'] and original_both[r['id']]:
                free = histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',old['positions'])
                bound = histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',old['positions'])
                d.update(status='ok', **direction_metrics(free,bound,terminal_distance(path,old['positions'])))
            else: d['reason'] = 'Raw, preparation, internal geometry, protonation or unchanged structural gate failed; diagnostic tables retained but not admitted.'
            r['direction'] = d
        rows.append(r)
    controls = [r for r in rows if r['kind']!='novel_design']
    # Negative controls are diagnostic; inability to prepare a presumed nonbinder does
    # not negate raw confidence separation, but prevents claiming full-panel validation.
    controls_pass = len(controls)==6 and all(r['prepared_eligible'] for r in controls)
    variants = []
    for old in v1['variants']:
        rs = [r for r in rows if r['id']==old['id']]
        parents = [r for r in rows if r['id']==old['parent_id']]
        variants.append({'id':old['id'], 'parent_id':old['parent_id'],
            'raw_both_seeds_pass':original_both[old['id']],
            'prepared_both_seeds_pass':len(rs)==2 and all(r['prepared_eligible'] for r in rs),
            'matched_parent_prepared_both_seeds_pass':len(parents)==2 and all(r['prepared_eligible'] for r in parents),
            'direction_pass_both_seeds':controls_pass and len(rs)==2 and all(r.get('direction',{}).get('direction_pass',False) for r in rs)})
    result = {'source_preparation_sha256':sha(out/'results.json'), 'source_v1_direction_sha256':sha(ROOT/'reports/expansion_his_direction/results.json'),
        'controls_prepared_pass':controls_pass,'rows':rows,'variants':variants,
        'interpretation':'Post-hoc uniform sensitivity analysis. Repaired diagnostic outcomes do not overwrite original outcomes or constitute independent replications. No experimental affinity, tag compatibility, mouse cross-reactivity or pH switch established.'}
    (out/'analysis.json').write_text(json.dumps(result,indent=2))
    for v in variants:
        print(v)
        for r in rows:
            if r['id']==v['id']: print(r['seed'], r.get('direction'))
    print('Controls', controls_pass, 'eligible', sum(r['prepared_eligible'] for r in rows), '/', len(rows))

if __name__=='__main__': main()
