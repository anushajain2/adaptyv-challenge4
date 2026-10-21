"""PROPKA on independently predicted, confirmed mutant interfaces; no new mutations."""
import copy
import csv
import hashlib
import json
from pathlib import Path
from prepare_controls import ROOT, read_structure
from prepare_target import ProteinSelect, save_pdb
from direction_check import pka_table,histidines,release_ratio,HIS_ATOMS
from summarize_direction import terminal_distance
from analyze_generation_comparison import design_pose_metrics,contacts_summary,screen_sample


def main():
    source=ROOT/'reports/expansion_his_confirmation/results.json'
    result=json.loads(source.read_text())
    assert result['controls_pass']
    m=json.loads((ROOT/'data/expansion_his_confirmation/manifest.json').read_text())
    prep=json.loads((ROOT/'reports/expansion_his_confirmation/preparation/results.json').read_text())
    assert len(prep['rows'])==len(m['cases'])*2
    assert prep['source_results_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
    lookup={c['id']:c for c in m['cases']};arms={a['id']:a for a in m['arms']}
    prepared={(r['id'],r['seed']):r for r in prep['rows']}
    out=ROOT/'reports/expansion_his_direction';out.mkdir(exist_ok=False)
    plan=json.loads((ROOT/'reports/expansion_his_panel/direction_plan.json').read_text())
    (out/'plan.json').write_text(json.dumps(plan,indent=2))
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    variants=[]
    for c in result['candidate_summaries']:
        case=lookup[c['id']]
        if not case.get('mutations') or not c['all_selected_seeds_full_screen']:continue
        positions=[int(s[1:-1]) for s in case['mutations']]
        rows=[]
        for prediction in c['selected_models']:
            seed=prediction['seed'];p=prepared[(c['id'],seed)]
            r={'id':c['id'],'seed':seed,'positions':positions,'status':'not_scored','direction_pass':False}
            try:
                assert p['status']=='ok' and p['displacement']['pass'], 'Prepared displacement or execution guard failed'
                path=ROOT/p['dest']/'prepared.pdb';s=read_structure(path)
                ref=read_structure(ROOT/'data/expansion_his_confirmation'/f"{c['id']}.pdb")
                scored={'metrics':prediction['metrics'],'pose':design_pose_metrics(s,ref),'geometry':p['prepared_geometry'],
                        'contacts':contacts_summary(s,arms[case['arm']],conserved)}
                check=screen_sample(scored,m['screen']);r['prepared_screen']=check
                assert check['pass'], 'Prepared model fails unchanged pose/contact/context criteria'
                assert all(HIS_ATOMS<=set(s[0]['B'][pos].child_dict) for pos in positions)
                folder=out/'models'/f"{c['id']}_seed{seed}";folder.mkdir(parents=True)
                (folder/'complex.pdb').write_bytes(path.read_bytes())
                save_pdb(s,folder/'free_binder.pdb',ProteinSelect({'B'}))
                free=histidines(pka_table(folder/'free_binder.pdb'),'B',positions)
                bound=histidines(pka_table(folder/'complex.pdb'),'B',positions)
                shifts=[b-f for f,b in zip(free,bound)]
                distance=terminal_distance(folder/'complex.pdb',positions)
                ratio=release_ratio(free,bound)
                r.update(status='ok',free_pkas=free,bound_pkas=bound,shifts=shifts,
                         artificial_terminus_distance_A=distance,his_only_Kd6_over_Kd7_4=ratio,
                         direction_pass=all(d<=-.5 for d in shifts) and any(6<=f<=7 for f in free) and ratio>=2 and distance>=10,
                         source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            except Exception as e:r.update(status='failed',error=f'{type(e).__name__}: {e}')
            rows.append(r)
            print(c['id'],seed,r['status'],r.get('shifts'),flush=True)
        variants.append({'id':c['id'],'parent_id':case['parent_id'],'mutations':case['mutations'],'sequence':case['binder_sequence'],
                         'comparison_only':case.get('comparison_only',False),'matched_parent_passes_both_seeds':next((x['all_selected_seeds_full_screen'] for x in result['candidate_summaries'] if x['id']==case['parent_id']),False),
                         'models':rows,'direction_pass_both_seeds':len(rows)==2 and all(r['direction_pass'] for r in rows)})
    record={'source_results_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'plan':plan,'variants':variants,
            'limitations':['Fixed conformation in free and bound calculation','194-residue target crop with artificial termini','No full proton linkage or unbound ensemble','Introduced-His independent-site proxy, not measured affinity','No experimental human or mouse binding established']}
    (out/'results.json').write_text(json.dumps(record,indent=2))
    (out/'direction_hypotheses.fasta').write_text(''.join(f">{v['id']} COMPUTATIONAL_HYPOTHESIS_ONLY\n{v['sequence']}\n" for v in variants if v['direction_pass_both_seeds']))
    print(json.dumps({'evaluated':len(variants),'direction_hypotheses':sum(v['direction_pass_both_seeds'] for v in variants)},indent=2))
if __name__=='__main__':main()
