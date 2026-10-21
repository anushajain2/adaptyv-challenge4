"""Prepare the approved bounded diagnostic panel and apical expansion."""
import copy
import hashlib
import itertools
import json
import re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import yaml
from direction_check import enumerate_pairs
from prepare_controls import ROOT, read_structure, sequence, save_pdb, save_input_cif
from structural_preparation import PROTOCOL
from repair_mutant_input_names import pinned_target_regex


def mutate_reference(parent, positions):
    s=copy.deepcopy(parent)
    assert len(set(positions))==len(positions)
    for p in positions:
        r=s[0]['B'][p]
        assert r.resname not in {'HIS','GLY','PRO','CYS'}
        r.resname='HIS'
        for a in list(r):
            if a.name not in {'N','CA','C','O'}:r.detach_child(a.id)
    return s


def salt_contact(s, pos):
    r=s[0]['B'][pos]
    if r.resname in {'ASP','GLU'}:
        own=[a for a in r if a.name in {'OD1','OD2','OE1','OE2'}]
        other=[a for a in s[0]['A'].get_atoms() if (a.parent.resname=='ARG' and a.name in {'NH1','NH2','NE'}) or (a.parent.resname=='LYS' and a.name=='NZ')]
    elif r.resname in {'ARG','LYS'}:
        own=[a for a in r if a.name in {'NH1','NH2','NE','NZ'}]
        other=[a for a in s[0]['A'].get_atoms() if a.parent.resname in {'ASP','GLU'} and a.name in {'OD1','OD2','OE1','OE2'}]
    else:return False
    return any(np.linalg.norm(a.coord-b.coord)<3.5 for a in own for b in other)


def select_sites(parent_id, rows):
    models=[read_structure(ROOT/r['raw_scoring_output']/'model.pdb') for r in rows]
    options=[enumerate_pairs(s,10000) for s in models]
    site_maps=[{s['position']:s for p in ps for s in p['sites']} for ps in options]
    common=set.intersection(*(set(m) for m in site_maps))
    audit=[]
    for p in sorted(set.union(*(set(m) for m in site_maps))):
        salt=any(salt_contact(s,p) for s in models)
        reasons=[]
        if p not in common:reasons.append('not peripheral contact in both source poses')
        if salt:reasons.append('possible charged anchor within 3.5 A in a source pose')
        if parent_id=='apical_long_015' and p==51:reasons.append('preserve Asp51 except diagnostic')
        scores=[m[p]['heuristic'] for m in site_maps if p in m]
        audit.append({'position':p,'sites':[m.get(p) for m in site_maps], 'excluded_reasons':reasons,'worst_pose_heuristic':min(scores)})
    eligible=sorted((r for r in audit if not r['excluded_reasons']),key=lambda r:(-r['worst_pose_heuristic'],r['position']))
    chosen=[r['position'] for r in eligible[:5]]
    assert len(chosen)==5, (parent_id,chosen)
    pairs=[]
    for a,b in itertools.combinations(sorted(chosen),2):
        ds=[float(np.linalg.norm(s[0]['B'][a]['CA'].coord-s[0]['B'][b]['CA'].coord)) for s in models]
        if all(5<=d<=16 for d in ds):
            score=sum(next(r['worst_pose_heuristic'] for r in eligible if r['position']==p) for p in (a,b))
            pairs.append({'positions':[a,b],'CA_distances_A':ds,'score':score})
    pairs.sort(key=lambda p:(-p['score'],p['positions']))
    return [(p,) for p in chosen]+[tuple(p['positions']) for p in pairs[:4]], {'parent':parent_id,'sources':[{'seed':r['seed'],'sample':r['sample'],'sha256':hashlib.sha256((ROOT/r['raw_scoring_output']/'model.pdb').read_bytes()).hexdigest()} for r in rows], 'site_audit':audit,'chosen_sites':chosen,'chosen_pairs':pairs[:4]}


def finish(out,m):
    m['prepared_utc']=datetime.now(timezone.utc).isoformat()
    m['structural_preparation']=PROTOCOL
    m['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(m,indent=2))
    report=ROOT/'reports'/out.name; report.mkdir(exist_ok=False)
    (report/'plan.json').write_text(json.dumps(m,indent=2))


def panel():
    src=ROOT/'data/confirmation';out=ROOT/'data/his_panel';out.mkdir(exist_ok=False)
    m=json.loads((src/'manifest.json').read_text())
    originals=m['cases'];cases=[];audits=[]
    parents=['apical_long_015','apical_short_008']
    for c in originals:
        if c['id'] in parents or c['arm']=='control':
            cases.append(copy.deepcopy(c))
            for ext in ('.cif','.npz','.pdb'):(out/(c['id']+ext)).write_bytes((src/(c['id']+ext)).read_bytes())
    for parent in parents:
        report='mutant_validation' if parent=='apical_long_015' else 'confirmation'
        result=json.loads((ROOT/f'reports/{report}/results.json').read_text())
        rows=[r for r in result['rows'] if r['id']==parent and r['confidence_selected']]
        variants,audit=select_sites(parent,rows);audits.append(audit)
        if parent=='apical_long_015':variants += [(51,),(70,),(51,70)]
        variants=list(dict.fromkeys(tuple(sorted(p)) for p in variants))
        original=next(c for c in originals if c['id']==parent)
        s=read_structure(src/f'{parent}.pdb');seq=original['binder_sequence']
        for ps in variants:
            mutations=[f'{seq[p-1]}{p}H' for p in ps];name=parent+'_'+'_'.join(mutations)
            assert re.fullmatch(pinned_target_regex(),name)
            mutant=mutate_reference(s,ps);newseq=sequence(mutant[0]['B'])
            save_pdb(mutant,out/f'{name}.pdb');save_input_cif(mutant,out/f'{name}.cif',m['target_sequence'],newseq)
            np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(newseq)))
            cases.append({**original,'id':name,'binder_sequence':newseq,'parent_id':parent,'mutations':mutations,
                          'diagnostic_only':ps in [(51,),(70,),(51,70)] and parent=='apical_long_015'})
    (out/'contact_audit.json').write_text(json.dumps(audits,indent=2))
    m.update(cases=cases,selected_ids=[c['id'] for c in cases if c['kind']=='novel_design'],
             scope='Bounded His diagnostic panel including both parents; no experimental claims',
             refolding={**m['refolding'],'seeds':[1009],'diffusion_samples_per_seed':5},
             max_prediction_seconds=1500,experiment_cap_usd=1.5,total_cap_usd=30,
             next_round_allocation={'panel_plus_confirmation_usd':3,'apical_expansion_usd':6},
             prospective_decision={'raw_screen':'Unchanged structural thresholds; target-integrity-aware confidence selection.',
                'preparation':'Apply once to every selected candidate/control; report separately without replacing raw calls.',
                'confirmation':'If selected controls pass, confirm at most six raw-screen-passing new variants, both parents and all three controls. Rank variants by interface ipTM with at most three per parent. Exclude failed double-mutant comparator. Seeds 1103 and 1201, five samples each. Parents retained even if their first-stage screen fails. No alternate seeds.',
                'direction':'Only variants passing both fresh selected seeds with controls passing receive PROPKA direction triage on independently predicted interfaces. No raw-model gate relaxation.'})
    finish(out,m)
    print(json.dumps({'panel_cases':len(cases),'predictions':len(cases)*5,'audits':[{k:a[k] for k in ('parent','chosen_sites','chosen_pairs')} for a in audits]},indent=2))


def expansion():
    src=ROOT/'data/broad_campaign';out=ROOT/'data/apical_expansion';out.mkdir(exist_ok=False)
    m=json.loads((src/'manifest.json').read_text())
    arms=[]
    for i,old in enumerate(a for a in m['arms'] if a['id'] in {'apical_short','apical_long'}):
        a=copy.deepcopy(old);a.update(id='expansion_'+a['id'],num_designs=100,generation_seed=20261201+i,inverse_fold_seed=20261301+i)
        arms.append(a);(out/f"{a['id']}.yaml").write_bytes((src/f"{old['id']}.yaml").read_bytes())
    for name in ('target_6wrw.cif','spatial_map.csv'):(out/name).write_bytes((src/name).read_bytes())
    for c in m['controls']:
        for ext in ('.cif','.npz','.pdb'):(out/(c['id']+ext)).write_bytes((src/(c['id']+ext)).read_bytes())
    m.update(arms=arms,scope='200 new apical miniproteins, 100 per length range; one-sample triage only',
             refolding={**m['refolding'],'seeds':[1301]},
             selection_rule=json.loads((ROOT/'data/confirmation/manifest.json').read_text())['selection_rule'],
             prior_campaign='reports/broad_campaign/results.json',
             diversity_policy='Exclude exact prior sequences; shortlist by existing 80% sequence-identity rule across both campaigns. Report structural diversity after deep confirmation.',
             budget_allocation={'this_expansion_usd':6,'separate_panel_plus_confirmation_usd':3})
    finish(out,m)
    print('Prepared 200-design apical expansion')

if __name__=='__main__':
    panel();expansion()
