"""Prepare a bounded, prospectively specified His panel on three confirmed parents."""
import argparse,copy,hashlib,itertools,json,re
from pathlib import Path
import numpy as np
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
from prepare_controls import ROOT,read_structure,sequence,save_pdb,save_input_cif
from prepare_next_round import mutate_reference,salt_contact,finish
from repair_mutant_input_names import pinned_target_regex

PARENTS=['expansion_apical_long_092','expansion_apical_short_018','expansion_apical_short_026']
COMPARISONS=['apical_short_008','apical_short_008_A41H']
PLAN={'parents':PARENTS,'maximum_variants_per_parent':8,'site_policy':{'max_sites':5,'min_free_sasa_A2':25,'buried_sasa_range_A2':[2,110],'max_contact_distance_A':5,'salt_anchor_distance_A':3.5,'high_contact_anchor':'Exclude residues in the top 20% by buried SASA AND contact atom-pair count among contacting residues in either source pose. Geometric proxy, not energy.','strata':'Prefer up to three peripheral sites (buried SASA <=65 in both poses) and two partially buried sites (>65 in either pose); fill vacancies by worst-pose 3*basic_atoms-oxygen_atoms.','pair_policy':'At most three pairs, CA separation 5-16 A in both poses; rank summed worst-pose site heuristic, lexical tie break.'},'screen_seed':1709,'confirmation_seeds':[1801,1901],'samples_per_seed':5,'confirmation':'Require all selected structural controls to pass; at most two new passing variants per parent by selected interface ipTM descending, lexical tie break. Keep all three parents and both A41H comparison sequences. No substitute samples or seed retries.','direction':'Both selected fresh-seed mutant models and controls must pass raw checks; same frozen preparation and strict direction criteria as prior round. Report matched-parent eligibility separately.','budget':{'screen_usd':2,'confirmation_usd':2,'total_new_reservation_usd':4,'project_gross_cap_usd':30}}


def site_metrics(s):
    bound=copy.deepcopy(s[0]);free=copy.deepcopy(s[0]);free.detach_child('A')
    ShrakeRupley(n_points=100).compute(bound,level='R');ShrakeRupley(n_points=100).compute(free,level='R')
    target=[a for a in bound['A'].get_atoms() if a.element!='H'];xyz=np.array([a.coord for a in target]);rows=[]
    for r in bound['B']:
        own=[a for a in r if a.element!='H'];dist=np.linalg.norm(np.array([a.coord for a in own])[:,None]-xyz[None],axis=2)
        nearby=[target[i] for i in np.where(dist.min(axis=0)<5)[0]]
        basic=sum(a.name in {'NZ','NH1','NH2','NE'} and a.parent.resname in {'LYS','ARG'} for a in nearby)
        oxy=sum(a.element=='O' for a in nearby)
        rows.append({'position':r.id[1],'residue':seq1(r.resname),'free_sasa_A2':float(free['B'][r.id].sasa),'buried_sasa_A2':float(free['B'][r.id].sasa-r.sasa),'min_target_distance_A':float(dist.min()),'contact_atom_pairs_lt5A':int((dist<5).sum()),'nearby_basic_atoms':basic,'nearby_oxygen_atoms':oxy,'heuristic':3*basic-oxy,'salt_contact':bool(salt_contact(s,r.id[1]))})
    contacting=[r for r in rows if r['min_target_distance_A']<=5 and r['buried_sasa_A2']>=2]
    n=max(1,int(np.ceil(.2*len(contacting))))
    bsa={r['position'] for r in sorted(contacting,key=lambda r:(-r['buried_sasa_A2'],r['position']))[:n]}
    counts={r['position'] for r in sorted(contacting,key=lambda r:(-r['contact_atom_pairs_lt5A'],r['position']))[:n]}
    for r in rows:r['high_contact_anchor_proxy']=r['position'] in bsa&counts
    return {r['position']:r for r in rows}


def select(parent,selected):
    models=[read_structure(ROOT/r['raw_scoring_output']/'model.pdb') for r in selected];metrics=[site_metrics(s) for s in models];audit=[]
    for pos in sorted(metrics[0]):
        ms=[m[pos] for m in metrics];reasons=[]
        if any(r['residue'] in 'GPC H'.replace(' ','') for r in ms):reasons.append('excluded residue identity')
        if any(r['salt_contact'] for r in ms):reasons.append('possible charged anchor')
        if any(r['high_contact_anchor_proxy'] for r in ms):reasons.append('high-contact anchor proxy')
        if any(not (r['min_target_distance_A']<=5 and 2<=r['buried_sasa_A2']<=110 and r['free_sasa_A2']>=25) for r in ms):reasons.append('outside site geometry in a source pose')
        audit.append({'position':pos,'models':ms,'excluded_reasons':reasons,'stratum':'peripheral' if all(r['buried_sasa_A2']<=65 for r in ms) else 'partially_buried','worst_pose_heuristic':min(r['heuristic'] for r in ms)})
    eligible=sorted([r for r in audit if not r['excluded_reasons']],key=lambda r:(-r['worst_pose_heuristic'],r['position']))
    chosen=[r for r in eligible if r['stratum']=='peripheral'][:3]+[r for r in eligible if r['stratum']=='partially_buried'][:2]
    chosen+=( [r for r in eligible if r not in chosen][:5-len(chosen)] )
    pairs=[]
    for a,b in itertools.combinations(sorted(chosen,key=lambda r:r['position']),2):
        ps=[a['position'],b['position']];ds=[float(np.linalg.norm(s[0]['B'][ps[0]]['CA'].coord-s[0]['B'][ps[1]]['CA'].coord)) for s in models]
        if all(5<=d<=16 for d in ds):pairs.append({'positions':ps,'CA_distances_A':ds,'score':a['worst_pose_heuristic']+b['worst_pose_heuristic']})
    pairs.sort(key=lambda p:(-p['score'],p['positions']))
    variants=[(r['position'],) for r in chosen]+[tuple(r['positions']) for r in pairs[:3]]
    return variants,{'parent':parent,'sources':[{'seed':r['seed'],'sample':r['sample'],'source':r['raw_scoring_output']+'/model.pdb','sha256':hashlib.sha256((ROOT/r['raw_scoring_output']/'model.pdb').read_bytes()).hexdigest()} for r in selected],'site_audit':audit,'chosen_sites':[r['position'] for r in chosen],'chosen_pairs':pairs[:3]}


def choose_confirmation(result,m):
    if not result['controls_pass']:return None
    lookup={c['id']:c for c in m['cases']};chosen=[]
    for parent in PARENTS:
        eligible=[r for r in result['candidate_summaries'] if lookup[r['id']].get('parent_id')==parent and r['all_selected_seeds_candidate_screen']]
        eligible.sort(key=lambda r:(-min(x['metrics']['design_to_target_iptm'] for x in r['selected_models']),r['id']))
        chosen.extend(r['id'] for r in eligible[:2])
    return chosen


def main():
    p=argparse.ArgumentParser();p.add_argument('--confirm',action='store_true');args=p.parse_args()
    if args.confirm:
        src=ROOT/'data/expansion_his_panel';m=json.loads((src/'manifest.json').read_text());rpath=ROOT/'reports/expansion_his_panel/results.json';r=json.loads(rpath.read_text());chosen=choose_confirmation(r,m)
        if not chosen:raise SystemExit('No new eligible variants or controls fail: no automatic confirmation.')
        out=ROOT/'data/expansion_his_confirmation';out.mkdir(exist_ok=False)
        ids=chosen+PARENTS+COMPARISONS;m['cases']=[c for c in m['cases'] if c['id'] in ids or c['arm']=='control']
        for c in m['cases']:
            for ext in ('.pdb','.cif','.npz'):
                name=c['id']+ext;assert hashlib.sha256((src/name).read_bytes()).hexdigest()==m['input_sha256'][name];(out/name).write_bytes((src/name).read_bytes())
        m.update(source_run=r['run'],source_results_sha256=hashlib.sha256(rpath.read_bytes()).hexdigest(),selected_ids=ids,selection_details={'chosen_variants':chosen,'max_per_parent':2},refolding={**m['refolding'],'seeds':PLAN['confirmation_seeds']},scope='Fresh-seed confirmation of selected new His variants; A41H comparison retained')
        finish(out,m);print(json.dumps({'chosen':chosen,'cases':len(m['cases'])},indent=2));return
    src=ROOT/'data/expansion_validation';m=json.loads((src/'manifest.json').read_text());result=json.loads((ROOT/'reports/expansion_validation/results.json').read_text())
    assert result['controls_pass'];out=ROOT/'data/expansion_his_panel';out.mkdir(exist_ok=False);cases=[];audits=[]
    for c in m['cases']:
        if c['id'] in PARENTS or c['arm']=='control':
            cases.append(copy.deepcopy(c))
            for ext in ('.pdb','.cif','.npz'):(out/(c['id']+ext)).write_bytes((src/(c['id']+ext)).read_bytes())
    for parent in PARENTS:
        r=next(c for c in result['candidate_summaries'] if c['id']==parent);assert r['all_selected_seeds_full_screen'];variants,audit=select(parent,r['selected_models']);audits.append(audit)
        original=next(c for c in m['cases'] if c['id']==parent);s=read_structure(src/f'{parent}.pdb');seq=original['binder_sequence']
        for ps in variants:
            mutations=[f'{seq[i-1]}{i}H' for i in ps];name=parent+'_'+'_'.join(mutations);assert re.fullmatch(pinned_target_regex(),name)
            mutant=mutate_reference(s,ps);newseq=sequence(mutant[0]['B']);save_pdb(mutant,out/f'{name}.pdb');save_input_cif(mutant,out/f'{name}.cif',m['target_sequence'],newseq)
            np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(newseq)))
            cases.append({**original,'id':name,'binder_sequence':newseq,'parent_id':parent,'mutations':mutations})
    oldsrc=ROOT/'data/his_confirmation';old=json.loads((oldsrc/'manifest.json').read_text())
    m['arms'] += [a for a in old['arms'] if a['id']=='apical_short']
    for c in old['cases']:
        if c['id'] not in COMPARISONS:continue
        cases.append({**c,'comparison_only':True})
        for ext in ('.pdb','.cif','.npz'):(out/(c['id']+ext)).write_bytes((oldsrc/(c['id']+ext)).read_bytes())
    for k in ['input_sha256','prepared_utc','repair','prior_campaign','diversity_policy','selection_details','budget_allocation','source_results_sha256','controls']:m.pop(k,None)
    m.update(cases=cases,source_run=result['run'],scope='Bounded His variants on three confirmed apical parents; A41H comparison',selected_ids=[c['id'] for c in cases if c['kind']=='novel_design'],prospective_plan=PLAN,refolding={**m['refolding'],'seeds':[1709],'diffusion_samples_per_seed':5},experiment_cap_usd=2,max_prediction_seconds=2100,total_cap_usd=30)
    (out/'contact_audit.json').write_text(json.dumps(audits,indent=2));finish(out,m)
    report=ROOT/'reports/expansion_his_panel';direction=json.loads((ROOT/'reports/his_panel/direction_plan.json').read_text());direction['matched_parent_policy']='Report parent pass in both fresh seeds separately; variant-only passes are not evidence of improved affinity.';(report/'direction_plan.json').write_text(json.dumps(direction,indent=2))
    print(json.dumps({'cases':len(cases),'new_variants':sum(c.get('parent_id') in PARENTS for c in cases),'predictions':len(cases)*5,'sites':[{k:a[k] for k in ['parent','chosen_sites','chosen_pairs']} for a in audits]},indent=2))

if __name__=='__main__':main()
