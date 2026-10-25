"""Freeze the final anchor-preserving panel before any mutant prediction."""
import copy,itertools,json
import numpy as np
from novelty_search import ROOT,sha,save,now
from prepare_expansion_his import site_metrics
from prepare_next_round import mutate_reference
from prepare_controls import read_structure,sequence,save_pdb,save_input_cif
from receptor_preparation_v2 import PROTOCOL
PARENT='diversity_long_037'
COHORT='Every confidence-selected final His panel candidate and control, including failures.'

def main():
    src=ROOT/'data/diversity_confirmation';out=ROOT/'data/final_his_panel'
    assert not out.exists()
    m=json.loads((src/'manifest.json').read_text())
    paths=[ROOT/f'reports/diversity_confirmation/review_bundle/{PARENT}_seed{s}/raw.pdb' for s in [2309,2411]]
    models=[read_structure(p) for p in paths];metrics=[site_metrics(s) for s in models];audit=[]
    for pos in sorted(metrics[0]):
        rows=[a[pos] for a in metrics];reasons=[]
        if any(r['residue'] in 'GPCHDEKR' for r in rows):reasons.append('preserve charged, glycine, proline, cysteine and existing histidine')
        if pos<=5 or pos>90:reasons.append('preserve five residues at each terminus')
        if any(r['salt_contact'] or r['high_contact_anchor_proxy'] for r in rows):reasons.append('anchor proxy')
        if any(not (r['min_target_distance_A']<=5 and 2<=r['buried_sasa_A2']<=65 and r['free_sasa_A2']>=25) for r in rows):reasons.append('not peripheral in both poses')
        audit.append({'position':pos,'models':rows,'excluded_reasons':reasons,'worst_pose_heuristic':min(r['heuristic'] for r in rows)})
    eligible=sorted([r for r in audit if not r['excluded_reasons']],key=lambda r:(-r['worst_pose_heuristic'],r['position']))[:12]
    if not any(all(r['nearby_basic_atoms']>0 for r in a['models']) for a in eligible):raise SystemExit('No reproducible basic-neighbour site; do not launch.')
    pairs=[]
    for a,b in itertools.combinations(eligible,2):
        if not any(all(r['nearby_basic_atoms']>0 for r in c['models']) for c in [a,b]):continue
        ps=sorted([a['position'],b['position']]);ds=[float(np.linalg.norm(s[0]['B'][ps[0]]['CA'].coord-s[0]['B'][ps[1]]['CA'].coord)) for s in models]
        if all(5<=d<=16 for d in ds):pairs.append({'positions':ps,'CA_distances_A':ds,'score':a['worst_pose_heuristic']+b['worst_pose_heuristic']})
    pairs.sort(key=lambda p:(-p['score'],p['positions']))
    variants=[(r['position'],) for r in eligible]+[tuple(p['positions']) for p in pairs[:2]]
    assert 1<=len(variants)<=16
    out.mkdir();cases=[]
    for c in m['cases']:
        if c['id']==PARENT or c['arm']=='control':
            cases.append(copy.deepcopy(c))
            for ext in ['.pdb','.cif','.npz']:(out/(c['id']+ext)).write_bytes((src/(c['id']+ext)).read_bytes())
    original=next(c for c in cases if c['id']==PARENT);s=read_structure(src/f'{PARENT}.pdb');seq=original['binder_sequence']
    for positions in variants:
        mutations=[f'{seq[p-1]}{p}H' for p in positions];name=PARENT+'_'+'_'.join(mutations)
        mutant=mutate_reference(s,positions);newseq=sequence(mutant[0]['B'])
        save_pdb(mutant,out/f'{name}.pdb');save_input_cif(mutant,out/f'{name}.cif',m['target_sequence'],newseq)
        np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(newseq)))
        cases.append({**original,'id':name,'binder_sequence':newseq,'parent_id':PARENT,'mutations':mutations})
    save(out/'contact_audit.json',{'sources':{str(p.relative_to(ROOT)):sha(p) for p in paths},'site_audit':audit,'chosen_sites':[r['position'] for r in eligible],'chosen_pairs':pairs[:2]})
    for key in ['input_sha256','prior_campaign','diversity_policy','budget_allocation','selected_ids','source_results_sha256']:m.pop(key,None)
    protocol=copy.deepcopy(PROTOCOL);protocol['cohort']=COHORT
    plan={'max_variants':16,'site_policy':'Neutral peripheral sites in both parent poses; preserve charged residues, high-contact anchors and five terminal residues. Up to twelve singles and two pairs; pairs include a site with basic neighbours in both poses, CA spacing 5-16 A. Geometry is a hypothesis, not interaction energy.',
      'seeds':[2503,2609],'samples_per_seed':5,'selection':'Unchanged target-integrity-aware confidence selection. Both selected models and all selected controls must pass; no alternative-sample rescue or seed retries.',
      'direction':'Unchanged v2 preparation of ALL selected models. Direction requires raw pass, matched parent pass, preparation movement/internal geometry and prepared screen pass, valid PROPKA groups; each introduced His shift <= -0.5, at least one free pKa 6-7, His-only ratio >=2, artificial target termini >=10 A, consistently in both poses.',
      'tag_followup':'At most two two-pose direction passes, ranked by worst-pose His-only ratio descending then worst-pose ipTM descending then ID. Matched bare/tagged parent and variants, plus all three controls; fresh seeds 2707/2801. No qualifying variant means no tag job.',
      'budget':{'panel_usd':1.5,'conditional_tag_usd':1.5,'total_new_usd':3,'project_cap_usd':30},'stop':'One panel and at most one conditional tag job; freeze submission afterwards. No broad generation.'}
    m.update(cases=cases,selected_ids=[c['id'] for c in cases if c['kind']=='novel_design'],source_run='diversity-confirm-024-6wrw',scope='Final bounded anchor-preserving His panel',prepared_utc=now(),prospective_plan=plan,structural_preparation=protocol,refolding={**m['refolding'],'seeds':plan['seeds']},max_prediction_seconds=1500,experiment_cap_usd=1.5)
    m['input_sha256']={p.name:sha(p) for p in out.iterdir() if p.is_file()};save(out/'manifest.json',m)
    report=ROOT/'reports/final_his_panel';report.mkdir();save(report/'plan.json',m)
    print(json.dumps({'variants':len(variants),'cases':len(cases),'predictions':len(cases)*10,'sites':[r['position'] for r in eligible],'pairs':pairs[:2]},indent=2))
if __name__=='__main__':main()
