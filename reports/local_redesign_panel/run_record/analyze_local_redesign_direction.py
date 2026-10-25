"""Prospective local-redesign comparison; strict and sensitivity outcomes stay separate."""
import csv
import json
import numpy as np
from novelty_search import ROOT,sha,save
from prepare_controls import read_structure
from analyze_controls import pose_metrics
from analyze_generation_comparison import design_pose_metrics,contacts_summary,screen_sample
from analyze_receptor_preparation_v2 import direction_metrics
from direction_check import histidines
from summarize_direction import terminal_distance
from audit_preparation_motion import motion_audit


def his_environment(structure, position):
    his=structure[0]['B'][position]; nitrogens=[his[n] for n in ['ND1','NE2']];rows=[]
    for chain in structure[0]:
        for residue in chain:
            if chain.id=='B' and residue.id[1]==position:continue
            for atom in residue:
                if atom.element not in {'N','O'}:continue
                distances=[float(np.linalg.norm(n.coord-atom.coord)) for n in nitrogens]
                if min(distances)>6:continue
                category='polar_atom'
                if (residue.resname=='ARG' and atom.name in {'NE','NH1','NH2'}) or (residue.resname=='LYS' and atom.name=='NZ'):category='basic_cationic_group'
                if (residue.resname=='ASP' and atom.name in {'OD1','OD2'}) or (residue.resname=='GLU' and atom.name in {'OE1','OE2'}):category='acidic_carboxyl_group'
                rows.append({'chain':chain.id,'position':residue.id[1],'native_target_position':residue.id[1]+189 if chain.id=='A' else None,
                             'residue':residue.resname,'atom':atom.name,'category':category,'distance_to_ND1_NE2_A':distances})
    return {'position':position,'neighbours':sorted(rows,key=lambda r:min(r['distance_to_ND1_NE2_A'])),
            'interpretation':'Distance inventory only. Nitrogen proximity does not establish an acceptor hydrogen bond, protonation microstate, repulsion energy or pH-dependent affinity.'}


def main():
    folder=ROOT/'reports/local_redesign_panel';rawpath=folder/'results.json';preppath=folder/'preparation/results.json'
    raw=json.loads(rawpath.read_text());prep=json.loads(preppath.read_text());m=json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text())
    assert prep['source_results_sha256']==sha(rawpath)
    cases={c['id']:c for c in m['cases']};arms={a['id']:a for a in m['arms']}
    selected={(r['id'],r['seed'],r['sample']):r for r in raw['rows'] if r['confidence_selected']}
    assert set(selected)=={(r['id'],r['seed'],r['sample']) for r in prep['rows']}
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    rows=[]
    for p in sorted(prep['rows'],key=lambda r:(r['id'],r['seed'])):
        previous=selected[p['id'],p['seed'],p['sample']];c=cases[p['id']]
        row={'id':p['id'],'seed':p['seed'],'sample':p['sample'],'kind':c['kind'],'strict_eligible':False,'sensitivity_eligible':False,
             'raw_pass':previous.get('screen',{}).get('pass',previous.get('control_pass')),'direction':{'status':'excluded','direction_pass':False}}
        try:
            assert p['status']=='ok','Preparation execution failed'
            path=ROOT/p['dest']/'prepared.pdb';assert sha(path)==p['prepared_sha256'] and sha(ROOT/p['source'])==p['source_sha256']
            structure=read_structure(path);reference=read_structure(ROOT/'data/local_redesign_panel'/f"{c['id']}.pdb")
            if c['kind']=='novel_design':
                check=screen_sample({'metrics':previous['metrics'],'pose':design_pose_metrics(structure,reference),
                    'geometry':p['prepared_geometry'],'contacts':contacts_summary(structure,arms[c['arm']],conserved)},m['screen'])
            elif c['kind']=='positive':
                pose=pose_metrics(structure,reference)
                check={'pass':previous['metrics']['design_to_target_iptm']>=.5 and pose['binder_CA_fold_rmsd_A']<=2 and
                    pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and pose['binder_interface_CA_rmsd_after_target_fit_A']<=2 and pose['native_contact_recall']>=.5,'pose':pose}
            else:check={'pass':previous['control_pass'],'note':'Original confidence-separation check retained; preparation cannot establish nonbinding.'}
            motion=motion_audit(read_structure(ROOT/p['source']),structure)
            common=bool(row['raw_pass'] and check['pass'] and p['prepared_internal_geometry']['pass'] and all(v['status']=='ok' for v in p['propka'].values()))
            row.update(prepared_screen=check,strict_eligible=common and p['displacement']['pass'],
                sensitivity_eligible=common and motion['interface_sensitivity_pass'],motion=motion,displacement=p['displacement'],
                internal_geometry=p['prepared_internal_geometry'],propka=p['propka'],source=str(path.relative_to(ROOT)),source_sha256=sha(path))
            positions=c.get('histidine_positions',[])
            if positions:
                row['histidine_environments']=[his_environment(structure,pos) for pos in positions]
                if all(v['status']=='ok' for v in p['propka'].values()):
                    free=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',positions)
                    bound=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',positions)
                    row['direction']={'status':'diagnostic',**direction_metrics(free,bound,terminal_distance(path,positions))}
        except Exception as e:row['error']=f'{type(e).__name__}: {e}'
        rows.append(row)
    byid={name:sorted([r for r in rows if r['id']==name],key=lambda r:r['seed']) for name in cases}
    controls=[r for r in rows if r['kind']!='novel_design'];strict_controls=len(controls)==6 and all(r['strict_eligible'] for r in controls)
    sensitivity_controls=len(controls)==6 and all(r['sensitivity_eligible'] for r in controls)
    variants=[]
    for c in m['cases']:
        if not c.get('local_redesign'):continue
        rs=byid[c['id']];base=byid[c['parent_id']];original=byid[c['unmodified_parent_id']]
        strict=len(rs)==2 and all(r['strict_eligible'] and r['direction']['direction_pass'] for r in rs)
        sensitive=len(rs)==2 and all(r['sensitivity_eligible'] and r['direction']['direction_pass'] for r in rs)
        baseline_ok=len(base)==2 and all(r['strict_eligible'] and r['direction']['status']=='diagnostic' for r in base)
        parent_ok=len(original)==2 and all(r['strict_eligible'] for r in original)
        ratios=[]
        if len(rs)==2 and baseline_ok and all(r['direction']['status']=='diagnostic' for r in rs):
            assert [r['seed'] for r in rs]==[r['seed'] for r in base]
            ratios=[r['direction']['his_only_Kd6_over_Kd7_4']/b['direction']['his_only_Kd6_over_Kd7_4'] for r,b in zip(rs,base)]
        improved=baseline_ok and bool(ratios) and (all(x>=1.2 for x in ratios) or not all(r['direction']['direction_pass'] for r in base))
        variants.append({'id':c['id'],'sequence':c['binder_sequence'],'baseline_id':c['parent_id'],'mutations':c['mutations'],
            'strict_structure_and_direction_both':strict,'sensitivity_structure_and_direction_both':sensitive,
            'baseline_strict_eligible_both':baseline_ok,'unmodified_parent_strict_eligible_both':parent_ok,
            'paired_His_only_ratio_relative_to_baseline':ratios,'prospective_improvement_heuristic':bool(improved),
            'fully_qualified_for_followup':bool(raw['controls_pass'] and strict_controls and strict and baseline_ok and parent_ok and improved),
            'models':rs})
    eligible=[v for v in variants if v['fully_qualified_for_followup']]
    eligible.sort(key=lambda v:(-min(r['direction']['his_only_Kd6_over_Kd7_4'] for r in v['models']),
        -min(selected[r['id'],r['seed'],r['sample']]['metrics']['design_to_target_iptm'] for r in v['models']),v['id']))
    result={'raw_results_sha256':sha(rawpath),'preparation_results_sha256':sha(preppath),'raw_controls_pass':raw['controls_pass'],
        'strict_prepared_controls_pass':strict_controls,'sensitivity_prepared_controls_pass':sensitivity_controls,
        'rows':rows,'variants':variants,'followup_ids':[v['id'] for v in eligible[:2]],
        'limitations':['Fixed-conformation PROPKA comparison; no unbound ensemble','His-only independent-site ratio is a heuristic, not measured affinity or full proton linkage',
            'Same predictor fresh seeds measure repeatability, not independent experimental validation','No experimental pH, binding, mouse or tag claim','No historical result or fallback submission changed']}
    save(folder/'direction.json',result)
    print(json.dumps({'raw_controls':raw['controls_pass'],'strict_controls':strict_controls,'sensitivity_controls':sensitivity_controls,
                     'strict_direction_variants':sum(v['strict_structure_and_direction_both'] for v in variants),'followup_ids':result['followup_ids']},indent=2))


if __name__=='__main__':main()
