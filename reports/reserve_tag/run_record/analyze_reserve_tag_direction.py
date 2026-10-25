"""Continuous held-out readouts alongside unchanged legacy qualification."""
import csv
import json
import numpy as np
from novelty_search import ROOT,sha,save
from prepare_controls import read_structure,geometry
from analyze_tag_comparison import core_structure,core_screen
from analyze_controls import pose_metrics,contact_pairs
from analyze_generation_comparison import design_pose_metrics,contacts_summary
from audit_preparation_motion import motion_audit
from analyze_receptor_preparation_v2 import direction_metrics
from direction_check import histidines
from summarize_direction import terminal_distance
from audit_local_redesign_mechanism import group_terms
from analyze_local_redesign_direction import his_environment
from prepare_reserve_tag import RESERVE,ABLATION


def main():
    out=ROOT/'reports/reserve_tag';rawpath=out/'results.json';preppath=out/'preparation/results.json'
    raw=json.loads(rawpath.read_text());prep=json.loads(preppath.read_text());m=json.loads((ROOT/'data/reserve_tag/manifest.json').read_text())
    assert prep['source_results_sha256']==sha(rawpath)
    cases={c['id']:c for c in m['cases']};arms={a['id']:a for a in m['arms']}
    selected={(r['id'],r['seed'],r['sample']):r for r in raw['rows'] if r['confidence_selected']}
    assert set(selected)=={(r['id'],r['seed'],r['sample']) for r in prep['rows']}
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    full=read_structure(ROOT/'data/raw/11QC.cif');rows=[]
    for pr in sorted(prep['rows'],key=lambda p:(p['id'],p['seed'])):
        previous=selected[pr['id'],pr['seed'],pr['sample']];c=cases[pr['id']]
        row={k:previous[k] for k in ['id','base_id','condition','kind','seed','sample','core_length']}
        row.update(raw_pass=previous.get('screen',{}).get('pass',previous.get('control_pass')),strict_eligible=False,sensitivity_eligible=False,
            direction={'status':'excluded','direction_pass':False},raw_tag_clear=previous['tag']['pass'])
        try:
            assert pr['status']=='ok','Preparation failed'
            path=ROOT/pr['dest']/'prepared.pdb';assert sha(path)==pr['prepared_sha256'] and sha(ROOT/pr['source'])==pr['source_sha256']
            s=read_structure(path);core=core_structure(s,c['core_length']);ref=read_structure(ROOT/c['reference'])
            if c['kind']=='novel_design':
                check=core_screen({**previous,'pose':design_pose_metrics(core,ref),'geometry':geometry(core,{i:i+189 for i in range(1,195)},full),
                    'contacts':contacts_summary(core,arms[c['arm']],conserved)},m)
            elif c['kind']=='positive':
                pose=pose_metrics(core,ref);check={'pass':previous['core_mean_cross_PAE_A']<m['core_confidence']['max_mean_cross_PAE_A'] and
                    pose['binder_CA_fold_rmsd_A']<=2 and pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and
                    pose['binder_interface_CA_rmsd_after_target_fit_A']<=2 and pose['native_contact_recall']>=.5,'pose':pose}
            else:check={'pass':previous['control_pass'],'note':'Raw negative confidence check retained; minimization cannot recompute confidence.'}
            rawcore=core_structure(read_structure(ROOT/pr['source']),c['core_length']);motion=motion_audit(rawcore,core)
            valid_propka=all(v['status']=='ok' for v in pr['propka'].values())
            common=bool(row['raw_pass'] and check['pass'] and pr['prepared_internal_geometry']['pass'] and valid_propka)
            tags=[list(pair) for pair in sorted(contact_pairs(s)) if pair[1]>c['core_length']]
            g=pr['prepared_geometry'];prepared_tag_clear=not tags and all(g[k]==0 for k in g if 'clash_pairs_lt2A' in k)
            row.update(strict_eligible=common and pr['displacement']['pass'],sensitivity_eligible=common and motion['interface_sensitivity_pass'],
                prepared_screen=check,displacement=pr['displacement'],full_displacement=pr['full_displacement'],motion=motion,
                internal_geometry=pr['prepared_internal_geometry'],propka=pr['propka'],prepared_tag_contacts=tags,
                prepared_full_geometry=g,prepared_tag_clear=prepared_tag_clear,source=str(path.relative_to(ROOT)),source_sha256=sha(path))
            positions=c['histidine_positions']
            if positions and valid_propka:
                assert positions==[41]
                free=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',positions)
                bound=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',positions)
                row['direction']={'status':'diagnostic',**direction_metrics(free,bound,terminal_distance(path,positions))}
                f=group_terms(path.with_name('free_binder.pdb'),41,free[0]);b=group_terms(path,41,bound[0])
                row['mechanism']={'free':f,'bound':b,'delta_pka_contributions':{key:b['contributions_pka_units'][key]-value for key,value in f['contributions_pka_units'].items()},'environment':his_environment(s,41)}
        except Exception as e:row['error']=f'{type(e).__name__}: {e}'
        rows.append(row)
    controls={}
    for condition in ['bare','tagged']:
        rs=[r for r in rows if r['condition']==condition and r['kind']!='novel_design'];positive=[r for r in rs if r['kind']=='positive']
        controls[condition]={'raw_pass':raw['controls'][condition]['pass'],
            'strict_full_panel_pass':len(rs)==6 and all(r['strict_eligible'] for r in rs),
            'sensitivity_full_panel_pass':len(rs)==6 and all(r['sensitivity_eligible'] for r in rs),
            'positive_strict_pass':len(positive)==2 and all(r['strict_eligible'] for r in positive),
            'positive_sensitivity_pass':len(positive)==2 and all(r['sensitivity_eligible'] for r in positive)}
    summaries=[]
    for c in m['cases']:
        if c['kind']!='novel_design':continue
        rs=[r for r in rows if r['id']==c['id']];control=controls[c['condition']]
        strict=len(rs)==2 and all(r['strict_eligible'] for r in rs)
        direction=len(rs)==2 and bool(c['histidine_positions']) and all(r['direction']['direction_pass'] for r in rs)
        tag_clear=len(rs)==2 and all(r['raw_tag_clear'] and r.get('prepared_tag_clear',False) for r in rs)
        summaries.append({'id':c['id'],'base_id':c['base_id'],'condition':c['condition'],'strict_eligible_both':strict,
            'sensitivity_eligible_both':len(rs)==2 and all(r['sensitivity_eligible'] for r in rs),
            'legacy_direction_both':direction,'tag_clear_both':tag_clear,
            'legacy_fully_qualified':bool(control['raw_pass'] and control['strict_full_panel_pass'] and strict and direction and tag_clear),
            'models':rs})
    paired=[]
    for base in [RESERVE,'apical_short_008_A41H']:
        for seed in m['refolding']['seeds']:
            rs={r['condition']:r for r in rows if r['base_id']==base and r['seed']==seed}
            p={'base_id':base,'seed':seed,'full_structure_and_controls_eligible':False}
            if len(rs)==2 and all(r['direction']['status']=='diagnostic' for r in rs.values()):
                bare,tagged=rs['bare'],rs['tagged']
                p.update(tagged_minus_bare_shift=tagged['direction']['shifts'][0]-bare['direction']['shifts'][0],
                    tagged_over_bare_His_only_ratio=tagged['direction']['his_only_Kd6_over_Kd7_4']/bare['direction']['his_only_Kd6_over_Kd7_4'],
                    full_structure_and_controls_eligible=all(r['strict_eligible'] and controls[k]['raw_pass'] and controls[k]['strict_full_panel_pass'] for k,r in rs.items()))
            paired.append(p)
    save(out/'direction.json',{'raw_results_sha256':sha(rawpath),'preparation_results_sha256':sha(preppath),
        'controls':controls,'rows':rows,'summaries':summaries,'paired_tag_effects':paired,'automatic_promotion':False,
        'interpretation':'Continuous held-out sensitivity, alongside unchanged legacy decisions. Diagnostic numbers do not rescue excluded models. Full suffix chemistry; only designed His41 contributes to the His-only proxy. No measured affinity or pH release.'})
    print(json.dumps({'controls':controls,'summaries':[{k:v for k,v in r.items() if k!='models'} for r in summaries]},indent=2))


if __name__=='__main__':main()
