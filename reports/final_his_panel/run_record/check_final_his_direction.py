"""Apply unchanged raw, preparation, control and pKa gates to final His panel."""
import csv,json
from novelty_search import ROOT,sha,save
from prepare_controls import read_structure
from analyze_controls import pose_metrics
from analyze_generation_comparison import design_pose_metrics,contacts_summary,screen_sample
from analyze_receptor_preparation_v2 import direction_metrics
from direction_check import histidines
from summarize_direction import terminal_distance

def main():
    root=ROOT/'reports/final_his_panel';raw=json.loads((root/'results.json').read_text());prep=json.loads((root/'preparation/results.json').read_text())
    manifest=json.loads((ROOT/'data/final_his_panel/manifest.json').read_text());cases={c['id']:c for c in manifest['cases']};arms={a['id']:a for a in manifest['arms']}
    assert prep['source_results_sha256']==sha(root/'results.json')
    selected={(r['id'],r['seed'],r['sample']):r for r in raw['rows'] if r['confidence_selected']}
    assert set(selected)=={(p['id'],p['seed'],p['sample']) for p in prep['rows']} and len(selected)==len(prep['rows'])
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    summaries={c['id']:c for c in raw['candidate_summaries']};rows=[]
    for p in sorted(prep['rows'],key=lambda p:(p['id'],p['seed'])):
        old=selected[(p['id'],p['seed'],p['sample'])];case=cases[p['id']]
        r={'id':p['id'],'seed':p['seed'],'sample':p['sample'],'kind':case['kind'],'prepared_eligible':False,'direction':{'status':'excluded','direction_pass':False}}
        try:
            assert sha(ROOT/p['source'])==p['source_sha256']
            assert p['status']=='ok','Preparation execution failed'
            path=ROOT/p['dest']/'prepared.pdb';assert sha(path)==p['prepared_sha256']
            s=read_structure(path);ref=read_structure(ROOT/'data/final_his_panel'/f"{p['id']}.pdb")
            if case['kind']=='novel_design':
                check=screen_sample({'metrics':old['metrics'],'pose':design_pose_metrics(s,ref),'geometry':p['prepared_geometry'],'contacts':contacts_summary(s,arms[case['arm']],conserved)},manifest['screen'])
            elif case['kind']=='positive':
                q=pose_metrics(s,ref);check={'pass':old['metrics']['design_to_target_iptm']>=.5 and q['binder_CA_fold_rmsd_A']<=2 and q['binder_interface_CA_rmsd_after_target_fit_A'] is not None and q['binder_interface_CA_rmsd_after_target_fit_A']<=2 and q['native_contact_recall']>=.5,'pose':q}
            else:check={'pass':old['control_pass'],'note':'Preserve raw presumed-negative confidence check; preparation does not recompute confidence.'}
            r.update(prepared_screen=check,displacement=p['displacement'],internal_geometry=p['prepared_internal_geometry'],propka=p['propka'],source=str(path.relative_to(ROOT)),source_sha256=sha(path))
            assert old.get('screen',{}).get('pass',old.get('control_pass')),'Raw selected screen failed'
            assert p['displacement']['pass'] and p['prepared_internal_geometry']['pass'],'Prepared movement/internal geometry failed'
            assert check['pass'],'Prepared structural screen failed'
            assert all(v['status']=='ok' for v in p['propka'].values()),'PROPKA incomplete/protonation group failure'
            r['prepared_eligible']=True
            if case.get('mutations'):
                assert raw['controls_pass'] and summaries[case['id']]['all_selected_seeds_full_screen'],'Both raw poses and controls required'
                assert summaries[case['parent_id']]['all_selected_seeds_full_screen'],'Matched parent raw check failed'
                positions=[int(x[1:-1]) for x in case['mutations']]
                free=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',positions)
                bound=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',positions)
                r['direction']={'status':'ok',**direction_metrics(free,bound,terminal_distance(path,positions))}
        except Exception as e:r['exclusion']=f'{type(e).__name__}: {e}'
        rows.append(r)
    controls=[r for r in rows if r['kind']!='novel_design'];cp=len(controls)==6 and all(r['prepared_eligible'] for r in controls)
    variants=[]
    for c in manifest['cases']:
        if not c.get('mutations'):continue
        rs=[r for r in rows if r['id']==c['id']];parents=[r for r in rows if r['id']==c['parent_id']]
        parentpass=len(parents)==2 and all(r['prepared_eligible'] for r in parents)
        dp=len(rs)==2 and all(r['direction']['direction_pass'] for r in rs)
        variants.append({'id':c['id'],'sequence':c['binder_sequence'],'parent_id':c['parent_id'],'mutations':c['mutations'],'raw_pass_both_seeds':summaries[c['id']]['all_selected_seeds_full_screen'],'matched_parent_prepared_pass_both_seeds':parentpass,'direction_pass_both_seeds':dp,'fully_qualified_direction_pass_both_seeds':cp and parentpass and dp,'models':rs})
    eligible=[v for v in variants if v['fully_qualified_direction_pass_both_seeds']]
    eligible.sort(key=lambda v:(-min(r['direction']['his_only_Kd6_over_Kd7_4'] for r in v['models']),-min(r['metrics']['design_to_target_iptm'] for r in summaries[v['id']]['selected_models']),v['id']))
    result={'source_raw_sha256':sha(root/'results.json'),'source_preparation_sha256':sha(root/'preparation/results.json'),'raw_controls_pass':raw['controls_pass'],'prepared_controls_pass':cp,'rows':rows,'variants':variants,'tag_selected_ids':[v['id'] for v in eligible[:2]],'limitations':['Fixed bound conformation used for free/bound comparison','His-only independent-site ratio is not measured affinity','Cropped target and no unbound ensemble','No experimental switch, mouse or tag claim']}
    save(ROOT/'reports/final_his_direction/results.json',result)
    print(json.dumps({'prepared_controls_pass':cp,'eligible_prepared':sum(r['prepared_eligible'] for r in rows),'raw_variant_passes':sum(v['raw_pass_both_seeds'] for v in variants),'direction_passes':sum(v['direction_pass_both_seeds'] for v in variants),'tag_selected_ids':result['tag_selected_ids']},indent=2))
if __name__=='__main__':main()
