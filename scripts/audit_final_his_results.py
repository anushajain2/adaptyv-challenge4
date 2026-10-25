"""Verify final-panel cohort, exclusions, pKa evidence and archived confidence."""
import json
import numpy as np
from novelty_search import ROOT,sha,save,now
from direction_check import histidines
from analyze_receptor_preparation_v2 import direction_metrics
from summarize_direction import terminal_distance

def main():
    out=ROOT/'reports/final_his_panel';raw=json.loads((out/'results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text());direction=json.loads((ROOT/'reports/final_his_direction/results.json').read_text());bundle=json.loads((out/'review_bundle/manifest.json').read_text())
    m=json.loads((ROOT/'data/final_his_panel/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    assert len(raw['rows'])==160 and len(m['cases'])==16
    expected={(c['id'],s,i) for c in m['cases'] for s in [2503,2609] for i in range(5)}
    assert {(r['id'],r['seed'],r['sample']) for r in raw['rows']}==expected
    selected={(r['id'],r['seed'],r['sample']) for r in raw['rows'] if r['confidence_selected']}
    assert {(p['id'],p['seed'],p['sample']) for p in prep['rows']}==selected
    assert {(p['id'],p['seed'],p['sample']) for p in bundle['models']}==selected
    assert prep['source_results_sha256']==sha(out/'results.json')
    for p in bundle['models']:
        for name,digest in p['files'].items():assert sha(out/'review_bundle'/f"{p['id']}_seed{p['seed']}"/name)==digest
    assert len(bundle['all_sample_confidence'])==32
    for p in bundle['all_sample_confidence']:
        path=ROOT/p['path'];assert sha(path)==p['sha256']
        with np.load(path,allow_pickle=False) as data:
            assert data['pae'].shape[0]==5 and data['pae'].shape[1]==data['pae'].shape[2]
            assert np.isfinite(data['pae']).all()
    count=0
    for v in direction['variants']:
        assert v['id'] in cases and len(v['models'])<=2
        for r in v['models']:
            d=r['direction']
            if d['status']!='ok':assert not d['direction_pass'];continue
            assert v['raw_pass_both_seeds'] and r['prepared_eligible']
            path=ROOT/r['source'];assert sha(path)==r['source_sha256']
            ps=[int(x[1:-1]) for x in v['mutations']]
            f=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',ps)
            b=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',ps)
            check=direction_metrics(f,b,terminal_distance(path,ps))
            for key in check:assert d[key]==check[key],(v['id'],key)
            count+=1
        assert v['fully_qualified_direction_pass_both_seeds']==bool(direction['prepared_controls_pass'] and v['matched_parent_prepared_pass_both_seeds'] and v['direction_pass_both_seeds'])
    assert len(direction['tag_selected_ids'])<=2
    assert all(v['fully_qualified_direction_pass_both_seeds'] for v in direction['variants'] if v['id'] in direction['tag_selected_ids'])
    paths=[p for folder in ['final_his_panel','final_his_direction'] for p in (ROOT/'reports'/folder).rglob('*') if p.is_file() and p.name!='round_evidence_audit.json']
    save(out/'round_evidence_audit.json',{'status':'PASS','checked_utc':now(),'predictions':160,'selected_models':len(selected),'full_case_confidence_arrays':32,'direction_models_recomputed':count,'no_unselected_pose_rescue':True,'source_and_bundle_hashes_match':True,'tag_selected_ids':direction['tag_selected_ids'],'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}})
    print('Final evidence verified:',len(selected),'selected models;',count,'interpretable pKa models.')
if __name__=='__main__':main()
