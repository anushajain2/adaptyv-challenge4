import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT,sha

def test_final_direction_requires_both_poses_parent_and_controls():
    d=json.loads((ROOT/'reports/final_his_direction/results.json').read_text())
    raw=json.loads((ROOT/'reports/final_his_panel/results.json').read_text())
    assert len(raw['rows'])==160 and len(d['variants'])==12
    summaries={c['id']:c for c in raw['candidate_summaries']}
    for v in d['variants']:
        assert v['raw_pass_both_seeds']==summaries[v['id']]['all_selected_seeds_full_screen']
        assert v['direction_pass_both_seeds']==(len(v['models'])==2 and all(r['direction']['direction_pass'] for r in v['models']))
        assert v['fully_qualified_direction_pass_both_seeds']==bool(d['prepared_controls_pass'] and v['matched_parent_prepared_pass_both_seeds'] and v['direction_pass_both_seeds'])
        if not v['raw_pass_both_seeds']:assert all(r['direction']['status']=='excluded' for r in v['models'])
    assert len(d['tag_selected_ids'])<=2
    assert set(d['tag_selected_ids'])<={v['id'] for v in d['variants'] if v['fully_qualified_direction_pass_both_seeds']}

def test_final_evidence_includes_failed_selected_models_and_confidence_arrays():
    out=ROOT/'reports/final_his_panel'
    audit=json.loads((out/'round_evidence_audit.json').read_text())
    assert audit['status']=='PASS' and audit['predictions']==160
    for path,digest in audit['files'].items():assert sha(ROOT/path)==digest
    bundle=json.loads((out/'review_bundle/manifest.json').read_text())
    assert len(bundle['all_sample_confidence'])==32
    raw=json.loads((out/'results.json').read_text())
    assert {(p['id'],p['seed'],p['sample']) for p in bundle['models']}=={(p['id'],p['seed'],p['sample']) for p in raw['rows'] if p['confidence_selected']}
