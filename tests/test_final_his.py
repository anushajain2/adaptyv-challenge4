import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT,sha
from launch_final_his import PHASES,resource_bound
from prepare_controls import read_structure,sequence

def test_final_panel_preserves_anchors_and_uses_two_fresh_seeds():
    p=ROOT/'data/final_his_panel';m=json.loads((p/'manifest.json').read_text());a=json.loads((p/'contact_audit.json').read_text())
    assert m['refolding']['seeds']==[2503,2609] and m['refolding']['diffusion_samples_per_seed']==5
    assert m['screen']==json.loads((ROOT/'data/diversity_confirmation/manifest.json').read_text())['screen']
    parent=next(c for c in m['cases'] if c['id']=='diversity_long_037')['binder_sequence']
    variants=[c for c in m['cases'] if c.get('mutations')]
    assert len(variants)==12 and len(m['cases'])==16
    allowed=set(a['chosen_sites'])
    for c in variants:
        mutated={i+1 for i,(x,y) in enumerate(zip(parent,c['binder_sequence'])) if x!=y}
        assert mutated<allowed and all(parent[i-1] not in 'GPCHDEKR' and 5<i<=90 for i in mutated)
        assert sequence(read_structure(p/(c['id']+'.pdb'))[0]['B'])==c['binder_sequence']
    for name,digest in m['input_sha256'].items():assert sha(p/name)==digest

def test_final_round_resource_caps_fit_authorized_three_dollars():
    assert sum(p[2] for p in PHASES.values())==3
    assert all(resource_bound(p)<PHASES[p][2] for p in PHASES)
    text=(ROOT/'scripts/modal_final_his.py').read_text()
    assert 'timeout=1650, startup_timeout=120, retries=0' in text
    assert 'max_containers=1' in text
