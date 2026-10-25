import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT,sha
from prepare_controls import read_structure,sequence
from prepare_reserve_tag import RESERVE,ABLATION
from launch_reserve_tag import BOUND
from analyze_tag_comparison import core_arrays,cross_pae


def test_reserve_panel_identity_and_full_construct_masks():
    src=ROOT/'data/reserve_tag';m=json.loads((src/'manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    assert len(cases)==14 and m['refolding']['seeds']==[3509,3607] and BOUND<2
    assert m['core_confidence']==json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())['core_confidence']
    for c in cases.values():
        assert c['binder_sequence']==c['core_sequence']+(m['suffix'] if c['condition']=='tagged' else '')
        assert sequence(read_structure(src/f"{c['id']}.pdb")[0]['B'])==c['binder_sequence']
        with np.load(src/f"{c['id']}.npz") as data:
            mask=data['design_mask'];assert len(mask)==194+len(c['binder_sequence']) and not mask[:194].any() and mask[194:].all()
        assert sha(ROOT/c['reference'])==c['reference_sha256']
        assert c['histidine_positions']==([41] if c['base_id'] in [RESERVE,'apical_short_008_A41H'] else [])
    for name,digest in m['input_sha256'].items():assert sha(src/name)==digest


def test_ablation_changes_only_his41_and_tags_are_identical():
    m=json.loads((ROOT/'data/reserve_tag/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    for condition in ['bare','tagged']:
        a=cases[RESERVE+'_'+condition]['binder_sequence'];b=cases[ABLATION+'_'+condition]['binder_sequence']
        assert [(i+1,x,y) for i,(x,y) in enumerate(zip(a,b)) if x!=y]==[(41,'H','A')]
    assert 'histidine_positions' in cases[ABLATION+'_tagged'] and not cases[ABLATION+'_tagged']['histidine_positions']


def test_core_confidence_does_not_include_suffix_arrays():
    n=194+58;pae=np.ones((n+64,n+64))*30;pae[:n,:n]=5;plddt=np.ones(n+64)
    cp,_=core_arrays(pae,plddt,58)
    assert cross_pae(cp,58)==5 and cross_pae(pae,122)>5
