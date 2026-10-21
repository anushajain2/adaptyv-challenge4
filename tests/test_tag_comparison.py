import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_controls import ROOT,read_structure,sequence
from prepare_assay_tag_review import SUFFIX
from analyze_tag_comparison import core_arrays,cross_pae,select_core,core_structure,core_screen
from prepare_tag_models import displacement

def test_suffix_cannot_change_core_confidence_or_sample_selection():
    rng=np.random.default_rng(7);n=194+58;pae=rng.uniform(0,10,(n+64,n+64));plddt=np.ones(n+64)
    a,b=core_arrays(pae,plddt,58);expected=cross_pae(a,58)
    pae[n:,:]=99;pae[:,n:]=99
    assert cross_pae(core_arrays(pae,plddt,58)[0],58)==expected
    rows=[{'sample':i,'target_CA_rmsd_to_input_A':1,'core_mean_cross_PAE_A':x} for i,x in enumerate([14,12,12])]
    assert select_core(rows)==1
    rows[1]['target_CA_rmsd_to_input_A']=3;assert select_core(rows)==2

def test_tag_inputs_preserve_sequence_target_and_core_coordinates():
    m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())
    assert len(m['cases'])==14 and len(m['refolding']['seeds'])==2
    for c in m['cases']:
        s=read_structure(ROOT/'data/tag_comparison'/f"{c['id']}.pdb");ref=read_structure(ROOT/c['reference'])
        assert sequence(s[0]['B'])==c['core_sequence']+(SUFFIX if c['condition']=='tagged' else '')
        core=core_structure(s,c['core_length'])
        for chain in ('A','B'):
            for r in ref[0][chain]:
                for a in r:np.testing.assert_array_equal(core[0][chain][r.id][a.name].coord,a.coord)
        with np.load(ROOT/'data/tag_comparison'/f"{c['id']}.npz") as data:
            mask=data['design_mask'];assert len(mask)==194+len(c['binder_sequence'])
            assert not mask[:194].any() and mask[194:].all()

def test_flexible_tag_movement_does_not_hide_core_movement():
    import copy
    s=read_structure(ROOT/'data/tag_comparison/apical_short_008_A41H_tagged.pdb')
    moved=copy.deepcopy(s)
    for r in moved[0]['B']:
        if r.id[1]>58:
            for a in r:a.coord+=np.array([2.,0,0])
    assert displacement(s,moved,58)['pass']
    assert not displacement(s,moved)['pass']
    moved[0]['B'][41]['CA'].coord+=np.array([2.,0,0])
    assert not displacement(s,moved,58)['pass']

def test_whole_chain_iptm_is_replaced_without_dropping_geometry_checks():
    import copy
    m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())
    previous=json.loads((ROOT/'reports/expansion_his_confirmation/results.json').read_text())
    row=copy.deepcopy(next(r for r in previous['rows'] if r['id']=='apical_short_008_A41H' and r['seed']==1901 and r['confidence_selected']))
    row['metrics']['design_to_target_iptm']=.01;row['core_mean_cross_PAE_A']=14.
    assert core_screen(row,m)['pass']
    row['geometry']['modeled_target_clash_pairs_lt2A']=1
    assert not core_screen(row,m)['pass']
    row['geometry']['modeled_target_clash_pairs_lt2A']=0;row['core_mean_cross_PAE_A']=25.
    assert not core_screen(row,m)['pass']
