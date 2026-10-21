import copy
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_controls import ROOT,read_structure
from receptor_preparation_v2 import displacement,internal_geometry
from analyze_receptor_preparation_v2 import direction_metrics
from prepare_assay_tag_review import SUFFIX


def test_backbone_movement_guard_uses_original_frame_and_both_chains():
    source=read_structure(ROOT/'data/expansion_his_confirmation/apical_short_008_A41H.pdb')
    assert displacement(source,source)['pass']
    for chain in ('A','B'):
        moved=copy.deepcopy(source)
        for atom in moved[0][chain].get_atoms():atom.coord+=np.array([1.1,0,0])
        assert not displacement(source,moved)['pass']


def test_internal_clash_screen_excludes_covalent_pairs_and_detects_overlap(tmp_path):
    def atom(n,name,x,y):
        return f'ATOM  {n:5d} {name:^4} ALA A   1    {x:8.3f}{y:8.3f}{0.:8.3f}  1.00 20.00          {name[0]:>2}\n'
    p=tmp_path/'tiny.pdb'
    p.write_text(atom(1,'N',0,0)+atom(2,'CA',1.45,0)+atom(3,'C',2.3,1.2)+atom(4,'O',3.5,1.2)+'END\n')
    assert internal_geometry(p)['pass']
    p.write_text(atom(1,'N',0,0)+atom(2,'CA',1.45,0)+atom(3,'C',2.3,1.2)+atom(4,'O',.5,0)+'END\n')
    assert not internal_geometry(p)['pass']


def test_original_direction_heuristics_not_relaxed():
    assert direction_metrics([6.4],[5.3],20)['direction_pass']
    result=direction_metrics([7.2],[6.1],20)
    assert not result['direction_pass'] and not result['checks']['any_free_pka_6_to_7']
    assert not direction_metrics([6.4],[5.3],9)['direction_pass']


def test_documented_suffix_does_not_alter_binder_core():
    path=ROOT/'data/assay_tag_review/manifest.json'
    data=json.loads(path.read_text())
    assert data['suffix']==SUFFIX and len(data['cases'])==14
    for case in data['cases']:
        assert case['tagged_sequence']==case['bare_sequence']+SUFFIX
        assert case['tagged_sequence'][:case['bare_length']]==case['bare_sequence']
