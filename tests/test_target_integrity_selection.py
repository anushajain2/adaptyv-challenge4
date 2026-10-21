import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from target_integrity_selection import select_sample, target_ca_rmsd


def row(i, rmsd, confidence, binder_pass=False):
    return {'sample': i, 'target_CA_rmsd_to_input_A': rmsd,
            'metrics': {'iptm': confidence, 'ptm': confidence},
            'screen': {'pass': binder_pass}}


def test_deformed_high_confidence_sample_is_excluded():
    assert select_sample([row(0, .4, .63), row(1, 5.3, .65)]) == 0


def test_no_geometry_cherry_picking():
    assert select_sample([row(0, .4, .3, True), row(1, .5, .6, False)]) == 1


def test_no_eligible_sample_has_no_fallback():
    assert select_sample([row(0, 2.01, .7), row(1, math.nan, .9)]) is None


def test_boundary_and_tie_are_deterministic():
    assert select_sample([row(4, 2., .6), row(2, 1., .6)]) == 2


def test_negative_case_uses_same_selection():
    assert select_sample([row(0, .4, .12), row(1, .5, .15)]) == 1


def test_target_fit_is_rigid_transform_invariant_and_ignores_binder():
    import copy
    import numpy as np
    from prepare_controls import ROOT, read_structure
    ref = read_structure(ROOT / 'data/controls_matched/2ds25_5_6wrw.pdb')
    pred = copy.deepcopy(ref)
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    for atom in pred.get_atoms():
        atom.coord = atom.coord @ rotation + [12., 27., -4.]
    for atom in pred[0]['B'].get_atoms():
        atom.coord += 500
    assert target_ca_rmsd(pred, ref) < 1e-5
