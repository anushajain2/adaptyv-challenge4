import copy
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_controls import contact_pairs, pose_metrics
from audit_control_samples import restore_predicted_atoms, sample_structures
from prepare_controls import make_complex, read_structure, save_input_cif, sequence
from score_matched_controls import interface_pae, run_ipsae


def reference():
    return read_structure(ROOT / "data/controls/3ds18_crystal.pdb")


def test_interface_pose_does_not_hide_displacement_by_realignment():
    ref = reference()
    moved = copy.deepcopy(ref)
    for atom in moved[0]["B"].get_atoms():
        atom.coord = atom.coord + np.array([12., 0., 0.])
    scores = pose_metrics(moved, ref)
    assert scores["binder_CA_fold_rmsd_A"] < 1e-4
    assert scores["binder_interface_CA_rmsd_after_target_fit_A"] == pytest.approx(12, abs=1e-4)
    assert scores["binder_CA_pose_rmsd_after_target_fit_A"] == pytest.approx(12, abs=1e-4)


def test_all_sample_reader_checks_coordinate_correspondence(tmp_path):
    s = reference()
    xyz = np.array([a.coord for a in s.get_atoms()])
    (tmp_path / "refold_cif").mkdir()
    (tmp_path / "fold_out_npz").mkdir()
    save_input_cif(s, tmp_path/"refold_cif/test.cif", sequence(s[0]["A"]), sequence(s[0]["B"]))
    coords = np.stack([xyz, xyz+10])
    fields = {key: np.array([1., 0.]) for key in
              ("iptm", "ptm", "design_to_target_iptm", "interaction_pae",
               "min_design_to_target_pae", "design_ptm")}
    path = tmp_path/"fold_out_npz/test.npz"
    np.savez(path, coords=coords, atom_resolved_mask=np.ones((1, len(xyz)), bool), **fields)
    samples = list(sample_structures(tmp_path, "test"))
    assert len(samples) == 2 and samples[0][3] and not samples[1][3]
    assert np.allclose(next(samples[1][1].get_atoms()).coord, xyz[0]+10)
    coords[0] = coords[0][::-1]
    np.savez(path, coords=coords, atom_resolved_mask=np.ones((1, len(xyz)), bool), **fields)
    with pytest.raises(AssertionError, match="atom order mismatch"):
        list(sample_structures(tmp_path, "test"))


def test_contact_pae_respects_both_directions_and_empty_interface():
    s = reference()
    n = len(s[0]["A"])
    pae = np.full((n+len(s[0]["B"]),)*2, 30.)
    for a, b in contact_pairs(s):
        pae[a-1, n+b-1], pae[n+b-1, a-1] = 2., 6.
    scores = interface_pae(s, pae)
    assert scores["contact_pae_both_directions_A"] == 4
    assert scores["contact_pae_A_to_B_A"] == 2
    for atom in s[0]["B"].get_atoms():
        atom.coord += np.array([1000., 0., 0.])
    assert interface_pae(s, pae)["contact_pae_both_directions_A"] is None


def test_official_ipsae_distinguishes_synthetic_uncertain_and_confident_pairs(tmp_path):
    s = reference()
    n = len(s[0]["A"])
    size = n+len(s[0]["B"])
    pae = np.full((size, size), 30.)
    high = run_ipsae(s, pae, np.full(size, .9), tmp_path/"uncertain")
    assert high["ipsae_min"] == 0
    assert next(read_structure(tmp_path/"uncertain/model.pdb").get_atoms()).bfactor == 90
    for a, b in contact_pairs(s):
        pae[a-1, n+b-1] = pae[n+b-1, a-1] = 1.
    low = run_ipsae(s, pae, np.full(size, .9), tmp_path/"confident")
    assert 0 < low["ipsae_min"] <= low["ipsae_max"] <= 1


def test_matched_panel_has_identical_target_coordinates_within_each_receptor():
    manifest = json.loads((ROOT/"data/controls_matched/manifest.json").read_text())
    for receptor in ("6wrv", "6wrw"):
        cases = [c for c in manifest["cases"] if c["receptor_id"] == receptor]
        assert {c["binder_id"] for c in cases} == {"3ds18", "2ds25_5", "scrambled", "caplacizumab"}
        targets = [np.array([a.coord for a in read_structure(ROOT/c["pdb"])[0]["A"].get_atoms()])
                   for c in cases]
        assert all(np.array_equal(targets[0], xyz) for xyz in targets)


def test_predicted_sidechain_restoration_is_checked_not_imputed():
    ref = reference()
    s = make_complex([ref[0]["B"][1]], [ref[0]["B"][1]])
    order = json.loads((ROOT/"configs/boltzgen_atom_order.json").read_text())["ref_atoms"]
    residues = [s[0]["A"][1], s[0]["B"][1]]
    coords = np.array([r[name].coord for r in residues for name in order[r.resname]])
    mapping = np.zeros((len(coords), 2), bool)
    n = len(order[residues[0].resname])
    mapping[:n, 0] = True
    mapping[n:, 1] = True
    expected = residues[1]["CB"].coord.copy()
    residues[1].detach_child("CB")
    assert restore_predicted_atoms(s, coords, mapping) == 1
    assert np.array_equal(residues[1]["CB"].coord, expected)
    with pytest.raises(AssertionError, match="canonical atom order mismatch"):
        restore_predicted_atoms(s, coords[::-1], mapping)
