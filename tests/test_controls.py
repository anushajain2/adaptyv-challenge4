import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_controls import read_structure, sequence, geometry


def controls():
    return json.loads((ROOT / "data/controls/manifest.json").read_text())


def test_control_sequences_and_masks_cover_unresolved_residues():
    for case in controls()["cases"]:
        cif = MMCIF2Dict(str(ROOT / case["cif"]))
        sequences = cif["_entity_poly.pdbx_seq_one_letter_code_can"]
        assert sequences == [case["target_sequence"], case["binder_sequence"]]
        mask = np.load(ROOT / Path(case["cif"]).with_suffix(".npz"))["design_mask"]
        n = len(sequences[0])
        assert len(mask) == sum(map(len, sequences))
        assert not mask[:n].any() and mask[n:].all()
        for label, entity, residue, index in zip(
            cif["_atom_site.label_asym_id"], cif["_atom_site.label_entity_id"],
            cif["_atom_site.label_comp_id"], cif["_atom_site.label_seq_id"]
        ):
            from Bio.SeqUtils import seq1
            assert int(entity) == ("A", "B").index(label)+1
            assert seq1(residue) == sequences[int(entity)-1][int(index)-1]


def test_scrambled_control_preserves_composition():
    cases = {c["id"]: c for c in controls()["cases"]}
    a = cases["3ds18_current"]["binder_sequence"]
    b = cases["3ds18_scrambled"]["binder_sequence"]
    assert a != b and Counter(a) == Counter(b)


def test_published_geometry_and_atv_author_chain_pair():
    manifest = controls()
    assert manifest["atv_source_chains"] == {"target": "D", "binder": "A"}
    cases = {c["id"]: c for c in manifest["cases"]}
    for name in ["3ds18_crystal", "atv_permuted"]:
        g = cases[name]["geometry"]
        assert g["target_fit_rmsd_A"] < 1.1
        assert len(g["native_target_contacts_lt5A"]) >= 10
        assert all(g[k] == 0 for k in g if k.endswith("clash_pairs_lt2A"))
    atv = read_structure(ROOT / cases["atv_permuted"]["pdb"])
    assert len(atv[0]["A"]) == 142
    assert len(cases["atv_permuted"]["target_sequence"]) == 144
    mapping = cases["atv_permuted"]["target_native_map"]
    assert mapping["50"] == 379 and mapping["51"] == 194


def test_geometry_counter_detects_deliberate_overlap():
    case = next(c for c in controls()["cases"] if c["id"] == "3ds18_crystal")
    s = read_structure(ROOT / case["pdb"])
    binder_atom = next(s[0]["B"].get_atoms())
    binder_atom.coord = s[0]["A"][1]["CA"].coord.copy()
    result = geometry(s, {int(k): v for k, v in case["target_native_map"].items()},
                      read_structure(ROOT / "data/raw/11QC.cif"))
    assert result["modeled_target_clash_pairs_lt2A"] > 0


def test_pose_metrics_ignore_global_rigid_motion():
    import copy
    import pytest
    from analyze_controls import pose_metrics
    reference = read_structure(ROOT / "data/controls/3ds18_crystal.pdb")
    moved = copy.deepcopy(reference)
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    for atom in moved.get_atoms():
        atom.coord = atom.coord @ rotation + np.array([20., -5., 7.])
    result = pose_metrics(moved, reference)
    assert result["native_contact_recall"] == 1
    assert result["binder_interface_CA_rmsd_after_target_fit_A"] == pytest.approx(0, abs=1e-4)
    assert result["binder_CA_pose_rmsd_after_target_fit_A"] == pytest.approx(0, abs=1e-4)
