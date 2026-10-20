import copy
import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from prepare_controls import read_structure, sequence
from launch_generation_comparison import gross_cost
from analyze_generation_comparison import screen_sample, design_pose_metrics


def test_arms_differ_only_in_sheet_constraint_and_have_no_binder_template():
    folder = ROOT/"data/generation_comparison"
    a = yaml.safe_load((folder/"unconstrained.yaml").read_text())
    b = yaml.safe_load((folder/"sheet.yaml").read_text())
    assert b["entities"][1]["protein"].pop("secondary_structure") == {"sheet":"33..37"}
    assert a == b
    s = read_structure(folder/"target_6wrw.cif")
    assert [c.id for c in s[0]] == ["A"]
    ref = read_structure(ROOT/"data/controls_matched/2ds25_5_6wrw.pdb")
    assert sequence(s[0]["A"]) == sequence(ref[0]["A"])
    assert np.allclose([a.coord for a in s[0]["A"].get_atoms()],
                       [a.coord for a in ref[0]["A"].get_atoms()])
    entries = a["entities"][0]["file"]["binding_types"]
    assert not any("binding" in e["chain"] and "not_binding" in e["chain"] for e in entries)
    merged = {k:v for e in entries for k,v in e["chain"].items()}
    assert not (set(merged["binding"].split(",")) & set(merged["not_binding"].split(",")))


def test_no_file_spec_silently_overrides_hotspots_with_exclusions():
    for path in (ROOT/"configs").glob("*.yaml"):
        spec = yaml.safe_load(path.read_text())
        for entity in spec["entities"]:
            for entry in entity.get("file",{}).get("binding_types",[]):
                assert not {"binding","not_binding"} <= set(entry["chain"]), path.name


def test_budget_counts_gross_not_credit_reduced_bill():
    assert gross_cost({"metered_cost":"3","billed_cost":"0",
                       "metered_cost_breakdown":{"gpu":"4","cpu":"1"}}) == 5
    assert gross_cost({"metered_cost":"6","metered_cost_breakdown":{"gpu":"4"}}) == 6


def test_structural_screen_rejects_wrong_pose_even_with_good_fold_and_confidence():
    plan = json.loads((ROOT/"data/generation_comparison/manifest.json").read_text())["screen"]
    row = {"metrics":{"design_to_target_iptm":.8},
           "pose":{"binder_CA_fold_rmsd_A":.5,"binder_interface_CA_rmsd_after_target_fit_A":.5,
                   "native_contact_recall":.9},
           "geometry":{"target_fit_rmsd_A":.5,"modeled_target_clash_pairs_lt2A":0,
                       "omitted_same_subunit_clash_pairs_lt2A":0,"partner_subunit_clash_pairs_lt2A":0,
                       "observed_glycans_clash_pairs_lt2A":0},
           "contacts":{"hotspots_contacted_crop":[19,21],"fraction_pairs_in_reference_region":.9}}
    assert screen_sample(row,plan)["pass"]
    moved = copy.deepcopy(row)
    moved["pose"]["binder_interface_CA_rmsd_after_target_fit_A"] = 12
    assert "interface_pose" in screen_sample(moved,plan)["failed"]
    invalid_context = copy.deepcopy(row)
    invalid_context["geometry"]["target_fit_rmsd_A"] = 3
    assert not screen_sample(invalid_context,plan)["pass"]
    clash = copy.deepcopy(row)
    clash["geometry"]["partner_subunit_clash_pairs_lt2A"] = 1
    assert not screen_sample(clash,plan)["pass"]


def test_inverse_fold_placeholder_sidechains_cannot_define_design_contacts():
    ref = read_structure(ROOT/"data/controls_matched/2ds25_5_6wrw.pdb")
    placeholder = copy.deepcopy(ref)
    for r in placeholder[0]["B"]:
        for atom in r:
            if atom.name not in {"N","CA","C","O"}: atom.coord = np.zeros(3)
    metrics = design_pose_metrics(ref,placeholder)
    assert metrics["native_contact_recall"] == 1
    assert metrics["binder_interface_CA_rmsd_after_target_fit_A"] < 1e-4
