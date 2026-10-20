"""Compare blind control refolding with experimental geometry and negative controls."""
import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import Superimposer
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

from prepare_controls import ROOT, geometry, read_structure, residues, sequence


def contact_pairs(s, cutoff=5.):
    target = [(r.id[1], a.coord) for r in s[0]["A"] for a in r if a.element != "H"]
    xyz = np.array([x[1] for x in target])
    tree = cKDTree(xyz)
    pairs = set()
    for r in s[0]["B"]:
        for neighbors in tree.query_ball_point(
                np.array([a.coord for a in r if a.element != "H"]), cutoff):
            pairs.update((target[i][0], r.id[1]) for i in neighbors)
    return pairs


def pose_metrics(predicted, reference):
    """Use only observed reference residues; preserve target-relative binder placement."""
    common = [r.id[1] for r in reference[0]["A"]
              if "CA" in r and (" ", r.id[1], " ") in predicted[0]["A"]]
    sup = Superimposer()
    sup.set_atoms([reference[0]["A"][i]["CA"] for i in common],
                  [predicted[0]["A"][i]["CA"] for i in common])
    target_rms = float(sup.rms)
    bind = [r.id[1] for r in reference[0]["B"]
            if "CA" in r and (" ", r.id[1], " ") in predicted[0]["B"]]
    xyz = np.array([predicted[0]["B"][i]["CA"].coord for i in bind])
    ref = np.array([reference[0]["B"][i]["CA"].coord for i in bind])
    rot, trans = sup.rotran
    ligand_rms = float(np.sqrt(((xyz @ rot+trans-ref)**2).sum(axis=1).mean()))
    native_contacts = contact_pairs(reference)
    interface_bind = sorted({b for _, b in native_contacts} & set(bind))
    predicted_interface = np.array([predicted[0]["B"][i]["CA"].coord for i in interface_bind])
    reference_interface = np.array([reference[0]["B"][i]["CA"].coord for i in interface_bind])
    interface_rms = float(np.sqrt(((predicted_interface @ rot+trans-reference_interface)**2)
                                 .sum(axis=1).mean())) if interface_bind else None
    sup.set_atoms([reference[0]["B"][i]["CA"] for i in bind],
                  [predicted[0]["B"][i]["CA"] for i in bind])
    recovered = contact_pairs(predicted) & native_contacts
    return {
        "target_CA_rmsd_to_reference_A": target_rms,
        "binder_CA_pose_rmsd_after_target_fit_A": ligand_rms,
        "binder_interface_CA_rmsd_after_target_fit_A": interface_rms,
        "binder_CA_fold_rmsd_A": float(sup.rms),
        "native_contacts": len(native_contacts), "recovered_contacts": len(recovered),
        "native_contact_recall": len(recovered)/len(native_contacts) if native_contacts else None,
    }


def analyze(run):
    manifest = json.loads((run / "inputs/manifest.json").read_text())
    cases = {c["id"]: c for c in manifest["cases"]}
    full = read_structure(ROOT / "data/raw/11QC.cif")
    rows = []
    for seed in manifest["refolding"]["seeds"]:
        dest = run / f"seed-{seed}"
        summary = json.loads((dest / "confidence_summary.json").read_text())
        preflight = json.loads((dest / "preflight.json").read_text())
        assert {p["id"] for p in preflight} == set(cases)
        for metrics in summary:
            case = cases[metrics["id"]]
            predicted = read_structure(dest / "refold_cif" / f"{case['id']}.cif")
            reference_input = read_structure(ROOT / case["pdb"])
            observed = {r.id[1] for r in predicted[0]["A"]}
            expected_observed = {r.id[1] for r in reference_input[0]["A"]}
            assert observed == expected_observed, case["id"]
            assert all(seq1(r.resname) == case["target_sequence"][r.id[1]-1]
                       for r in predicted[0]["A"]), case["id"]
            assert sequence(predicted[0]["B"]) == case["binder_sequence"], case["id"]
            selected = metrics["selected_sample"]
            iptm = np.asarray(metrics["design_to_target_iptm"]).reshape(-1)
            row = {"id": case["id"], "kind": case["kind"], "seed": seed,
                   "selected_sample": selected,
                   "selected_interface_iptm": float(iptm[selected]),
                   "all_interface_iptm": iptm.tolist(),
                   "unresolved_target_positions_omitted_from_export": sorted(
                       set(range(1, len(case["target_sequence"])+1))-observed),
                   "geometry": geometry(predicted,
                       {int(k): v for k, v in case["target_native_map"].items()
                        if int(k) in observed}, full)}
            row["full_context_fit_valid"] = row["geometry"]["target_fit_rmsd_A"] <= 2
            if case["kind"].startswith("positive"):
                # Current-target cases are compared to the experimental binding mode,
                # with corresponding native receptor numbering, not a forced pose.
                ref_id = {"3ds18_current": "3ds18_crystal"}.get(case["id"], case["id"])
                reference = read_structure(ROOT / cases[ref_id]["pdb"])
                row["pose"] = pose_metrics(predicted, reference)
                row["interface_geometry_recovered"] = (
                    row["pose"]["binder_interface_CA_rmsd_after_target_fit_A"] <= 2
                    and row["pose"]["native_contact_recall"] >= .5)
                row["diagnostic_recovery"] = (
                    row["selected_interface_iptm"] >= .6
                    and row["pose"]["binder_CA_pose_rmsd_after_target_fit_A"] <= 5
                    and row["pose"]["native_contact_recall"] >= .3)
            rows.append(row)
    assert len(rows) == len(cases)*len(manifest["refolding"]["seeds"])
    initial = [r for r in rows if r["id"] == "3ds18_current"]
    negatives = [r for r in rows if r["kind"] == "negative_hypothesis"]
    max_neg = max(r["selected_interface_iptm"] for r in negatives)
    min_pos = min(r["selected_interface_iptm"] for r in initial)
    calibrated = all(r["diagnostic_recovery"] for r in initial) and min_pos-max_neg >= .2
    return {
        "run": run.name, "rows": rows,
        "decision": {
            "3ds18_current_recovers_both_seeds": all(r["diagnostic_recovery"] for r in initial),
            "minimum_positive_minus_maximum_negative_iptm": min_pos-max_neg,
            "small_followup_supported": calibrated,
            "interface_geometry_flag_scope": "Descriptive interface-only recovery (CA RMSD <=2 A and contact recall >=0.5); does not override the prospective generation gate.",
            "gate_scope": "Prospective diagnostic heuristic, not validated proof of binding or universal pipeline accuracy.",
        },
        "limitations": [
            "These published structures may occur in predictor training data.",
            "Negative controls are hypothesized nonbinders, not experimentally verified nonbinders.",
            "Five diffusion samples per seed; confidence-selected structure evaluated for pose.",
            f"Target-only templates; MSA condition: {manifest['refolding']['msa']}.",
            "ATV uses observed Fc 237–444, with a cropped circularly permuted receptor; glycosylation omitted in prediction.",
            "atv_current pose reference is the transferred experimental pose; it is not an independent crystal structure.",
        ],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--out", type=Path, default=ROOT / "reports/controls")
    args = p.parse_args()
    result = analyze(args.run)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(result, indent=2))
    for row in result["rows"]:
        print(row["id"], row["seed"], round(row["selected_interface_iptm"], 3),
              row.get("diagnostic_recovery"), row.get("pose", {}))
    print(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    main()
