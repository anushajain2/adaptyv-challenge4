"""Score full PAE for every sample with the pinned official ipSAE implementation."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from analyze_controls import contact_pairs
from audit_control_samples import audit, sample_structures
from prepare_controls import ROOT, save_pdb

IPSAE = ROOT / "third_party/ipsae/ipsae.py"
IPSAE_COMMIT = "6174cf9e71cb1bd660cc805856a18c4871a6dec3"


def interface_pae(structure, pae):
    """Mean directional PAE over actual contacting residue pairs, in token order."""
    ordered = [(c.id, r.id[1]) for c in structure[0] for r in c if "CA" in r]
    assert pae.shape == (len(ordered), len(ordered))
    indices = {key: i for i, key in enumerate(ordered)}
    pairs = sorted(contact_pairs(structure))
    if not pairs:
        return {"contact_residue_pairs": 0, "contact_pae_A_to_B_A": None,
                "contact_pae_B_to_A_A": None, "contact_pae_both_directions_A": None}
    a = np.array([indices[("A", t)] for t, _ in pairs])
    b = np.array([indices[("B", u)] for _, u in pairs])
    ab, ba = float(pae[a, b].mean()), float(pae[b, a].mean())
    return {"contact_residue_pairs": len(pairs), "contact_pae_A_to_B_A": ab,
            "contact_pae_B_to_A_A": ba, "contact_pae_both_directions_A": (ab+ba)/2}


def run_ipsae(structure, pae, plddt, dest):
    dest.mkdir(parents=True, exist_ok=True)
    ordered_residues = [r for c in structure[0] for r in c if "CA" in r]
    assert len(ordered_residues) == len(plddt) == len(pae)
    scale = 100 if np.max(plddt) <= 1 else 1
    # The upstream selected CIF can carry another sample's B factors. Export
    # this sample's actual pLDDT for viewers as well as the ipSAE input array.
    for residue, confidence in zip(ordered_residues, plddt):
        for atom in residue:
            atom.bfactor = float(confidence * scale)
    save_pdb(structure, dest / "model.pdb")
    np.savez_compressed(dest / "pae_model.npz", pae=pae)
    np.savez_compressed(dest / "plddt_model.npz", plddt=plddt)
    # Do not invent pair-chain ipTM values: the upstream tool can run without
    # its optional confidence JSON. Use the original metrics for ipTM instead.
    result = subprocess.run(
        [sys.executable, str(IPSAE), str(dest/"pae_model.npz"),
         str(dest/"model.pdb"), "10", "10"],
        capture_output=True, text=True, timeout=60)
    (dest / "ipsae.log").write_text(result.stdout + result.stderr)
    result.check_returncode()
    lines = (dest / "model_10_10.txt").read_text().splitlines()
    header = next(line.split() for line in lines if line.startswith("Chn1"))
    parsed = [dict(zip(header, line.split())) for line in lines
              if line.startswith(("A ", "B "))]
    directions = {row["Chn1"]+"_"+row["Chn2"]: float(row["ipSAE"])
                  for row in parsed if row["Type"] == "asym"}
    assert set(directions) == {"A_B", "B_A"}
    return {"ipsae_directional": directions, "ipsae_min": min(directions.values()),
            "ipsae_max": max(directions.values())}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.run / "inputs/manifest.json").read_text())
    cases = {c["id"]: c for c in manifest["cases"]}
    result = audit(args.run, restore_sidechains=True)
    lookup = {(r["id"], r["seed"], r["sample"]): r for r in result["rows"]}
    for seed in manifest["refolding"]["seeds"]:
        folder = args.run / f"seed-{seed}"
        for case in cases.values():
            name = case["id"]
            with np.load(folder/"fold_out_npz"/f"{name}.npz", allow_pickle=False) as data:
                paes, plddts = data["pae"].copy(), data["plddt"].copy()
                asym = data["asym_id"][0]
                n = len(case["target_sequence"])
                assert np.all(asym[:n] == 0) and np.all(asym[n:] == 1), name
                assert data["token_pad_mask"].all(), name
            for i, structure, metrics, _ in sample_structures(folder, name, restore_sidechains=True):
                # Check scalar PAE against the exported matrix to catch sample,
                # chain, unit, and token-order mistakes in the export path.
                mean_interchain = float((paes[i, :n, n:].mean() + paes[i, n:, :n].mean())/2)
                assert abs(mean_interchain - metrics["interaction_pae"]) < .03, (name, i)
                assert len(structure[0]["A"]) == n
                assert len(structure[0]["B"]) == len(case["binder_sequence"])
                dest = ROOT / "outputs/matched_scoring" / args.run.name / f"seed-{seed}" / name / str(i)
                score = run_ipsae(structure, paes[i], plddts[i], dest)
                score.update(interface_pae(structure, paes[i]))
                score["binder_mean_plddt"] = float(plddts[i, n:].mean())
                score["target_mean_plddt"] = float(plddts[i, :n].mean())
                row = lookup[(name, seed, i)]
                row.update({k: case[k] for k in ("binder_id", "receptor_id", "role")})
                row["scores"] = score
                row["raw_scoring_output"] = str(dest.relative_to(ROOT))
            print(f"Scored seed {seed}: {name}", flush=True)
    result["analysis_plan"] = manifest["analysis_plan"]
    result["predicted_atom_export"] = (
        "Restored predicted sidechains omitted by input atom_resolved_mask, using saved "
        "coordinates and pinned canonical atom order. Every retained atom was independently "
        "checked against that order before restoration. No imputation or new prediction.")
    result["ipsae_provenance"] = {
        "repository": "https://github.com/DunbrackLab/IPSAE", "commit": IPSAE_COMMIT,
        "sha256": hashlib.sha256(IPSAE.read_bytes()).hexdigest(),
        "pae_cutoff_A": 10, "distance_cutoff_A": 10,
        "primary_score": "min(A_to_B, B_to_A)", "source_unmodified": True,
        "note": "The optional ipTM column in raw ipSAE output is unset; use original NPZ ipTM.",
    }
    # Train one provisional scalar threshold using only 3DS18 and matched
    # negatives. Evaluate the held-out positive separately; never retune.
    thresholds = {}
    for receptor in ("6wrv", "6wrw"):
        rows = [r for r in result["rows"] if r["receptor_id"] == receptor]
        negatives = [r for r in rows if r["kind"] == "negative_hypothesis"]
        threshold = max(r["scores"]["ipsae_min"] for r in negatives)
        checks = {}
        for binder in ("3ds18", "2ds25_5"):
            subset = [r for r in rows if r["binder_id"] == binder]
            checks[binder] = {
                "samples_above_negative_max": sum(r["scores"]["ipsae_min"] > threshold for r in subset),
                "interface_recovered_and_above": sum(
                    r["scores"]["ipsae_min"] > threshold and r["interface_recovered"] for r in subset),
                "whole_fold_and_interface_and_score": sum(
                    r["scores"]["ipsae_min"] > threshold and r["interface_recovered"]
                    and r["pose"]["binder_CA_fold_rmsd_A"] <= 2.5 for r in subset),
                "ranked_samples": [
                    {"seed": seed, "sample": best["sample"], "ipsae_min": best["scores"]["ipsae_min"],
                     "above_negative_max": best["scores"]["ipsae_min"] > threshold,
                     "interface_recovered": best["interface_recovered"],
                     "binder_fold_rmsd_A": best["pose"]["binder_CA_fold_rmsd_A"]}
                    for seed in manifest["refolding"]["seeds"]
                    for best in [max((r for r in subset if r["seed"] == seed),
                                     key=lambda r: r["scores"]["ipsae_min"])]],
            }
        thresholds[receptor] = {"negative_max": threshold, "positive_checks": checks,
                               "warning": "Descriptive threshold on two presumed negatives, not validated specificity"}
    result["provisional_threshold_checks"] = thresholds
    result["generation_gate_passed"] = any(
        all(all(x["above_negative_max"] and x["interface_recovered"]
                and x["binder_fold_rmsd_A"] <= 2.5 for x in check["ranked_samples"])
            for check in receptor["positive_checks"].values())
        for receptor in thresholds.values())
    (args.output / "results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(thresholds, indent=2))
    print("Generation gate passed:", result["generation_gate_passed"])


if __name__ == "__main__":
    main()
