"""Prepare a matched two-receptor panel and freeze its scoring plan before inference."""
import copy
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
from scipy.spatial import cKDTree

from prepare_controls import (
    ROOT, fit_mapped, make_complex, read_structure, residues, save_input_cif,
    save_pdb, sequence,
)


def main():
    out = ROOT / "data/controls_matched"
    out.mkdir(exist_ok=True)
    old = json.loads((ROOT / "data/controls/manifest.json").read_text())
    old_cases = {c["id"]: c for c in old["cases"]}
    original = read_structure(ROOT / old_cases["3ds18_crystal"]["pdb"])
    raw = read_structure(ROOT / "data/raw/6WRW.cif")
    target = [raw[0]["A"][i] for i in range(190, 384)]
    tree = cKDTree(np.array([a.coord for r in target for a in r]))
    counts = {c.id: int(tree.query_ball_point(
        np.array([a.coord for r in residues(c) for a in r]), 5, return_length=True).sum())
        for c in raw[0] if len(residues(c)) == 88}
    binder = max(counts, key=counts.get)
    second = make_complex(target, residues(raw[0][binder]))
    assert sequence(second[0]["A"]) == sequence(original[0]["A"])
    negatives = {
        "scrambled": read_structure(ROOT / old_cases["3ds18_scrambled"]["pdb"]),
        "caplacizumab": read_structure(ROOT / old_cases["caplacizumab_current"]["pdb"]),
    }
    receptors = {"6wrv": original, "6wrw": second}
    binders = {"3ds18": original, "2ds25_5": second, **negatives}
    records = []
    for receptor_id, receptor in receptors.items():
        for binder_id, source in binders.items():
            # Only target coordinates are used as templates. Placement here just
            # preserves the experimentally defined reference for the positives.
            fitted, _ = fit_mapped(source, {i: i for i in range(1, 195)}, receptor)
            s = make_complex(residues(receptor[0]["A"]), residues(fitted[0]["B"]))
            name = f"{binder_id}_{receptor_id}"
            path = out / f"{name}.pdb"
            save_pdb(s, path)
            t, b = sequence(s[0]["A"]), sequence(s[0]["B"])
            save_input_cif(s, path.with_suffix(".cif"), t, b)
            np.savez_compressed(path.with_suffix(".npz"),
                                design_mask=np.array([False]*len(t) + [True]*len(b)))
            positive = binder_id in ("3ds18", "2ds25_5")
            records.append({
                "id": name, "kind": "positive" if positive else "negative_hypothesis",
                "binder_id": binder_id, "receptor_id": receptor_id,
                "role": "held_out_positive" if binder_id == "2ds25_5" else "calibration",
                "reference_id": (f"{binder_id}_{'6wrv' if binder_id == '3ds18' else '6wrw'}"
                                 if positive else name),
                "target_sequence": t, "binder_sequence": b,
                "target_native_map": {i: i+189 for i in range(1, 195)},
                "target_chain": "A", "binder_chain": "B",
                "pdb": str(path.relative_to(ROOT)),
                "cif": str(path.with_suffix(".cif").relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "note": "Published positive or presumed negative; never a novel candidate",
            })
    manifest = {
        "scope": "Matched calibration, not competition candidates",
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "max_prediction_seconds": 1200,
        "source_6wrw": {"url": "https://files.rcsb.org/download/6WRW.cif",
                       "sha256": hashlib.sha256((ROOT/"data/raw/6WRW.cif").read_bytes()).hexdigest(),
                       "receptor_chain": "A", "binder_chain": binder,
                       "binder_chain_contact_counts": counts},
        "refolding": {**old["refolding"], "export_pae": True},
        "analysis_plan": {
            "primary_ipsae_cutoffs_A": {"pae": 10, "distance": 10},
            "ipsae_statistic": "Minimum of the two directional ipSAE values; also report maximum",
            "interface_pae": "Both directions over residue pairs with heavy-atom distance <5 A; report contact count",
            "recovery_definition": "Native contact recall >=0.5 AND receptor-fitted interface CA RMSD <=2 A",
            "whole_binder_check": "Report binder self-fit CA RMSD; <=2.5 A is a separate provisional fold criterion",
            "decision": "Proceed only if a ranking rule learned on 3DS18 separates matched negatives, recovers interfaces on both seeds, and transfers to held-out 2DS25.5 without changing cutoffs. Whole-binder fold failures block claiming parent qualification.",
            "uncertainty": "Two presumed negatives and published positives do not estimate prospective hit rate. No post-hoc lowering of structural criteria.",
            "compute_planning_ceiling_usd": 5,
        },
        "cases": records,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"cases": len(records), "samples": len(records)*10,
                      "source": manifest["source_6wrw"]}, indent=2))


if __name__ == "__main__":
    main()
