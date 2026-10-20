"""Compare sequence-only Chai controls with experimental interfaces and folds."""
import argparse
import json
from pathlib import Path

import numpy as np

from analyze_controls import pose_metrics
from prepare_controls import ROOT, geometry, make_complex, read_structure, residues, sequence
from score_matched_controls import interface_pae, run_ipsae


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("--output", type=Path, default=ROOT/"reports/controls_chai")
    args = p.parse_args()
    manifest = json.loads((args.run/"inputs/manifest.json").read_text())
    cases = {c["binder_id"]: c for c in manifest["cases"] if c["receptor_id"] == "6wrv"}
    full = read_structure(ROOT/"data/raw/11QC.cif")
    rows = []
    for binder, case in cases.items():
        folder = args.run/binder
        files = json.loads((folder/"candidate_order.json").read_text())
        arrays = np.load(folder/"confidence_arrays.npz", allow_pickle=False)
        assert len(files) == len(arrays["pae"]) == 5
        for i, name in enumerate(files):
            raw = read_structure(folder/name)
            chains = [c for c in raw[0] if residues(c)]
            # PAE order follows FASTA order; require the output to agree exactly.
            assert len(chains) == 2
            assert sequence(chains[0]) == case["target_sequence"]
            assert sequence(chains[1]) == case["binder_sequence"]
            s = make_complex(residues(chains[0]), residues(chains[1]))
            pae, plddt = arrays["pae"][i], arrays["plddt"][i]
            size = len(case["target_sequence"]) + len(case["binder_sequence"])
            assert pae.shape == (size, size) and plddt.shape == (size,)
            assert np.isfinite(pae).all() and np.isfinite(plddt).all()
            dest = ROOT/"outputs/chai_scoring"/args.run.name/binder/str(i)
            scores = run_ipsae(s, pae, plddt, dest)
            scores.update(interface_pae(s, pae))
            n = len(s[0]["A"])
            scores["binder_mean_plddt"] = float(plddt[n:].mean())
            with np.load(folder/f"scores.model_idx_{i}.npz", allow_pickle=False) as native_scores:
                native = {k: native_scores[k].tolist() for k in native_scores.files}
            row = {"binder_id": binder, "kind": case["kind"], "sample": i, "seed": 42,
                   "scores": scores, "chai_scores": native,
                   "geometry": geometry(s, {int(k): v for k, v in case["target_native_map"].items()}, full)}
            row["full_context_fit_valid"] = row["geometry"]["target_fit_rmsd_A"] <= 2
            if case["kind"] == "positive":
                reference = read_structure(ROOT/next(
                    c["pdb"] for c in manifest["cases"] if c["id"] == case["reference_id"]))
                row["pose"] = pose_metrics(s, reference)
                row["interface_recovered"] = (
                    row["pose"]["native_contact_recall"] >= .5 and
                    row["pose"]["binder_interface_CA_rmsd_after_target_fit_A"] <= 2)
            rows.append(row)
    result = {
        "run": args.run.name, "settings": json.loads((args.run/"settings.json").read_text()),
        "limitation": "One seed; different predictor and template conditioning. Published structures may occur in training. Presumed negatives are not measured nonbinders.",
        "rows": rows,
        "summary": {
            binder: {"samples": 5,
                     "interface_recovered": sum(r.get("interface_recovered", False)
                                                for r in rows if r["binder_id"] == binder),
                     "ipsae_min_range": [min(r["scores"]["ipsae_min"] for r in rows if r["binder_id"] == binder),
                                          max(r["scores"]["ipsae_min"] for r in rows if r["binder_id"] == binder)]}
            for binder in cases},
    }
    args.output.mkdir(exist_ok=True, parents=True)
    (args.output/"results.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
