"""Sequence-only Chai-1 cross-check: four controls, one seed, five samples each."""
import argparse
import gc
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    import numpy as np
    import torch
    from chai_lab.chai1 import run_inference

    manifest = json.loads((args.output/"inputs/manifest.json").read_text())
    cases = [c for c in manifest["cases"] if c["receptor_id"] == "6wrv"]
    assert len(cases) == 4
    settings = {"num_trunk_recycles": 3, "num_diffn_timesteps": 200,
                "num_diffn_samples": 5, "num_trunk_samples": 1, "seed": 42,
                "use_esm_embeddings": True, "use_msa_server": False,
                "use_templates_server": False, "device": "cuda:0", "low_memory": True}
    (args.output/"settings.json").write_text(json.dumps({
        "chai_lab_version": "0.6.1", "settings": settings,
        "input": "Sequences only; no receptor/binder template or contact constraints",
        "comparison_limit": "Different predictor AND conditioning from target-templated BoltzGen",
        "cases": [c["id"] for c in cases]}, indent=2))
    for case in cases:
        name = case["binder_id"]
        fasta = args.output/"inputs"/f"{name}.fasta"
        fasta.write_text(f">protein|name=target\n{case['target_sequence']}\n"
                         f">protein|name=binder\n{case['binder_sequence']}\n")
        dest = args.output/name
        print(f"Starting {name}", flush=True)
        result = run_inference(fasta, output_dir=dest, **settings)
        size = len(case["target_sequence"]) + len(case["binder_sequence"])
        assert tuple(result.pae.shape) == (5, size, size)
        assert tuple(result.plddt.shape) == (5, size)
        np.savez_compressed(dest/"confidence_arrays.npz",
                            pae=result.pae.cpu().numpy(), plddt=result.plddt.cpu().numpy())
        (dest/"candidate_order.json").write_text(json.dumps(
            [p.name for p in result.cif_paths], indent=2))
        del result
        gc.collect()
        torch.cuda.empty_cache()
        print(f"Finished {name}", flush=True)


if __name__ == "__main__":
    main()
