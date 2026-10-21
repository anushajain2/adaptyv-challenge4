"""Pinned generation preflight, seeded stages, and matched refolding assembly."""
import argparse
import json
import subprocess
from pathlib import Path


def config_path(folder, stage):
    import yaml
    steps = yaml.safe_load((folder/"steps.yaml").read_text())["steps"]
    return folder/next(s["config_file"] for s in steps if s["name"] == stage)


def preflight(out):
    import hydra
    import numpy as np
    from omegaconf import OmegaConf
    from boltzgen.data import const
    source = out/"inputs"
    manifest = json.loads((source/"manifest.json").read_text())
    checks = []
    for arm in manifest["arms"]:
        name = arm["id"]
        dest = out/name
        subprocess.run(["boltzgen", "configure", str(source/f"{name}.yaml"),
                        "--output", str(dest), "--protocol", "protein-anything",
                        "--num_designs", "10", "--budget", "10", "--devices", "1",
                        "--num_workers", "0", "--steps", "design", "inverse_folding"], check=True)
        cfg = OmegaConf.load(config_path(dest, "design"))
        data = hydra.utils.instantiate(cfg.data)
        feat = data.predict_set.get_sample(source/f"{name}.yaml")
        mask = feat["design_mask"].numpy().astype(bool)
        assert len(mask) == 282 and not mask[:194].any() and mask[194:].all()
        assert not feat.get("exception", False)
        ss = feat["ss_type"].numpy()
        positions = np.flatnonzero(ss[194:] == const.ss_type_ids["SHEET"])+1
        expected = list(range(33,38)) if name == "sheet" else []
        assert positions.tolist() == expected, (name, positions)
        bind = feat["binding_type"].numpy()
        hotspots = np.flatnonzero(bind[:194] == const.binding_type_ids["BINDING"])+1
        assert hotspots.tolist() == manifest["hotspots_crop"], hotspots
        assert len(data.predict_set) == 10
        checks.append({"arm":name,"design_tokens":int(mask.sum()),
                       "target_tokens":194,"sheet_positions":positions.tolist(),
                       "hotspots_crop":hotspots.tolist(),"design_count":len(data.predict_set)})
    (out/"generation_preflight.json").write_text(json.dumps(checks, indent=2))


def run_stage(out, arm, stage):
    import hydra
    from omegaconf import OmegaConf
    from pytorch_lightning import seed_everything
    manifest = json.loads((out/"inputs/manifest.json").read_text())
    settings = next(a for a in manifest["arms"] if a["id"] == arm)
    seed = settings["generation_seed" if stage == "design" else "inverse_fold_seed"]
    seed_everything(seed, workers=True)
    cfg = OmegaConf.load(config_path(out/arm, stage))
    hydra.utils.instantiate(cfg).run(cfg)


def assemble(out):
    import numpy as np
    from prepare_controls import make_complex, read_structure, residues, sequence, save_input_cif, save_pdb, fit_mapped
    source = out/"inputs"
    manifest = json.loads((source/"manifest.json").read_text())
    target = read_structure(source/"target_6wrw.cif")
    refs = out/"references"
    refs.mkdir()
    cases = []
    for arm in manifest["arms"]:
        folder = out/arm["id"]/"intermediate_designs_inverse_folded"
        files = sorted(p for p in folder.glob("*.cif") if not p.stem.endswith("_native"))
        assert len(files) == 10, (arm["id"], len(files))
        for index, path in enumerate(files):
            raw = read_structure(path)
            chain_sequences = {c.id:sequence(c) for c in raw[0]}
            targets = [c for c,s in chain_sequences.items() if s == manifest["target_sequence"]]
            assert len(targets) == 1 and len(chain_sequences) == 2
            binder_id = next(c for c in chain_sequences if c != targets[0])
            normalized = make_complex(residues(raw[0][targets[0]]), residues(raw[0][binder_id]))
            assert len(normalized[0]["B"]) == 88
            fitted, rms = fit_mapped(normalized, {i:i for i in range(1,195)}, target)
            # Fix the same original receptor conformation for every refolding.
            ref = make_complex(residues(target[0]["A"]), residues(fitted[0]["B"]))
            name = f"{arm['id']}_{index:02d}"
            binder_seq = sequence(ref[0]["B"])
            assert set(binder_seq) <= set("ACDEFGHIKLMNPQRSTVWY")
            save_pdb(ref, refs/f"{name}.pdb")
            save_input_cif(ref, refs/f"{name}.cif", manifest["target_sequence"], binder_seq)
            np.savez_compressed(refs/f"{name}.npz", design_mask=np.array([False]*194+[True]*88))
            cases.append({"id":name,"arm":arm["id"],"kind":"novel_design",
                          "target_sequence":manifest["target_sequence"],"binder_sequence":binder_seq,
                          "target_native_map":manifest["target_native_map"],
                          "generated_source":str(path.relative_to(out)),
                          "generation_target_fit_rmsd_A":rms})
    for control in manifest["controls"]:
        case = dict(control, arm="control")
        cases.append(case)
        for ext in (".cif", ".npz", ".pdb"):
            (refs/(case["id"]+ext)).write_bytes((source/(case["id"]+ext)).read_bytes())
    (out/"cases.json").write_text(json.dumps(cases, indent=2))
    for seed in manifest["refolding"]["seeds"]:
        dest = out/f"seed-{seed}"
        dest.mkdir()
        for case in cases:
            for ext in (".cif", ".npz"):
                (dest/(case["id"]+ext)).write_bytes((refs/(case["id"]+ext)).read_bytes())
        (dest/"expected.json").write_text(json.dumps({c["id"]:c for c in cases}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("operation", choices=["preflight", "stage", "assemble"])
    p.add_argument("output", type=Path)
    p.add_argument("--arm", choices=["unconstrained", "sheet"])
    p.add_argument("--stage", choices=["design", "inverse_folding"])
    args = p.parse_args()
    if args.operation == "preflight": preflight(args.output)
    elif args.operation == "assemble": assemble(args.output)
    else: run_stage(args.output, args.arm, args.stage)
