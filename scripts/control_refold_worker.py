"""Run the pinned BoltzGen folding task directly, with no binder templates."""
import argparse
import json
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("design_dir", type=Path)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--target-msa", type=Path)
    args = p.parse_args()
    import hydra
    import numpy as np
    from omegaconf import OmegaConf
    from pytorch_lightning import seed_everything
    from boltzgen.cli.boltzgen import ARTIFACTS, get_artifact_path, config_dir

    seed_everything(args.seed, workers=True)
    if args.target_msa:
        # Explicit diagnostic adapter: inject target MSA at Input construction.
        # This is separate from the unchanged upstream single-sequence baseline.
        from boltzgen.task.predict import data_from_generated as generated
        from boltzgen.data.parse.a3m import process_a3m
        target_msa = process_a3m(args.target_msa, max_seqs=512)
        assert len(target_msa.sequences) > 1
        original_input = generated.Input

        def with_target_msa(**kwargs):
            assert kwargs["msa"] == {}
            assert len(kwargs["structure"].chains) == 2
            kwargs["msa"] = {0: target_msa}
            return original_input(**kwargs)

        generated.Input = with_target_msa
    artifact_args = argparse.Namespace(cache=None, models_token=None, force_download=False)
    checkpoint = get_artifact_path(artifact_args, ARTIFACTS["folding"][0])
    moldir = get_artifact_path(artifact_args, ARTIFACTS["moldir"][0], repo_type="dataset")
    cfg = OmegaConf.load(config_dir / "fold.yaml")
    cfg = OmegaConf.merge(cfg, OmegaConf.from_dotlist([
        f"output={args.design_dir}", f"data.design_dir={args.design_dir}",
        f"checkpoint={checkpoint}", f"data.cfg.moldir={moldir}",
        "data.cfg.num_workers=0", "data.target_templates=true",
        "trainer.devices=1", "override.use_kernels=true",
        "recycling_steps=3", "sampling_steps=200", "diffusion_samples=5",
    ]))
    if args.target_msa:
        cfg.data.cfg.max_seqs = 512
    OmegaConf.save(cfg, args.design_dir / "fold_config.yaml")
    task = hydra.utils.instantiate(cfg)
    expected = json.loads((args.design_dir / "expected.json").read_text())
    checks = []
    # Fail explicitly before prediction if a residue mask or chain role is wrong.
    for path in task.data.predict_set.generated_paths:
        features = task.data.predict_set.get_sample(args.design_dir, path.stem)
        info = expected[path.stem]
        n_target = len(info["target_sequence"])
        n_binder = len(info["binder_sequence"])
        mask = features["chain_design_mask"].cpu().numpy().astype(bool)
        assert len(mask) == n_target+n_binder, (path.stem, len(mask), n_target+n_binder)
        assert not mask[:n_target].any() and mask[n_target:].all(), path.stem
        assert not features["exception"], path.stem
        assert not features["template_mask"][:, n_target:].any(), path.stem
        assert features["template_mask"][:, :n_target].all(), path.stem
        msa_depth = int(features["msa"].shape[0])
        if args.target_msa:
            assert msa_depth > 1, "Target MSA was silently dropped"
        template_info = {}
        for key, value in features.items():
            if "template" in key and hasattr(value, "shape"):
                template_info[key] = list(value.shape)
        checks.append({"id": path.stem, "target_tokens": n_target, "binder_tokens": n_binder,
                       "design_mask_verified": True, "binder_template_excluded": True,
                       "msa_depth": msa_depth, "template_features": template_info})
    assert len(checks) == len(expected)
    (args.design_dir / "preflight.json").write_text(json.dumps(checks, indent=2))
    task.run(cfg)
    from boltzgen.data import const
    summary = []
    for case in expected:
        path = args.design_dir / const.folding_dirname / f"{case}.npz"
        data = np.load(path, allow_pickle=False)
        metrics = {k: data[k].tolist() for k in [
            "iptm", "ptm", "design_to_target_iptm", "target_ptm",
            "design_ptm", "min_design_to_target_pae", "interaction_pae"] if k in data}
        metrics["selected_sample"] = int(np.argmax(.8*data["iptm"]+.2*data["ptm"]))
        summary.append({"id": case, "seed": args.seed, **metrics})
    (args.design_dir / "confidence_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
