"""Bounded prospective target-integrity mutant_validation with the cached pinned image."""
import json
import os
from pathlib import Path

import modal
import sys

ROOT = Path(__file__).resolve().parents[1]
PIN = "247b9bbd8b68a60aba854c2968d6a0cddd21ad6d"
if modal.is_local() and not any(x in sys.argv for x in ("--help", "-h")):
    if os.environ.get("TFR1_LAUNCH_RESERVED") != "1":
        raise RuntimeError("Use launch_expansion_his.py to reserve the budget before image creation.")
# Same image layers as modal_pilot.py; standalone so remote import needs no local modules.
image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git", "build-essential")
    .uv_pip_install("torch==2.5.1")
    .run_commands(
        "git clone https://github.com/HannesStark/boltzgen /opt/boltzgen",
        f"cd /opt/boltzgen && git checkout {PIN} && pip install -e .",
    )
    .run_commands("boltzgen download all")
    # Upstream A3M parser imports Redis for optional taxonomy lookup even when unused.
    .uv_pip_install("redis==5.2.1")
)
volume = modal.Volume.from_name("tfr1-pilot-results", create_if_missing=True)

app = modal.App("tfr1-expansion-his-confirmation", image=image)


@app.function(gpu="L40S", cpu=(2, 4), memory=(16384, 32768),
              timeout=2250, startup_timeout=120, retries=0,
              max_containers=1, scaledown_window=2, volumes={"/results": volume})
def run_controls(files: dict[str, bytes], run_id: str):
    import hashlib
    import signal
    import subprocess
    import time
    import traceback
    out = Path("/results") / run_id
    out.mkdir(parents=True, exist_ok=False)
    source = out / "inputs"
    source.mkdir()
    for name, data in files.items():
        if Path(name).name != name:
            raise ValueError("Flat input filenames required")
        (source / name).write_bytes(data)
    manifest = json.loads((source / "manifest.json").read_text())
    max_seconds = int(manifest.get("max_prediction_seconds", 2100))
    if not 60 <= max_seconds <= 2100:
        raise ValueError("Invalid prediction time limit")
    (out / "input_hashes.json").write_text(json.dumps(
        {k: hashlib.sha256(v).hexdigest() for k, v in files.items()}, indent=2))
    (out / "requirements-cloud.txt").write_text(subprocess.check_output(
        ["python", "-m", "pip", "freeze"], text=True))
    started = time.monotonic()
    statuses = []
    try:
        for seed in manifest["refolding"]["seeds"]:
            dest = out / f"seed-{seed}"
            dest.mkdir()
            for case in manifest["cases"]:
                name = case["id"]
                for suffix in (".cif", ".npz"):
                    (dest / (name+suffix)).write_bytes(files[name+suffix])
            (dest / "expected.json").write_text(json.dumps(
                {c["id"]: c for c in manifest["cases"]}, indent=2))
            cmd = ["python", str(source / "control_refold_worker.py"),
                   str(dest), "--seed", str(seed), "--samples", str(manifest["refolding"]["diffusion_samples_per_seed"])]
            if manifest["refolding"].get("target_msa"):
                cmd += ["--target-msa", str(source / manifest["refolding"]["target_msa"])]
            if manifest["refolding"].get("export_pae"):
                cmd += ["--export-pae"]
            log_path = out / f"seed-{seed}.log"
            with log_path.open("w") as log:
                proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                        env=dict(os.environ, PYTHONUNBUFFERED="1"),
                                        start_new_session=True)
                last = ""
                while proc.poll() is None:
                    if time.monotonic()-started > max_seconds:
                        os.killpg(proc.pid, signal.SIGTERM)
                        try:
                            proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            os.killpg(proc.pid, signal.SIGKILL)
                            proc.wait()
                        raise TimeoutError(f"Control run exceeded {max_seconds}-second compute bound")
                    try:
                        proc.wait(timeout=25)
                    except subprocess.TimeoutExpired:
                        pass
                    tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-3:])
                    if tail != last:
                        print(f"[seed {seed}] {tail}", flush=True)
                        last = tail
                    volume.commit()
            statuses.append({"seed": seed, "returncode": proc.returncode})
            if proc.returncode:
                raise RuntimeError(f"Control prediction failed for seed {seed}")
    except Exception:
        (out / "error.txt").write_text(traceback.format_exc())
        raise
    finally:
        result = {"run_id": run_id, "boltzgen_commit": PIN, "statuses": statuses,
                  "elapsed_seconds": time.monotonic()-started}
        (out / "run_status.json").write_text(json.dumps(result, indent=2))
        volume.commit()
    return result


@app.local_entrypoint()
def main(run_id: str = "his-expansion-confirm-021-6wrw", input_dir: str = "data/expansion_his_confirmation"):
    import re
    if not re.fullmatch(r"[A-Za-z0-9_-]+", run_id):
        raise ValueError("Invalid run id")
    reservation = json.loads((ROOT / "outputs" / f"{run_id}_launch.json").read_text())
    if reservation["run_id"] != run_id:
        raise ValueError("Run differs from budget reservation")
    path = (ROOT / input_dir).resolve()
    if not path.is_relative_to(ROOT.resolve()) or str(path.relative_to(ROOT)) != reservation["input_dir"]:
        raise ValueError("Input directory differs from reservation or lies outside workspace")
    files = {p.name: p.read_bytes() for p in path.iterdir()
             if p.suffix in {".cif", ".npz", ".json", ".a3m"}}
    worker = ROOT / "scripts/control_refold_worker.py"
    files[worker.name] = worker.read_bytes()
    import hashlib
    assert {k:hashlib.sha256(v).hexdigest() for k,v in files.items()} == reservation["input_sha256"]
    result = run_controls.remote(files, run_id)
    (ROOT / "outputs" / f"{run_id}_result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
