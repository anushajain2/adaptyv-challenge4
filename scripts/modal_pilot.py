"""Bounded Modal pilot adapted from hgbrian/biomodals/modal_boltzgen.py.

Changes: CPU-only build, exact upstream pin, complete file bundling, serial
protocols, persistent outputs, no retries, runtime bounds, local launch ledger.
"""
from pathlib import Path
import json
import os
import modal

# Guard happens during module import, BEFORE remote image builds can incur costs.
# --help stays available without a launch reservation.
import sys
if not any(x in sys.argv for x in ("--help", "-h")) and modal.is_local():
    if os.environ.get("TFR1_LAUNCH_RESERVED") != "1":
        raise RuntimeError("Use scripts/launch_pilot.py to reserve the budget before building the image.")

PIN = "247b9bbd8b68a60aba854c2968d6a0cddd21ad6d"
ROOT = Path(__file__).resolve().parents[1]
image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git", "build-essential")
    .uv_pip_install("torch==2.5.1")
    .run_commands(
        "git clone https://github.com/HannesStark/boltzgen /opt/boltzgen",
        f"cd /opt/boltzgen && git checkout {PIN} && pip install -e .",
    )
    .run_commands("boltzgen download all")
)
app = modal.App("tfr1-day1-3-pilot", image=image)
volume = modal.Volume.from_name("tfr1-pilot-results", create_if_missing=True)


@app.function(gpu="L40S", cpu=(2, 4), memory=(16384, 32768),
              timeout=14400, startup_timeout=1800, retries=0,
              max_containers=1, scaledown_window=2, volumes={"/results": volume})
def run_pilot(files: dict[str, bytes], hotspot: str, run_id: str, protocols: list[str]):
    import subprocess
    import time
    import traceback
    import signal
    import hashlib
    in_dir = Path("/tmp/tfr1-input")
    in_dir.mkdir(exist_ok=True)
    for name, content in files.items():
        if Path(name).name != name:
            raise ValueError("Only flat input paths allowed")
        (in_dir / name).write_bytes(content)
    out = Path("/results") / run_id
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    statuses = []
    (out / "input_hashes.json").write_text(json.dumps(
        {name:hashlib.sha256(data).hexdigest() for name,data in files.items()},indent=2))
    (out / "requirements-cloud.txt").write_text(
        subprocess.check_output(["python","-m","pip","freeze"],text=True))

    def logged_run(cmd,log_path,timeout):
        """Checkpoint progress and preserve partial outputs during long stages."""
        deadline=time.monotonic()+timeout
        with log_path.open("w") as log:
            proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,
                                  env=dict(os.environ,PYTHONUNBUFFERED="1"),
                                  start_new_session=True)
            last_tail=""
            while proc.poll() is None:
                remaining=deadline-time.monotonic()
                if remaining<=0:
                    os.killpg(proc.pid,signal.SIGTERM)
                    try:
                        proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid,signal.SIGKILL)
                        proc.wait()
                    raise TimeoutError(f"Stage timeout: {log_path.name}")
                try:
                    proc.wait(timeout=min(30,remaining))
                except subprocess.TimeoutExpired:
                    pass
                tail="\n".join(log_path.read_text(errors="replace").splitlines()[-3:])
                if tail!=last_tail:
                    print(f"[{log_path.name}]\n{tail}",flush=True)
                    last_tail=tail
                volume.commit()
            return proc.returncode
    try:
        # Validate BOTH modalities before launching either generation stage.
        for protocol in protocols:
            spec=in_dir/f"{hotspot}_{protocol}.yaml"
            check=logged_run(["boltzgen","check",str(spec)],
                             out/f"validation-{protocol}.log",300)
            if check:
                raise RuntimeError(f"Specification validation failed: {protocol}")
        # Budget is TOTAL designs, not the final number surviving filters.
        for protocol in protocols:
            spec = in_dir / f"{hotspot}_{protocol}.yaml"
            cmd = ["boltzgen", "run", str(spec), "--output", str(out / protocol),
                   "--protocol", protocol, "--num_designs", "25",
                   "--budget", "10", "--devices", "1"]
            remaining = max(1, 13800 - (time.monotonic() - start))
            returncode=logged_run(cmd,out/f"{protocol}.log",remaining)
            statuses.append({"protocol": protocol, "returncode": returncode,
                             "requested_designs": 25})
            volume.commit()
            if returncode:
                raise RuntimeError(f"BoltzGen failed: {protocol}")
    except Exception:
        (out / "error.txt").write_text(traceback.format_exc())
        raise
    finally:
        (out / "run_status.json").write_text(json.dumps({
            "statuses": statuses, "elapsed_seconds": time.monotonic()-start,
            "boltzgen_commit": PIN, "hotspot": hotspot,
        }, indent=2))
        volume.commit()
    # Caller downloads through Modal Volume CLI: avoids large RPC payloads.
    return {"volume": "tfr1-pilot-results", "run_id": run_id, "statuses": statuses}


@app.local_entrypoint()
def main(hotspot: str = "conserved_326_347", run_id: str = "pilot-001",
         protocol: str = "all"):
    import re
    import yaml
    if not re.fullmatch(r"[a-zA-Z0-9_-]+", run_id):
        raise ValueError("Unsafe run id")
    if not re.fullmatch(r"[a-zA-Z0-9_-]+", hotspot):
        raise ValueError("Unsafe hotspot name")
    if protocol not in {"all","protein-anything","nanobody-anything"}:
        raise ValueError("Invalid protocol")
    protocols=["protein-anything","nanobody-anything"] if protocol=="all" else [protocol]
    files = {}
    for protocol in protocols:
        path = ROOT / "configs" / f"{hotspot}_{protocol}.yaml"
        spec = yaml.safe_load(path.read_text())
        files[path.name] = path.read_bytes()
        for entity in spec["entities"]:
            if "file" in entity:
                name = entity["file"]["path"]
                if not isinstance(name, str) or Path(name).name != name:
                    raise ValueError("Use flattened spec paths")
                files[name] = (path.parent / name).read_bytes()
    # A prior run, including a failed attempt, needs cost reconciliation before another.
    ledger = ROOT / "outputs" / "cloud_launch.json"
    ledger.parent.mkdir(exist_ok=True)
    reservation = json.loads(ledger.read_text())
    if reservation["run_id"] != run_id:
        raise ValueError("Run id differs from budget reservation")
    result = run_pilot.remote(files, hotspot, run_id, protocols)
    (ledger.parent / "cloud_result.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"Download: .venv/bin/python scripts/download_pilot.py {run_id} outputs/pilot")
