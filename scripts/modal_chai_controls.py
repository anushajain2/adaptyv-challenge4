"""Bounded Chai-1 control check adapted from the Adaptyv/biomodals wrapper."""
import json
import os
import sys
from pathlib import Path

import modal

ROOT = Path(__file__).resolve().parents[1]
if modal.is_local() and not any(x in sys.argv for x in ("--help", "-h")):
    if os.environ.get("TFR1_LAUNCH_RESERVED") != "1":
        raise RuntimeError("Use launch_controls.py --engine chai")

image = modal.Image.debian_slim(python_version="3.12").uv_pip_install(
    "chai_lab==0.6.1", "torch==2.5.1")
volume = modal.Volume.from_name("tfr1-pilot-results", create_if_missing=False)
app = modal.App("tfr1-chai-control-check", image=image)


@app.function(gpu="L40S", cpu=(2, 4), memory=(16384, 32768),
              timeout=1350, startup_timeout=600, retries=0,
              max_containers=1, scaledown_window=2, volumes={"/results": volume})
def run_controls(files: dict[str, bytes], run_id: str):
    import hashlib
    import signal
    import subprocess
    import time
    import traceback
    out = Path("/results")/run_id
    out.mkdir(exist_ok=False)
    source = out/"inputs"
    source.mkdir()
    for name, data in files.items():
        if Path(name).name != name:
            raise ValueError("Flat filenames only")
        (source/name).write_bytes(data)
    (out/"input_hashes.json").write_text(json.dumps(
        {k: hashlib.sha256(v).hexdigest() for k, v in files.items()}, indent=2))
    (out/"requirements-cloud.txt").write_text(subprocess.check_output(
        ["python", "-m", "pip", "freeze"], text=True))
    started = time.monotonic()
    result = {"run_id": run_id, "completed": False, "timeout_seconds": 1200}
    try:
        with (out/"worker.log").open("w") as log:
            proc = subprocess.Popen(
                ["python", str(source/"chai_control_worker.py"), str(out)],
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                env=dict(os.environ, PYTHONUNBUFFERED="1"))
            last = ""
            while proc.poll() is None:
                if time.monotonic()-started > 1200:
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        proc.wait()
                    raise TimeoutError("Chai control run exceeded 20-minute limit")
                try:
                    proc.wait(timeout=25)
                except subprocess.TimeoutExpired:
                    pass
                tail = "\n".join((out/"worker.log").read_text(errors="replace").splitlines()[-3:])
                if tail != last:
                    print(tail, flush=True)
                    last = tail
                volume.commit()
        result["returncode"] = proc.returncode
        if proc.returncode:
            raise RuntimeError("Chai control worker failed; no automatic retry")
        result["completed"] = True
    except Exception:
        (out/"error.txt").write_text(traceback.format_exc())
        raise
    finally:
        result["elapsed_seconds"] = time.monotonic()-started
        (out/"run_status.json").write_text(json.dumps(result, indent=2))
        volume.commit()
    return result


@app.local_entrypoint()
def main(run_id: str, input_dir: str):
    reservation = json.loads((ROOT/"outputs/controls_launch.json").read_text())
    if (reservation["run_id"] != run_id or reservation["engine"] != "chai"
            or reservation["input_dir"] != input_dir):
        raise ValueError("Run differs from reserved inputs")
    files = {"manifest.json": (ROOT/input_dir/"manifest.json").read_bytes(),
             "chai_control_worker.py": (ROOT/"scripts/chai_control_worker.py").read_bytes()}
    result = run_controls.remote(files, run_id)
    print(json.dumps(result, indent=2))
