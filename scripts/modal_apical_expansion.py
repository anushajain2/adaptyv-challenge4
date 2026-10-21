"""Bounded 200-design campaign: one GPU, no retries, reserved before launch."""
import json
import os
import sys
from pathlib import Path

import modal

ROOT = Path(__file__).resolve().parents[1]
PIN = "247b9bbd8b68a60aba854c2968d6a0cddd21ad6d"
if modal.is_local() and not any(x in sys.argv for x in ("--help", "-h")):
    if os.environ.get("TFR1_LAUNCH_RESERVED") != "1":
        raise RuntimeError("Use launch_next_round.py before any image build.")
image = (modal.Image.debian_slim(python_version="3.12")
         .apt_install("git", "build-essential")
         .uv_pip_install("torch==2.5.1")
         .run_commands("git clone https://github.com/HannesStark/boltzgen /opt/boltzgen",
                       f"cd /opt/boltzgen && git checkout {PIN} && pip install -e .")
         .run_commands("boltzgen download all")
         .uv_pip_install("redis==5.2.1"))
volume = modal.Volume.from_name("tfr1-pilot-results", create_if_missing=True)
app = modal.App("tfr1-next-expansion", image=image)


@app.function(gpu="L40S", cpu=(2,4), memory=(16384,32768), timeout=6750,
              startup_timeout=600, retries=0, max_containers=1,
              scaledown_window=2, volumes={"/results":volume})
def run_comparison(files: dict[str,bytes], run_id: str):
    import hashlib
    import signal
    import subprocess
    import time
    import traceback
    started = time.monotonic()
    out = Path("/results")/run_id
    out.mkdir(exist_ok=False)
    source = out/"inputs"
    source.mkdir()
    for name, data in files.items():
        if Path(name).name != name: raise ValueError("Flat filenames required")
        (source/name).write_bytes(data)
    manifest = json.loads((source/"manifest.json").read_text())
    assert manifest["max_prediction_seconds"] == 6600
    statuses = []
    (out/"input_hashes.json").write_text(json.dumps(
        {k:hashlib.sha256(v).hexdigest() for k,v in files.items()}, indent=2))
    (out/"requirements-cloud.txt").write_text(subprocess.check_output(
        ["python","-m","pip","freeze"], text=True))

    def run(label, cmd):
        log_path = out/f"{label}.log"
        with log_path.open("w") as log:
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                    cwd=source, start_new_session=True,
                                    env=dict(os.environ,PYTHONUNBUFFERED="1"))
            last = ""
            while proc.poll() is None:
                if time.monotonic()-started >= 6600:
                    os.killpg(proc.pid,signal.SIGTERM)
                    try: proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid,signal.SIGKILL)
                        proc.wait()
                    raise TimeoutError("Experiment reached its 6600-second work limit")
                try: proc.wait(timeout=25)
                except subprocess.TimeoutExpired: pass
                tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-3:])
                if tail != last:
                    print(f"[{label}] {tail}", flush=True)
                    last = tail
                volume.commit()
        statuses.append({"stage":label,"returncode":proc.returncode,
                         "elapsed_total_seconds":time.monotonic()-started})
        if proc.returncode: raise RuntimeError(f"Failed stage: {label}")
    worker = str(source/"broad_campaign_worker.py")
    try:
        run("preflight",["python",worker,"preflight",str(out)])
        for arm in [a["id"] for a in manifest["arms"]]:
            for stage in ("design","inverse_folding"):
                run(f"{arm}-{stage}",["python",worker,"stage",str(out),"--arm",arm,"--stage",stage])
        run("assemble",["python",worker,"assemble",str(out)])
        for seed in manifest["refolding"]["seeds"]:
            run(f"seed-{seed}",["python",str(source/"control_refold_worker.py"),
                                str(out/f"seed-{seed}"),"--seed",str(seed),"--export-pae","--samples",str(manifest["refolding"]["diffusion_samples_per_seed"])])
    except Exception:
        (out/"error.txt").write_text(traceback.format_exc())
        raise
    finally:
        result = {"run_id":run_id,"boltzgen_commit":PIN,"statuses":statuses,
                  "elapsed_seconds":time.monotonic()-started}
        (out/"run_status.json").write_text(json.dumps(result,indent=2))
        volume.commit()
    return result


@app.local_entrypoint()
def main(run_id: str):
    import re
    assert re.fullmatch(r"[A-Za-z0-9_-]+",run_id)
    ledger = ROOT/"outputs"/f"{run_id}_launch.json"
    reservation = json.loads(ledger.read_text())
    assert reservation["run_id"] == run_id and 0 < reservation["reserved_usd"] <= 6
    folder = ROOT/"data/apical_expansion"
    files = {p.name:p.read_bytes() for p in folder.iterdir() if p.is_file()}
    for name in ("broad_campaign_worker.py","control_refold_worker.py",
                 "prepare_controls.py","prepare_target.py","direction_check.py","generation_comparison_worker.py"):
        files[name] = (ROOT/"scripts"/name).read_bytes()
    import hashlib
    assert {k:hashlib.sha256(v).hexdigest() for k,v in files.items()} == reservation["input_sha256"]
    result = run_comparison.remote(files,run_id)
    (ROOT/"outputs"/f"{run_id}_result.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
