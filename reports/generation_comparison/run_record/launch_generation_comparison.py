"""Reserve the authorized US$5 experiment inside the US$30 gross total."""
import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def gross_cost(record):
    return max(float(record["metered_cost"]),
               sum(float(v) for v in record.get("metered_cost_breakdown",{}).values()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id",default="generation-007-6wrw")
    p.add_argument("--reviewed-retry-of")
    p.add_argument("--retry-reason")
    args = p.parse_args()
    import re
    if not re.fullmatch(r"[A-Za-z0-9_-]+",args.run_id): raise SystemExit("Invalid run id")
    ledger = ROOT/"outputs"/f"{args.run_id}_launch.json"
    if ledger.exists(): raise SystemExit("Already attempted. No automatic relaunch or retry.")
    # Any repeat attempt must first be explicitly reconciled, not get a new $5 allowance.
    previous_attempts = list((ROOT/"outputs").glob("generation-*_launch.json"))
    prior = None
    if previous_attempts and not args.reviewed_retry_of:
        raise SystemExit("Prior generation experiment launch exists; reconcile it before another attempt.")
    if args.reviewed_retry_of:
        if not re.fullmatch(r"[A-Za-z0-9_-]+",args.reviewed_retry_of) or not args.retry_reason:
            raise SystemExit("A diagnosed retry requires a safe prior id and explanation")
        prior = json.loads((ROOT/"outputs"/f"{args.reviewed_retry_of}_launch.json").read_text())
        if prior.get("client_exit_code") != 1 or len(previous_attempts) != 1:
            raise SystemExit("Only the single diagnosed stopped preflight attempt may be retried here")
        status = json.loads((ROOT/"outputs/generation007_failed"/args.reviewed_retry_of/"run_status.json").read_text())
        if status["statuses"] != [s for s in status["statuses"] if s["stage"] == "preflight" and s["returncode"] == 1]:
            raise SystemExit("Retry adapter applies only to a preflight failure before generation")
    modal = str(ROOT/".venv/bin/modal")
    apps = json.loads(subprocess.check_output([modal,"app","list","--json"],text=True,timeout=45))
    if any(a.get("state") != "stopped" for a in apps):
        raise SystemExit("Active Modal applications found; reconcile outstanding spend before launch.")
    billing = json.loads(subprocess.check_output([modal,"billing","summary","--json"],text=True,timeout=45))
    old = json.loads((ROOT/"reports/controls_chai/billing_after.json").read_text())
    gross = max(gross_cost(billing),gross_cost(old))
    origin = prior["prior_gross_high_water_usd"] if prior else gross
    if prior:
        # Billing lags: reserve a conservative failed-attempt allowance including
        # maximum startup plus observed work time and $0.25 miscellaneous margin.
        failed_allowance = (600+status["elapsed_seconds"]+2)*(.000542+4*.0000131+32*.00000222)+.25
        gross = max(gross,origin+failed_allowance)
    remaining = origin+5-gross
    if origin+5 > 30 or remaining < 3.15:
        raise SystemExit("Insufficient experiment headroom for the bounded retry")
    folder = ROOT/"data/generation_comparison"
    files = {p.name:p.read_bytes() for p in folder.iterdir() if p.is_file()}
    for name in ("generation_comparison_worker.py","control_refold_worker.py",
                 "prepare_controls.py","prepare_target.py","direction_check.py"):
        files[name] = (ROOT/"scripts"/name).read_bytes()
    manifest = json.loads(files["manifest.json"])
    assert manifest["experiment_cap_usd"] == 5 and manifest["total_cap_usd"] == 30
    report = ROOT/"reports/generation_comparison/run_record"
    report.mkdir(parents=True,exist_ok=True)
    (report/"billing_before.json").write_text(json.dumps(billing,indent=2))
    (report/"apps_before.json").write_text(json.dumps(apps,indent=2))
    record = {"run_id":args.run_id,"started_utc":datetime.now(timezone.utc).isoformat(),
              "prior_gross_high_water_usd":gross,"reserved_usd":remaining,"authorized_total_usd":30,
              "experiment_origin_gross_usd":origin,"experiment_ceiling_gross_usd":origin+5,
              "reviewed_retry_of":args.reviewed_retry_of,"retry_reason":args.retry_reason,
              "failed_attempt_allowance_usd":failed_allowance if prior else 0,
              "input_sha256":{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
              "runtime_bound": {"work_seconds":3600,"function_timeout_seconds":3750,
                                "startup_timeout_seconds":600,"max_containers":1,"retries":0},
              "pricing": {"source":"https://modal.com/pricing","checked_date":"2026-10-21",
                          "L40S_per_second":0.000542,"max_cpu_cores":4,
                          "cpu_per_core_second":0.0000131,"max_memory_GiB":32,
                          "memory_per_GiB_second":0.00000222,
                          "conservative_runtime_resource_bound_usd":2.90,
                          "remaining_reservation_for_build_storage_egress_usd":2.10},
              "status":"launched_or_attempted"}
    with ledger.open("x") as f: json.dump(record,f,indent=2)
    (report/"launch_ledger.json").write_text(json.dumps(record,indent=2))
    (report/"manifest.json").write_bytes(files["manifest.json"])
    result = subprocess.run([modal,"run","--detach",str(ROOT/"scripts/modal_generation_comparison.py"),
                             "--run-id",args.run_id],cwd=ROOT,
                            env=dict(os.environ,TFR1_LAUNCH_RESERVED="1"))
    record["client_exit_code"] = result.returncode
    ledger.write_text(json.dumps(record,indent=2))
    (report/"launch_ledger.json").write_text(json.dumps(record,indent=2))
    raise SystemExit(result.returncode)


if __name__ == "__main__": main()
