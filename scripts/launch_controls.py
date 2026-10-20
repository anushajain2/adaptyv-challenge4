"""Reconcile live gross billing and reserve US$8 before the calibration run."""
import json
import argparse
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id", default="controls-001")
    p.add_argument("--reviewed-retry-reason")
    p.add_argument("--input-dir", default="data/controls")
    args = p.parse_args()
    import re
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.run_id):
        raise SystemExit("Invalid run id")
    input_dir = (ROOT / args.input_dir).resolve()
    if not input_dir.is_relative_to(ROOT) or not (input_dir / "manifest.json").exists():
        raise SystemExit("A prepared input directory inside the project is required")
    ledger = ROOT / "outputs/controls_launch.json"
    history = []
    if ledger.exists():
        prior = json.loads(ledger.read_text())
        if not args.reviewed_retry_reason or args.run_id == prior["run_id"]:
            raise SystemExit("Control launch already attempted; use a fresh id and diagnosed retry reason.")
        history = prior.pop("prior_attempts", []) + [prior]
    billing = json.loads(subprocess.check_output(
        [str(ROOT / ".venv/bin/modal"), "billing", "summary", "--json"], text=True, timeout=90))
    live = max(float(billing["metered_cost"]),
               sum(float(x) for x in billing["metered_cost_breakdown"].values()))
    previous = json.loads((ROOT / "reports/billing_final.json").read_text())
    gross = max(live, float(previous["metered_cost"]),
                sum(float(x) for x in previous["metered_cost_breakdown"].values()))
    if gross+8 > 30:
        raise SystemExit("Insufficient headroom under the authorized US$30 gross budget.")
    (ROOT / "reports/controls/billing_before.json").write_text(json.dumps(billing, indent=2))
    record = {"run_id": args.run_id, "authorized_total_usd": 30,
              "prior_gross_high_water_usd": gross, "reserved_usd": 8,
              "started_utc": datetime.now(timezone.utc).isoformat(),
              "status": "launched_or_attempted", "prior_attempts": history,
              "input_dir": str(input_dir.relative_to(ROOT)),
              "reviewed_retry_reason": args.reviewed_retry_reason,
              "note": "Existing user-confirmed US$30 gross workspace cap remains required."}
    with ledger.open("w" if history else "x") as f:
        json.dump(record, f, indent=2)
    result = subprocess.run(
        [str(ROOT / ".venv/bin/modal"), "run", "--detach",
         str(ROOT / "scripts/modal_controls.py"), "--run-id", args.run_id,
         "--input-dir", str(input_dir.relative_to(ROOT))],
        env=dict(os.environ, TFR1_LAUNCH_RESERVED="1"), cwd=ROOT)
    record["client_exit_code"] = result.returncode
    ledger.write_text(json.dumps(record, indent=2))
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
