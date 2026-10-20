"""Launch only after authentication and a dashboard US$30 gross usage cap."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace-usage-cap-confirmed",type=float,required=True,
                   help="Actual configured gross Workspace usage cap, before credits; must be <=30.")
    p.add_argument("--run-id",default="pilot-001")
    p.add_argument("--hotspot",default="conserved_326_347")
    p.add_argument("--protocol",choices=["all","protein-anything","nanobody-anything"],default="all")
    p.add_argument("--reviewed-retry",action="store_true",
                   help="After diagnosing a finished attempt and reconciling its cost; never retries automatically.")
    p.add_argument("--billing-record",type=Path,
                   help="Saved current Modal billing summary for a reviewed retry.")
    args=p.parse_args()
    if not 0<args.workspace_usage_cap_confirmed<=30:
        raise SystemExit("The approved total usage cap is US$30 including credits.")
    # Do not print credential values.
    if not Path.home().joinpath(".modal.toml").exists() and not (
        os.environ.get("MODAL_TOKEN_ID") and os.environ.get("MODAL_TOKEN_SECRET")
    ) and not (os.environ.get("MODAL_CONFIG_PATH") and Path(os.environ["MODAL_CONFIG_PATH"]).exists()):
        raise SystemExit("Authenticate first: .venv/bin/modal setup")
    ledger=ROOT/"outputs/cloud_launch.json"
    ledger.parent.mkdir(exist_ok=True)
    history=[]
    if ledger.exists():
        if not args.reviewed_retry or not args.billing_record:
            raise SystemExit("Prior attempt exists: diagnose it and provide a billing record before a reviewed retry.")
        prior=json.loads(ledger.read_text())
        if prior["run_id"]==args.run_id:
            raise SystemExit("Use a fresh run id so prior outputs are preserved.")
        billing=json.loads(args.billing_record.read_text())
        usage=max(float(billing["metered_cost"]),sum(
            float(v) for v in billing.get("metered_cost_breakdown",{}).values()))
        if usage+20>args.workspace_usage_cap_confirmed:
            raise SystemExit("Insufficient remaining budget for the US$20 planning reservation.")
        history=prior.pop("prior_attempts",[])+[prior]
    elif args.reviewed_retry:
        raise SystemExit("No prior attempt to review.")
    with ledger.open("w" if history else "x") as f:
        json.dump({"authorized_total_usd":30,"reserved_usd":20,
                   "workspace_usage_cap_confirmed_usd":args.workspace_usage_cap_confirmed,
                   "run_id":args.run_id,"status":"launched_or_attempted",
                   "hotspot":args.hotspot,
                   "protocol":args.protocol,
                   "prior_attempts":history,
                   "billing_record":str(args.billing_record) if args.billing_record else None,
                   "started_utc":datetime.now(timezone.utc).isoformat(),
                   "note":"Includes build attempt. Reconcile actual Modal charges before any relaunch."},f,indent=2)
    env=dict(os.environ,TFR1_LAUNCH_RESERVED="1")
    result=subprocess.run([str(ROOT/".venv/bin/modal"),"run","--detach",
                           str(ROOT/"scripts/modal_pilot.py"),"--run-id",args.run_id,
                           "--hotspot",args.hotspot,"--protocol",args.protocol],cwd=ROOT,env=env)
    record=json.loads(ledger.read_text())
    record["client_exit_code"]=result.returncode
    ledger.write_text(json.dumps(record,indent=2))
    raise SystemExit(result.returncode)


if __name__=="__main__":
    main()
