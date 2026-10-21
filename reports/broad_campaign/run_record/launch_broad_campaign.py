"""Reserve a bounded campaign within the approved US$30 gross total."""
import argparse
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime,timezone
from pathlib import Path
from launch_generation_comparison import gross_cost
ROOT=Path(__file__).resolve().parents[1]
WORKERS=('broad_campaign_worker.py','control_refold_worker.py','prepare_controls.py',
         'prepare_target.py','direction_check.py','generation_comparison_worker.py')

def bundled_files():
    files={p.name:p.read_bytes() for p in (ROOT/'data/broad_campaign').iterdir() if p.is_file()}
    files.update({name:(ROOT/'scripts'/name).read_bytes() for name in WORKERS})
    return files

def reservation_bound():
    return (6750+600+2)*(.000542+4*.0000131+32*.00000222)+.75

def main():
    p=argparse.ArgumentParser();p.add_argument('--run-id',default='broad-009-6wrw');args=p.parse_args()
    assert re.fullmatch(r'[A-Za-z0-9_-]+',args.run_id)
    ledger=ROOT/'outputs'/f'{args.run_id}_launch.json'
    if list((ROOT/'outputs').glob('broad-*_launch.json')):
        raise SystemExit('Campaign already attempted; reconcile existing run before any retry or expansion.')
    modal=str(ROOT/'.venv/bin/modal')
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45))
    if any(a.get('state')!='stopped' for a in apps):raise SystemExit('Outstanding apps must stop before launch.')
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    prior=json.loads((ROOT/'reports/generation_comparison/run_record/billing_after.json').read_text())
    gross=max(gross_cost(billing),gross_cost(prior))
    assert reservation_bound()<6 and gross+6<=30 and gross+16+7+3<=30
    files=bundled_files();manifest=json.loads(files['manifest.json'])
    assert sum(a['num_designs'] for a in manifest['arms'])==200
    assert manifest['max_prediction_seconds']==6600 and manifest['experiment_cap_usd']==6
    assert manifest['refolding']['diffusion_samples_per_seed']==1
    for name,sha in manifest['input_sha256'].items():assert hashlib.sha256(files[name]).hexdigest()==sha
    report=ROOT/'reports/broad_campaign/run_record';report.mkdir(parents=True,exist_ok=True)
    for name,data in [('billing_before',billing),('apps_before',apps)]:
        (report/f'{name}.json').write_text(json.dumps(data,indent=2))
    record={'run_id':args.run_id,'started_utc':datetime.now(timezone.utc).isoformat(),
        'prior_gross_high_water_usd':gross,'reserved_usd':6,'authorized_total_usd':30,
        'campaign_allocation_usd':16,'deep_evaluation_reserve_usd':7,'contingency_usd':3,
        'conservative_resource_and_misc_bound_usd':reservation_bound(),
        'runtime_bound':{'work_seconds':6600,'function_timeout_seconds':6750,'startup_timeout_seconds':600,'max_containers':1,'retries':0},
        'pricing_source':'https://modal.com/pricing','pricing_checked_date':'2026-10-21',
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},'status':'launched_or_attempted'}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    (report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
    (report/'manifest.json').write_bytes(files['manifest.json'])
    for name in WORKERS+('launch_broad_campaign.py','modal_broad_campaign.py','prepare_broad_campaign.py'):
        (report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    result=subprocess.run([modal,'run','--detach',str(ROOT/'scripts/modal_broad_campaign.py'),'--run-id',args.run_id],
        cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=result.returncode
    ledger.write_text(json.dumps(record,indent=2));(report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
    raise SystemExit(result.returncode)
if __name__=='__main__':main()
