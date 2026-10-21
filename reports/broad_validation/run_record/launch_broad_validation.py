"""Reserve fresh-seed shortlist evaluation inside the approved deep-evaluation allocation."""
import hashlib
import json
import os
import subprocess
from datetime import datetime,timezone
from pathlib import Path
from launch_generation_comparison import gross_cost
ROOT=Path(__file__).resolve().parents[1]

def main():
    run_id='validation-010-6wrw';folder=ROOT/'data/broad_validation';ledger=ROOT/'outputs/broad_validation_launch.json'
    if ledger.exists():raise SystemExit('Validation already attempted. No automatic retry.')
    manifest=json.loads((folder/'manifest.json').read_text())
    for name,sha in manifest['input_sha256'].items():
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==sha
    assert 0<len(manifest['selected_ids'])<=30 and manifest['refolding']['seeds']==[401,503]
    broad=json.loads((ROOT/'reports/broad_campaign/results.json').read_text())
    assert broad['controls_pass'] and broad['shortlisted_ids']==manifest['selected_ids']
    modal=str(ROOT/'.venv/bin/modal')
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45))
    if any(a.get('state')!='stopped' for a in apps):raise SystemExit('Outstanding jobs must stop first.')
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    after=json.loads((ROOT/'reports/broad_campaign/run_record/billing_after.json').read_text())
    initial=json.loads((ROOT/'outputs/broad-009-6wrw_launch.json').read_text())
    # Until billing settles, retain the full initial reservation rather than spend lagging headroom twice.
    gross=max(gross_cost(billing),gross_cost(after),initial['prior_gross_high_water_usd']+initial['reserved_usd'])
    assert gross+4<=30 and 4<=7
    files={p.name:p.read_bytes() for p in folder.iterdir() if p.suffix in {'.cif','.npz','.json','.a3m'}}
    worker=ROOT/'scripts/control_refold_worker.py';files[worker.name]=worker.read_bytes()
    report=ROOT/'reports/broad_validation/run_record';report.mkdir(parents=True,exist_ok=True)
    record={'run_id':run_id,'input_dir':'data/broad_validation','started_utc':datetime.now(timezone.utc).isoformat(),
        'prior_gross_high_water_usd':gross,'reserved_usd':4,'authorized_total_usd':30,'deep_evaluation_allocation_usd':7,
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
        'runtime_bound':{'work_seconds':3600,'function_timeout_seconds':3750,'startup_timeout_seconds':600,'max_containers':1,'retries':0},
        'conservative_resource_bound_usd':(3750+600+2)*(.000542+4*.0000131+32*.00000222),
        'status':'launched_or_attempted'}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    for name,value in [('billing_before',billing),('apps_before',apps),('launch_ledger',record)]:
        (report/f'{name}.json').write_text(json.dumps(value,indent=2))
    (report/'manifest.json').write_bytes(files['manifest.json'])
    for name in ('control_refold_worker.py','modal_broad_validation.py','launch_broad_validation.py','prepare_broad_validation.py'):
        (report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    result=subprocess.run([modal,'run','--detach',str(ROOT/'scripts/modal_broad_validation.py'),'--run-id',run_id,'--input-dir','data/broad_validation'],
        cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=result.returncode
    ledger.write_text(json.dumps(record,indent=2));(report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
    raise SystemExit(result.returncode)
if __name__=='__main__':main()
