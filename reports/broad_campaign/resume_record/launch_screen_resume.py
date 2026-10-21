"""Reconcile cancellation and reserve the remainder of the original $6 batch."""
import hashlib,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
from launch_generation_comparison import gross_cost
ROOT=Path(__file__).resolve().parents[1]
RATE=.000542+4*.0000131+32*.00000222

def main():
    ledger=ROOT/'outputs/screen_resume_launch.json'
    if ledger.exists():raise SystemExit('Continuation already attempted; no automatic retry.')
    modal=str(ROOT/'.venv/bin/modal');run_id='screen-resume-010-6wrw'
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45))
    if any(a.get('state')!='stopped' for a in apps):raise SystemExit('Active applications must be reconciled first.')
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    first=json.loads((ROOT/'outputs/broad-009-6wrw_launch.json').read_text())
    status=json.loads((ROOT/'outputs/broad009_complete/broad-009-6wrw/run_status.json').read_text())
    assert status['statuses'][-1]['stage']=='assemble' and all(s['returncode']==0 for s in status['statuses'])
    first_allowance=(status['elapsed_seconds']+600+2)*RATE+.25
    gross=max(gross_cost(billing),first['prior_gross_high_water_usd']+first_allowance)
    reservation=1.75;bound=(1900+120+2)*RATE+.25;ceiling=first['prior_gross_high_water_usd']+6
    assert bound<reservation and gross+reservation<=ceiling<=30
    folder=ROOT/'data/screen_resume'
    files={p.name:p.read_bytes() for p in folder.iterdir() if p.suffix in {'.cif','.npz','.json','.a3m'}}
    files['control_refold_worker.py']=(ROOT/'scripts/control_refold_worker.py').read_bytes()
    manifest=json.loads(files['manifest.json']);assert len(manifest['cases'])==131 and manifest['refolding']['seeds']==[307]
    for name,sha in manifest['input_sha256'].items():assert hashlib.sha256(files[name]).hexdigest()==sha
    record={'run_id':run_id,'input_dir':'data/screen_resume','started_utc':datetime.now(timezone.utc).isoformat(),
        'prior_gross_high_water_usd':gross,'reserved_usd':reservation,'authorized_total_usd':30,
        'shared_batch_ceiling_gross_usd':ceiling,'conservative_first_attempt_allowance_usd':first_allowance,
        'conservative_resume_bound_usd':bound,'reason':'Complete missing predictions after generic Modal cancellation. No model/preflight failure was reported. All 200 designs and 72 valid predictions preserved.',
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
        'runtime_bound':{'work_seconds':1800,'function_timeout_seconds':1900,'startup_timeout_seconds':120,'max_containers':1,'retries':0}}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    report=ROOT/'reports/broad_campaign/resume_record';report.mkdir(parents=True,exist_ok=True)
    for name,data in [('launch_ledger',record),('billing_before',billing),('apps_before',apps)]:
        (report/f'{name}.json').write_text(json.dumps(data,indent=2))
    (report/'manifest.json').write_bytes(files['manifest.json'])
    for n in ('prepare_screen_resume.py','launch_screen_resume.py','modal_screen_resume.py','control_refold_worker.py'):
        (report/n).write_bytes((ROOT/'scripts'/n).read_bytes())
    result=subprocess.run([modal,'run','--detach',str(ROOT/'scripts/modal_screen_resume.py'),'--run-id',run_id,'--input-dir','data/screen_resume'],
        cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=result.returncode;ledger.write_text(json.dumps(record,indent=2));(report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
    raise SystemExit(result.returncode)
if __name__=='__main__':main()
