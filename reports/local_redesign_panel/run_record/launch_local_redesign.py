"""Reserve the approved $6 round before launching the bounded $2 panel."""
import json
import os
import subprocess
import hashlib
from novelty_search import ROOT, sha, save, now
from launch_generation_comparison import gross_cost

RUN='local-redesign-027-6wrw'
RATE=.000542+4*.0000131+32*.00000222


def main():
    ledger=ROOT/'outputs'/f'{RUN}_launch.json';assert not ledger.exists(),'No automatic relaunch'
    src=ROOT/'data/local_redesign_panel';m=json.loads((src/'manifest.json').read_text())
    assert m['max_prediction_seconds']==1500 and m['experiment_cap_usd']==2
    for name,digest in m['input_sha256'].items():assert sha(src/name)==digest
    bound=(1650+120+2)*RATE+.25;assert bound<2
    modal=str(ROOT/'.venv/bin/modal')
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45));assert all(a['state']=='stopped' for a in apps)
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    previous=json.loads((ROOT/'reports/final_his_panel/run_record/budget_reconciliation.json').read_text())
    prior=max(previous['conservative_current_total_usd'],gross_cost(billing))+1
    assert prior+6<=30
    report=ROOT/'reports/local_redesign_panel/run_record';report.mkdir(exist_ok=False)
    save(report/'billing_before.json',billing);save(report/'apps_before.json',apps)
    files={p.name:p.read_bytes() for p in src.iterdir() if p.suffix in {'.cif','.npz','.json','.a3m'}}
    files['control_refold_worker.py']=(ROOT/'scripts/control_refold_worker.py').read_bytes()
    names=['control_refold_worker.py','modal_final_his.py','launch_local_redesign.py',
           'prepare_local_redesign.py','run_local_mpnn.py','build_local_redesign_panel.py',
           'audit_preparation_motion.py','analyze_his_panel.py','target_integrity_selection.py',
           'prepare_local_redesign_models.py','analyze_local_redesign_direction.py','receptor_preparation_v2.py',
           'analyze_generation_comparison.py','score_matched_controls.py','audit_control_samples.py']
    for name in names:(report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    record={'run_id':RUN,'input_dir':str(src.relative_to(ROOT)),'started_utc':now(),'reserved_usd':2,
        'authorized_total_usd':30,'round_reservation_usd':6,'prior_conservative_allowance_usd':prior,
        'resource_and_misc_bound_usd':bound,'pricing_checked_date':'2026-10-25','pricing_source':'https://modal.com/pricing',
        'runtime_bound':{'work_seconds':1500,'function_timeout_seconds':1650,'startup_timeout_seconds':120,'max_containers':1,'retries':0},
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
        'analysis_sha256':{p.name:sha(p) for p in report.glob('*.py')}}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    save(report/'launch_ledger.json',record)
    result=subprocess.run([modal,'run','--detach',str(ROOT/'scripts/modal_final_his.py'),'--run-id',RUN,'--input-dir',str(src.relative_to(ROOT))],cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=result.returncode;save(ledger,record);save(report/'launch_ledger.json',record)
    raise SystemExit(result.returncode)


if __name__=='__main__':main()
