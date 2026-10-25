"""Reserve one $2 confirmation/tag comparison within the existing $30 gross cap."""
import hashlib
import json
import os
import subprocess
from novelty_search import ROOT,sha,save,now
from launch_generation_comparison import gross_cost

RUN='reserve-tag-028-6wrw'
RATE=.000542+4*.0000131+32*.00000222
BOUND=(2250+120+2)*RATE+.4


def main():
    ledger=ROOT/'outputs'/f'{RUN}_launch.json';assert not ledger.exists(),'No automatic retry'
    source=ROOT/'data/reserve_tag';m=json.loads((source/'manifest.json').read_text())
    assert m['experiment_cap_usd']==2 and m['max_prediction_seconds']==2100 and BOUND<2
    for name,digest in m['input_sha256'].items():assert sha(source/name)==digest
    modal=str(ROOT/'.venv/bin/modal')
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45));assert all(a['state']=='stopped' for a in apps)
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    previous=json.loads((ROOT/'reports/local_redesign_panel/run_record/budget_reconciliation.json').read_text())
    prior=max(gross_cost(billing),previous['conservative_current_total_usd'])+1;assert prior+2<=30
    report=ROOT/'reports/reserve_tag/run_record';report.mkdir(exist_ok=False)
    save(report/'billing_before.json',billing);save(report/'apps_before.json',apps)
    names=['prepare_reserve_tag.py','launch_reserve_tag.py','analyze_reserve_tag.py','prepare_reserve_tag_models.py','analyze_reserve_tag_direction.py',
        'modal_tag_comparison.py','control_refold_worker.py','prepare_tag_models.py','prepare_tag_comparison.py','prepare_assay_tag_review.py',
        'analyze_tag_comparison.py','analyze_controls.py','analyze_generation_comparison.py','target_integrity_selection.py','audit_control_samples.py',
        'direction_check.py','score_matched_controls.py','audit_preparation_motion.py','audit_local_redesign_mechanism.py','analyze_local_redesign_direction.py','analyze_receptor_preparation_v2.py']
    for name in names:(report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    files={p.name:p.read_bytes() for p in source.iterdir() if p.is_file() and p.suffix in {'.cif','.npz','.json','.a3m'}}
    files['control_refold_worker.py']=(ROOT/'scripts/control_refold_worker.py').read_bytes()
    record={'run_id':RUN,'input_dir':'data/reserve_tag','started_utc':now(),'reserved_usd':2,'authorized_total_usd':30,
        'prior_conservative_allowance_usd':prior,'resource_and_misc_bound_usd':BOUND,'pricing_source':'https://modal.com/pricing','pricing_checked_date':'2026-10-26',
        'runtime_bound':{'work_seconds':2100,'function_timeout_seconds':2250,'startup_timeout_seconds':120,'max_containers':1,'retries':0},
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},'analysis_sha256':{p.name:sha(p) for p in report.glob('*.py')}}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    save(report/'launch_ledger.json',record)
    result=subprocess.run([modal,'run','--detach',str(ROOT/'scripts/modal_tag_comparison.py'),'--run-id',RUN,'--input-dir','data/reserve_tag'],cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=result.returncode;save(ledger,record);save(report/'launch_ledger.json',record)
    raise SystemExit(result.returncode)


if __name__=='__main__':main()
