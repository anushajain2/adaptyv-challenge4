"""Reserve one bounded $2 run inside the authorized $30 gross cap."""
import hashlib,json,os,subprocess
from datetime import datetime,timezone
from prepare_controls import ROOT
from launch_generation_comparison import gross_cost

RUN='tag-comparison-022-6wrw'
RATE=.000542+4*.0000131+32*.00000222
BOUND=(2250+120+2)*RATE+.4

def main():
    ledger=ROOT/'outputs'/f'{RUN}_launch.json';assert not ledger.exists(),'No automatic retry'
    source=ROOT/'data/tag_comparison';m=json.loads((source/'manifest.json').read_text())
    assert m['experiment_cap_usd']==2 and BOUND<2
    for name,digest in m['input_sha256'].items():assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest
    modal=str(ROOT/'.venv/bin/modal')
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45));assert all(a['state']=='stopped' for a in apps)
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    # Keep prior reservations as meter-lag allowance, even though prior apps stopped.
    origin=max(gross_cost(billing),7.96348782+12+4);assert origin+2<=30
    report=ROOT/'reports/tag_comparison/run_record';report.mkdir(exist_ok=False)
    for name,data in [('billing_before',billing),('apps_before',apps),('manifest',m)]:
        (report/f'{name}.json').write_text(json.dumps(data,indent=2))
    names=['prepare_tag_comparison.py','launch_tag_comparison.py','modal_tag_comparison.py','analyze_tag_comparison.py','control_refold_worker.py','receptor_preparation_v2.py','audit_control_samples.py','analyze_generation_comparison.py','prepare_controls.py','score_matched_controls.py','target_integrity_selection.py','direction_check.py']
    for name in names:(report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    files={p.name:p.read_bytes() for p in source.iterdir() if p.suffix in {'.cif','.npz','.json','.a3m'}}
    files['control_refold_worker.py']=(ROOT/'scripts/control_refold_worker.py').read_bytes()
    record={'run_id':RUN,'input_dir':'data/tag_comparison','started_utc':datetime.now(timezone.utc).isoformat(),
        'reserved_usd':2,'authorized_total_usd':30,'prior_conservative_gross_allowance_usd':origin,'resource_and_misc_bound_usd':BOUND,
        'pricing_source':'https://modal.com/pricing','pricing_checked_date':'2026-10-21',
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
        'analysis_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in report.glob('*.py')}}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    (report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
    p=subprocess.run([modal,'run','--detach',str(ROOT/'scripts/modal_tag_comparison.py'),'--run-id',RUN,'--input-dir','data/tag_comparison'],cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=p.returncode;ledger.write_text(json.dumps(record,indent=2));(report/'launch_ledger.json').write_text(json.dumps(record,indent=2));raise SystemExit(p.returncode)

if __name__=='__main__':main()
