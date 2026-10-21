"""Reserve at most $4 for this two-stage His panel within the $30 gross cap."""
import argparse,hashlib,json,os,subprocess
from datetime import datetime,timezone
from prepare_controls import ROOT
from launch_generation_comparison import gross_cost

PHASES={'panel':('his-expansion-020-6wrw','expansion_his_panel'),'confirmation':('his-expansion-confirm-021-6wrw','expansion_his_confirmation')}
RATE=.000542+4*.0000131+32*.00000222
BOUND=(2250+120+2)*RATE+.4
SOURCE_FILES=['prepare_expansion_his.py','launch_expansion_his.py','control_refold_worker.py','analyze_his_panel.py','prepare_selected_models.py','structural_preparation.py','check_expansion_his_direction.py','audit_expansion_his_direction.py','prepare_next_round.py','direction_check.py','target_integrity_selection.py','analyze_generation_comparison.py','prepare_controls.py','score_matched_controls.py','summarize_direction.py']


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=PHASES);a=p.parse_args();run,folder=PHASES[a.phase]
    ledger=ROOT/'outputs'/f'{run}_launch.json'
    if ledger.exists():raise SystemExit('Already attempted; no automatic relaunch or seed retries.')
    src=ROOT/'data'/folder;m=json.loads((src/'manifest.json').read_text());assert m['experiment_cap_usd']==2 and m['max_prediction_seconds']==2100 and BOUND<2
    for name,sha in m['input_sha256'].items():assert hashlib.sha256((src/name).read_bytes()).hexdigest()==sha,name
    modal=str(ROOT/'.venv/bin/modal');apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45));assert all(x['state']=='stopped' for x in apps),'Unaccounted active app'
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    # Retain the ENTIRE prior-round reservation as a conservative allowance for meter lag.
    origin=max(gross_cost(billing),9.9576284,7.96348782+12)
    previous=[]
    for rid,_ in PHASES.values():
        path=ROOT/'outputs'/f'{rid}_launch.json'
        if path.exists():previous.append(json.loads(path.read_text()))
    assert sum(x['reserved_usd'] for x in previous)+2<=4
    assert origin+sum(x['reserved_usd'] for x in previous)+2<=30
    if a.phase=='confirmation':
        assert len(previous)==1 and previous[0]['client_exit_code']==0
        assert json.loads((ROOT/'reports/expansion_his_panel/results.json').read_text())['controls_pass']
    files={p.name:p.read_bytes() for p in src.iterdir() if p.suffix in {'.cif','.npz','.json','.a3m'}};files['control_refold_worker.py']=(ROOT/'scripts/control_refold_worker.py').read_bytes()
    report=ROOT/'reports'/folder/'run_record';report.mkdir(exist_ok=False)
    for name,v in [('billing_before',billing),('apps_before',apps),('manifest',m),('direction_plan',json.loads((ROOT/'reports/expansion_his_panel/direction_plan.json').read_text()))]:(report/f'{name}.json').write_text(json.dumps(v,indent=2))
    for name in SOURCE_FILES+[f'modal_{folder}.py']:(report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    record={'run_id':run,'phase':a.phase,'input_dir':f'data/{folder}','started_utc':datetime.now(timezone.utc).isoformat(),'reserved_usd':2,'total_new_round_reservation_usd':4,'authorized_total_usd':30,'prior_conservative_gross_allowance_usd':origin,'resource_and_misc_bound_usd':BOUND,'pricing_source':'https://modal.com/pricing','pricing_checked_date':'2026-10-21','runtime_bound':{'work_seconds':2100,'function_timeout_seconds':2250,'startup_timeout_seconds':120,'max_containers':1,'retries':0},'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},'analysis_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in report.glob('*.py')}}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    (report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
    proc=subprocess.run([modal,'run','--detach',str(ROOT/'scripts'/f'modal_{folder}.py'),'--run-id',run,'--input-dir',f'data/{folder}'],cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=proc.returncode;ledger.write_text(json.dumps(record,indent=2));(report/'launch_ledger.json').write_text(json.dumps(record,indent=2));raise SystemExit(proc.returncode)

if __name__=='__main__':main()
