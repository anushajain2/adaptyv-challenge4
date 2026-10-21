"""Reserve the approved $3 diagnostic panel/confirmation and $6 apical expansion."""
import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from launch_generation_comparison import gross_cost
from prepare_controls import ROOT

PHASES={
 'panel':('his-panel-014-6wrw','his_panel',1.5,1500,1650,120),
 'expansion':('apical-015-6wrw','apical_expansion',6,6600,6750,600),
 'confirmation':('his-confirm-016-6wrw','his_confirmation',1.5,1500,1650,120),
}
WORKERS=('broad_campaign_worker.py','control_refold_worker.py','prepare_controls.py','prepare_target.py','direction_check.py','generation_comparison_worker.py')
RATE=.000542+4*.0000131+32*.00000222

def resource_bound(phase):
 _,_,_,_,timeout,startup=PHASES[phase]
 return (timeout+startup+2)*RATE+(.75 if phase=='expansion' else .25)

def main():
 p=argparse.ArgumentParser();p.add_argument('phase',choices=PHASES);args=p.parse_args()
 run,folder,reserve,work,timeout,startup=PHASES[args.phase]
 ledger=ROOT/'outputs'/f'{run}_launch.json'
 if ledger.exists():raise SystemExit('This phase was already attempted. No automatic retries.')
 src=ROOT/'data'/folder;m=json.loads((src/'manifest.json').read_text())
 for name,sha in m['input_sha256'].items():assert hashlib.sha256((src/name).read_bytes()).hexdigest()==sha
 assert m['max_prediction_seconds']==work and m['experiment_cap_usd']==reserve
 assert resource_bound(args.phase)<reserve
 prep=json.loads((ROOT/'reports/structural_preparation/results.json').read_text())
 assert len(prep['rows'])==12 and all(r['status']=='ok' for r in prep['rows']), 'Complete uniform preparation diagnostics before launch'
 modal=str(ROOT/'.venv/bin/modal')
 apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45))
 previous=[];allowed=[]
 for phase,(rid,_,allocation,*_) in PHASES.items():
  file=ROOT/'outputs'/f'{rid}_launch.json'
  if file.exists():
   record=json.loads(file.read_text());previous.append(record)
   if record.get('client_exit_code')==0:allowed.append(('tfr1-next-'+phase, None))
   elif record.get('app_id'):allowed.append(('tfr1-next-'+phase,record['app_id']))
 assert all(a['state']=='stopped' or any(a['description']==description and (app_id is None or a['app_id']==app_id) for description,app_id in allowed) for a in apps), 'Unaccounted active app'
 billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
 # Count full earlier reservations even if the meter has partially reconciled them.
 prior=max(gross_cost(billing),7.96348782)+sum(r['reserved_usd'] for r in previous)
 assert sum(r['reserved_usd'] for r in previous)+reserve<=9
 assert prior+reserve<=30
 files={p.name:p.read_bytes() for p in src.iterdir() if p.is_file() and (args.phase=='expansion' or p.suffix in {'.cif','.npz','.json','.a3m'})}
 workers=WORKERS if args.phase=='expansion' else ('control_refold_worker.py',)
 files.update({name:(ROOT/'scripts'/name).read_bytes() for name in workers})
 record={'run_id':run,'input_dir':f'data/{folder}','phase':args.phase,'started_utc':datetime.now(timezone.utc).isoformat(),
  'prior_gross_high_water_usd':prior,'reserved_usd':reserve,'authorized_total_usd':30,'next_round_cap_usd':9,
  'conservative_resource_and_misc_bound_usd':resource_bound(args.phase),
  'runtime_bound':{'work_seconds':work,'function_timeout_seconds':timeout,'startup_timeout_seconds':startup,'max_containers':1,'retries':0},
  'pricing_source':'https://modal.com/pricing','pricing_checked_date':'2026-10-21',
  'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},'status':'launched_or_attempted'}
 report=ROOT/'reports'/folder/'run_record';report.mkdir(parents=True,exist_ok=False)
 for name,value in [('billing_before',billing),('apps_before',apps),('manifest',m),('structural_preparation_protocol',prep['protocol'])]:
  (report/f'{name}.json').write_text(json.dumps(value,indent=2))
 for name in set(workers+('launch_next_round.py',f'modal_{folder}.py','prepare_next_round.py','analyze_his_panel.py','analyze_broad_campaign.py','analyze_apical_expansion.py','structural_preparation.py','target_integrity_selection.py')):
  (report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
 record['analysis_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in report.glob('*.py')}
 with ledger.open('x') as f:json.dump(record,f,indent=2)
 (report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
 cmd=[modal,'run','--detach',str(ROOT/'scripts'/f'modal_{folder}.py'),'--run-id',run]
 if args.phase!='expansion':cmd+=['--input-dir',f'data/{folder}']
 result=subprocess.run(cmd,cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
 record['client_exit_code']=result.returncode
 ledger.write_text(json.dumps(record,indent=2));(report/'launch_ledger.json').write_text(json.dumps(record,indent=2))
 raise SystemExit(result.returncode)
if __name__=='__main__':main()
