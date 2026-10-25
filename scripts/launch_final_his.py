"""Reserve final two-seed panel and conditional tag follow-up, at most $3 total."""
import argparse
import json
import os
import subprocess
from novelty_search import ROOT, sha, save, now
from launch_generation_comparison import gross_cost

PHASES={
    'panel':('final-his-025-6wrw','final_his_panel',1.5,1500,1650,.25),
    'tag':('final-his-tag-026-6wrw','final_his_tag',1.5,1500,1650,.25),
}

RATE=.000542+4*.0000131+32*.00000222
# Prior jobs stopped four days ago; reconcile unused reservations. Retain historical gross high-water plus $1 lag margin.
PRIOR_ALLOWANCE=12.23518765+1
WORKERS=('apical_expansion_worker.py','control_refold_worker.py','prepare_controls.py',
         'prepare_target.py','direction_check.py','generation_comparison_worker.py')


def resource_bound(phase):
    return (PHASES[phase][4]+120+2)*RATE+PHASES[phase][5]


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=PHASES);phase=p.parse_args().phase
    run,folder,reserve,work,timeout,misc=PHASES[phase]
    ledger=ROOT/'outputs'/f'{run}_launch.json'
    assert not ledger.exists(),'Already attempted; no automatic retry.'
    src=ROOT/'data'/folder;m=json.loads((src/'manifest.json').read_text())
    assert m['max_prediction_seconds']==work and m['experiment_cap_usd']==reserve
    assert resource_bound(phase)<reserve
    for name,digest in m['input_sha256'].items():assert sha(src/name)==digest
    modal=str(ROOT/'.venv/bin/modal')
    apps=json.loads(subprocess.check_output([modal,'app','list','--json'],text=True,timeout=45))
    assert all(a['state']=='stopped' for a in apps),'Active apps; reconcile first.'
    billing=json.loads(subprocess.check_output([modal,'billing','summary','--json'],text=True,timeout=45))
    origin=max(PRIOR_ALLOWANCE,gross_cost(billing))
    assert origin+3<=30,'Insufficient headroom for complete $3 round.'
    if phase=='tag':
        assert (ROOT/'reports/final_his_direction/results.json').exists()
        assert m['tag_selected_ids'] and len(m['tag_selected_ids'])<=2
    report=ROOT/'reports'/folder/'run_record';report.mkdir(parents=True,exist_ok=False)
    save(report/'billing_before.json',billing);save(report/'apps_before.json',apps)
    workers=('control_refold_worker.py',)
    files={p.name:p.read_bytes() for p in src.iterdir() if p.is_file() and (phase=='generation' or p.suffix in {'.cif','.npz','.json','.a3m'})}
    files.update({name:(ROOT/'scripts'/name).read_bytes() for name in workers})
    names=set(workers)|{'prepare_final_his.py','launch_final_his.py','modal_final_his.py',
        'analyze_apical_expansion.py','analyze_diversity_round.py','analyze_his_panel.py','analyze_generation_comparison.py',
        'target_integrity_selection.py','score_matched_controls.py','audit_control_samples.py','receptor_preparation_v2.py','prepare_final_his_models.py','check_final_his_direction.py'}
    for name in names:(report/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    import hashlib
    record={'run_id':run,'input_dir':str(src.relative_to(ROOT)),'started_utc':now(),'reserved_usd':reserve,
        'authorized_total_usd':30,'round_reservation_usd':3,'prior_conservative_allowance_usd':origin,
        'resource_and_misc_bound_usd':resource_bound(phase),'pricing_checked_date':'2026-10-25','pricing_source':'https://modal.com/pricing',
        'runtime_bound':{'work_seconds':work,'function_timeout_seconds':timeout,'startup_timeout_seconds':120,'max_containers':1,'retries':0},
        'input_sha256':{k:hashlib.sha256(v).hexdigest() for k,v in files.items()},
        'analysis_sha256':{p.name:sha(p) for p in report.glob('*.py')}}
    with ledger.open('x') as f:json.dump(record,f,indent=2)
    save(report/'launch_ledger.json',record)
    cmd=[modal,'run','--detach',str(ROOT/'scripts/modal_final_his.py'),'--run-id',run]
    cmd+=['--input-dir',str(src.relative_to(ROOT))]
    result=subprocess.run(cmd,cwd=ROOT,env=dict(os.environ,TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code']=result.returncode;save(ledger,record);save(report/'launch_ledger.json',record)
    raise SystemExit(result.returncode)


if __name__=='__main__':main()
