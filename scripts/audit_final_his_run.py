"""Verify a downloaded bounded run against the launch-time evidence."""
import argparse
import json
import shutil
from pathlib import Path
from novelty_search import ROOT, sha, save, now
from launch_final_his import PHASES


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=PHASES);p.add_argument('run',type=Path)
    a=p.parse_args();rid,folder,*_=PHASES[a.phase];run=a.run.resolve()
    assert run.name==rid
    dest=ROOT/'reports'/folder/'run_record';ledger=json.loads((dest/'launch_ledger.json').read_text())
    assert json.loads((run/'input_hashes.json').read_text())==ledger['input_sha256']
    for name,digest in ledger['input_sha256'].items():assert sha(run/'inputs'/name)==digest,name
    for name,digest in ledger['analysis_sha256'].items():
        assert sha(dest/name)==digest and sha(ROOT/'scripts'/name)==digest,name
    m=json.loads((run/'inputs/manifest.json').read_text())
    for name,digest in m['input_sha256'].items():assert sha(ROOT/'data'/folder/name)==digest,name
    status=json.loads((run/'run_status.json').read_text())
    assert all(s['returncode']==0 for s in status['statuses'])
    assert len(status['statuses'])==(7 if a.phase=='generation' else 2)
    assert status['boltzgen_commit']==m['boltzgen_commit']
    assert not (run/'error.txt').exists()
    for name in ('input_hashes.json','requirements-cloud.txt','run_status.json'):
        shutil.copyfile(run/name,dest/name)
    if a.phase=='generation':
        for name in ('generation_preflight.json','actual_feature_audit.json','cases.json'):
            shutil.copyfile(run/name,dest/name)
        expected=sum(x['num_designs'] for x in m['arms'])+len(m['controls'])
    else:expected=len(m['cases'])
    for seed in m['refolding']['seeds']:
        pre=json.loads((run/f'seed-{seed}/preflight.json').read_text())
        assert len(pre)==expected and all(x['binder_template_excluded'] and x['msa_depth']==1 for x in pre)
        shutil.copyfile(run/f'seed-{seed}/preflight.json',dest/f'seed-{seed}_preflight.json')
    save(dest/'verification.json',{'checked_utc':now(),'all_inputs_and_sources_match_launch':True,
        'all_cloud_stages_successful':True,'target_only_templates_and_empty_MSA':True,
        'files':{str(p.relative_to(ROOT)):sha(p) for p in dest.iterdir() if p.is_file() and p.name!='verification.json'}})
    print('Verified',a.phase,'cloud inputs, source snapshots, stages and prediction preflight.')


if __name__=='__main__':main()
