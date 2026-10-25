"""Audit the downloaded run against launch-time code, inputs and feature checks."""
import argparse
import json
import shutil
from pathlib import Path
from novelty_search import ROOT,sha,save,now
from launch_reserve_tag import RUN


def main():
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);run=parser.parse_args().run.resolve()
    assert run.name==RUN
    dest=ROOT/'reports/reserve_tag/run_record';ledger=json.loads((dest/'launch_ledger.json').read_text())
    assert json.loads((run/'input_hashes.json').read_text())==ledger['input_sha256']
    for name,digest in ledger['input_sha256'].items():assert sha(run/'inputs'/name)==digest,name
    for name,digest in ledger['analysis_sha256'].items():assert sha(dest/name)==digest and sha(ROOT/'scripts'/name)==digest,name
    m=json.loads((run/'inputs/manifest.json').read_text())
    assert sha(run/'inputs/manifest.json')==sha(ROOT/'data/reserve_tag/manifest.json')
    for name,digest in m['input_sha256'].items():assert sha(ROOT/'data/reserve_tag'/name)==digest
    status=json.loads((run/'run_status.json').read_text())
    assert status['statuses']==[{'seed':s,'returncode':0} for s in m['refolding']['seeds']]
    assert status['boltzgen_commit']=='247b9bbd8b68a60aba854c2968d6a0cddd21ad6d'
    assert not (run/'error.txt').exists() and status['elapsed_seconds']<=2100
    for name in ['input_hashes.json','requirements-cloud.txt','run_status.json']:
        shutil.copyfile(run/name,dest/name)
    for seed in m['refolding']['seeds']:
        pre=json.loads((run/f'seed-{seed}/preflight.json').read_text())
        assert {p['id'] for p in pre}=={c['id'] for c in m['cases']}
        assert all(p['binder_template_excluded'] and p['design_mask_verified'] and p['msa_depth']==1 for p in pre)
        shutil.copyfile(run/f'seed-{seed}/preflight.json',dest/f'seed-{seed}_preflight.json')
    save(dest/'verification.json',{'status':'PASS','checked_utc':now(),'all_inputs_and_sources_match_launch':True,
        'all_cloud_stages_successful':True,'target_only_templates_and_empty_MSA':True,'input_manifest_sha256':sha(run/'inputs/manifest.json')})
    print('Verified cloud inputs, pinned implementation, completion and template/MSA masks.')


if __name__=='__main__':main()
