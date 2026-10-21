"""Verify a downloaded phase against its prospective launch and archive provenance."""
import argparse,hashlib,json,shutil
from pathlib import Path
from prepare_controls import ROOT
from launch_expansion_his import PHASES


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=PHASES);p.add_argument('run',type=Path);a=p.parse_args();rid,folder=PHASES[a.phase];run=a.run.resolve();dest=ROOT/'reports'/folder/'run_record'
    assert run.name==rid;ledger=json.loads((dest/'launch_ledger.json').read_text());manifest=json.loads((run/'inputs/manifest.json').read_text())
    assert json.loads((run/'input_hashes.json').read_text())==ledger['input_sha256']
    for name,sha in manifest['input_sha256'].items():assert hashlib.sha256((ROOT/'data'/folder/name).read_bytes()).hexdigest()==sha,name
    assert hashlib.sha256((ROOT/'scripts/analyze_his_panel.py').read_bytes()).hexdigest()==ledger['analysis_sha256']['analyze_his_panel.py']
    for name,sha in ledger['analysis_sha256'].items():assert hashlib.sha256((dest/name).read_bytes()).hexdigest()==sha,name
    status=json.loads((run/'run_status.json').read_text());assert len(status['statuses'])==len(manifest['refolding']['seeds']) and all(x['returncode']==0 for x in status['statuses'])
    for name in ['input_hashes.json','requirements-cloud.txt','run_status.json']:shutil.copyfile(run/name,dest/name)
    for seed in manifest['refolding']['seeds']:
        pre=json.loads((run/f'seed-{seed}/preflight.json').read_text());assert len(pre)==len(manifest['cases']) and all(x['binder_template_excluded'] and x['msa_depth']==1 for x in pre)
        shutil.copyfile(run/f'seed-{seed}/preflight.json',dest/f'seed-{seed}_preflight.json')
    record={'all_input_hashes_match':True,'analysis_matches_prelaunch_snapshot':True,'archived_source_hashes_match':True,'preflight_case_count_and_masks_pass':True,'worker_statuses_pass':True}
    (dest/'verification.json').write_text(json.dumps(record,indent=2));print(json.dumps(record))

if __name__=='__main__':main()
