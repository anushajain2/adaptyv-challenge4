"""Prepare fresh-seed checks only for the completed campaign's eligible shortlist."""
import argparse
import hashlib
import json
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def prepare(run,out):
    result=json.loads((ROOT/'reports/broad_campaign/results.json').read_text())
    assert result['run']==run.name and result['controls_pass']
    chosen=result['shortlisted_ids'];assert 0<len(chosen)<=30
    ranked={r['id']:r for r in result['ranked_candidates']}
    assert all(ranked[n]['geometry_eligible'] for n in chosen)
    source_manifest=json.loads((run/'inputs/manifest.json').read_text())
    cases=json.loads((run/'cases.json').read_text())
    cases=[c for c in cases if c['id'] in chosen or c['arm']=='control']
    assert len(cases)==len(chosen)+3
    out.mkdir(parents=True,exist_ok=False)
    for c in cases:
        for ext in ('.cif','.npz','.pdb'):(out/(c['id']+ext)).write_bytes((run/'references'/(c['id']+ext)).read_bytes())
    m={**source_manifest,'scope':'Deep evaluation of previously selected one-sample leads; not a new generation batch',
       'prepared_utc':datetime.now(timezone.utc).isoformat(),'source_run':run.name,'cases':cases,
       'selected_ids':chosen,'selection_result_sha256':hashlib.sha256((ROOT/'reports/broad_campaign/results.json').read_bytes()).hexdigest(),
       'refolding':{**source_manifest['refolding'],'seeds':[401,503],'diffusion_samples_per_seed':5},
       'max_prediction_seconds':3600,'experiment_cap_usd':4,
       'screen':{**source_manifest['screen'],'replicates':'Both fresh seeds; confidence-selected sample from five samples per seed',
                 'decision':'Miniproteins passing both seeds are provisional computational parents. VHHs retain the missing-positive-control limitation. No affinity, mouse binding, or pH claim.'}}
    m['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(m,indent=2))
    print(json.dumps({'novel_candidates':len(chosen),'predictions':len(cases)*10,'output':str(out)},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--out',type=Path,default=ROOT/'data/broad_validation');a=p.parse_args();prepare(a.run,a.out)
