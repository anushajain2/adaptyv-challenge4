"""Fixed final-stage validation of at most six diverse initial expansion passes."""
import argparse
import hashlib
import json
from prepare_controls import ROOT
from prepare_next_round import finish
from analyze_broad_campaign import diverse_shortlist


def choose(result):
    if not result['controls_pass']:return []
    passing=[r for r in result['ranked_candidates'] if r['single_sample_full_screen']]
    return diverse_shortlist(passing,limit=6)


def main():
    p=argparse.ArgumentParser();p.add_argument('run');a=p.parse_args();run=(ROOT/a.run).resolve()
    report=ROOT/'reports/apical_expansion_fixed/results.json';r=json.loads(report.read_text())
    chosen=choose(r)
    if not chosen:raise SystemExit('No diverse full-screen passes with passing controls; no validation launch.')
    src=run/'references';out=ROOT/'data/expansion_validation';out.mkdir(exist_ok=False)
    m=json.loads((run/'inputs/manifest.json').read_text());ids=[x['id'] for x in chosen]
    cases=[c for c in json.loads((run/'cases.json').read_text()) if c['id'] in ids or c['arm']=='control']
    for c in cases:
        for ext in ('.cif','.npz','.pdb'):
            name=c['id']+ext;(out/name).write_bytes((src/name).read_bytes())
    m.update(cases=cases,selected_ids=ids,source_run=r['run'],source_results_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
             scope='Fresh-seed confirmation of at most six diverse full-screen expansion passes',
             refolding={**m['refolding'],'seeds':[1409,1511],'diffusion_samples_per_seed':5},
             max_prediction_seconds=1500,experiment_cap_usd=1.5,
             selection_details={'max_novel_sequences':6,'policy':'Require raw full-screen pass; round-robin length strata following frozen rank score, exclude prior exact sequences and >=80% identity to prior shortlisted leads; within-new-set 80% sequence identity diversity. No substitutions after results.', 'post_validation':'No automatic new generation or mutation search.'})
    finish(out,m)
    print(json.dumps({'selected_ids':ids,'predictions':len(cases)*10},indent=2))
if __name__=='__main__':main()
