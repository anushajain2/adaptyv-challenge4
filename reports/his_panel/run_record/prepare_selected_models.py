"""Apply one frozen preparation to every selected model, without score-based exceptions."""
import argparse
import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from prepare_controls import ROOT
from structural_preparation import prepare,PROTOCOL


def task(item):
    source,dest=item
    return prepare(ROOT/source,ROOT/dest)


def main():
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    data=json.loads(a.results.read_text())
    rows=([r for r in data['rows'] if r['confidence_selected']] if 'rows' in data else data['ranked_candidates']+data['controls'])
    plan=[{'id':r['id'],'seed':r['seed'],'sample':r['sample'],'source':r['raw_scoring_output']+'/model.pdb',
           'dest':str((a.output/'models'/f"{r['id']}_seed{r['seed']}").resolve().relative_to(ROOT))} for r in rows]
    a.output.mkdir(parents=True,exist_ok=False)
    record={'source_results_sha256':hashlib.sha256(a.results.read_bytes()).hexdigest(),'protocol':PROTOCOL,'plan':plan,'rows':[]}
    (a.output/'plan.json').write_text(json.dumps(record,indent=2))
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as executor:
        futures={executor.submit(task,(r['source'],r['dest'])):r for r in plan}
        for future in as_completed(futures):
            r=futures[future]
            record['rows'].append({**r,**future.result()})
            (a.output/'results.json').write_text(json.dumps(record,indent=2))
            print('Prepared',len(record['rows']),'/',len(plan),r['id'],flush=True)
    assert len(record['rows'])==len(plan)
if __name__=='__main__':main()
