"""Prepare every selected model once using unchanged v2 numerical settings."""
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from receptor_preparation_v2 import PROTOCOL, prepare
from novelty_search import ROOT, sha, save
from build_local_redesign_panel import COHORT
PROTOCOL['cohort']=COHORT


def main():
    source=ROOT/'reports/local_redesign_panel/results.json'
    out=ROOT/'reports/local_redesign_panel/preparation';out.mkdir(exist_ok=False)
    raw=json.loads(source.read_text());m=json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text())
    assert PROTOCOL==m['structural_preparation']
    selected=[r for r in raw['rows'] if r['confidence_selected']]
    plan=[{'id':r['id'],'seed':r['seed'],'sample':r['sample'],'source':r['raw_scoring_output']+'/model.pdb',
           'dest':str((out/'models'/f"{r['id']}_seed{r['seed']}").relative_to(ROOT))} for r in selected]
    result={'protocol':PROTOCOL,'source_results_sha256':sha(source),'source_hashes':{r['source']:sha(ROOT/r['source']) for r in plan},'plan':plan,'rows':[]}
    save(out/'plan.json',result)
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        for future in as_completed([pool.submit(prepare,row) for row in plan]):
            row=future.result();result['rows'].append(row);save(out/'results.json',result)
            print(len(result['rows']),'/',len(plan),row['id'],row['status'],'movement',row.get('displacement',{}).get('pass'),flush=True)
    assert len(result['rows'])==len(plan)


if __name__=='__main__':main()
