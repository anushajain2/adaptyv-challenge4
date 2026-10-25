"""Apply unchanged v2 physics to every selected confirmation model, failures included."""
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from receptor_preparation_v2 import PROTOCOL, prepare
from novelty_search import ROOT, sha, save

# Update cohort metadata in this process and spawned workers; numerical settings unchanged.
PROTOCOL['cohort']='Every confidence-selected final His panel candidate and control, including failures.'


def main():
    source=ROOT/'reports/final_his_panel/results.json'
    out=ROOT/'reports/final_his_panel/preparation';out.mkdir(exist_ok=False)
    data=json.loads(source.read_text())
    selected=[r for r in data['rows'] if r['confidence_selected']]
    plan=[{'id':r['id'],'seed':r['seed'],'sample':r['sample'],
           'source':r['raw_scoring_output']+'/model.pdb',
           'dest':str((out/'models'/f"{r['id']}_seed{r['seed']}").relative_to(ROOT))} for r in selected]
    frozen=json.loads((ROOT/'data/final_his_panel/manifest.json').read_text())['structural_preparation']
    assert PROTOCOL==frozen
    record={'protocol':PROTOCOL,'source_results_sha256':sha(source),'script_sha256':sha(Path(__file__)),
            'preparation_library_sha256':sha(ROOT/'scripts/receptor_preparation_v2.py'),
            'source_hashes':{r['source']:sha(ROOT/r['source']) for r in plan},'plan':plan,'rows':[]}
    save(out/'plan.json',record)
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(prepare,item) for item in plan]
        for future in as_completed(futures):
            row=future.result();record['rows'].append(row);save(out/'results.json',record)
            print(len(record['rows']),'/',len(plan),row['id'],row['status'],'movement',row.get('displacement',{}).get('pass'),flush=True)
    assert len(record['rows'])==len(plan)


if __name__=='__main__':main()
