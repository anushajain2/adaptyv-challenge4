"""Apply unchanged full-construct preparation once to every selected model."""
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor,as_completed
import prepare_tag_models as tag
from prepare_reserve_tag import COHORT
from novelty_search import ROOT,sha,save
tag.PROTOCOL['cohort']=COHORT


def main():
    source=ROOT/'reports/reserve_tag/results.json';out=ROOT/'reports/reserve_tag/preparation';out.mkdir(exist_ok=False)
    m=json.loads((ROOT/'data/reserve_tag/manifest.json').read_text());assert tag.PROTOCOL==m['preparation']
    selected=[r for r in json.loads(source.read_text())['rows'] if r['confidence_selected']]
    plan=[{'id':r['id'],'seed':r['seed'],'sample':r['sample'],'core_length':r['core_length'],
        'source':r['raw_scoring_output']+'/model.pdb','dest':str((out/'models'/f"{r['id']}_seed{r['seed']}").relative_to(ROOT))} for r in selected]
    record={'protocol':tag.PROTOCOL,'source_results_sha256':sha(source),'implementation_sha256':sha(ROOT/'scripts/prepare_tag_models.py'),
        'source_hashes':{r['source']:sha(ROOT/r['source']) for r in plan},'plan':plan,'rows':[]}
    save(out/'plan.json',record)
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        for future in as_completed([pool.submit(tag.prepare,item) for item in plan]):
            row=future.result();record['rows'].append(row);save(out/'results.json',record)
            print(len(record['rows']),'/',len(plan),row['id'],row['seed'],row['status'],'movement',row.get('displacement',{}).get('pass'),flush=True)
    assert len(record['rows'])==len(plan)


if __name__=='__main__':main()
