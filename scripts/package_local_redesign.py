"""Archive selected models and all confidence samples; independently check pKa records."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from novelty_search import ROOT,sha,save,now
from direction_check import histidines
from analyze_receptor_preparation_v2 import direction_metrics
from summarize_direction import terminal_distance
from run_local_mpnn import eligible_sequence


def main():
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);run=p.parse_args().run.resolve()
    out=ROOT/'reports/local_redesign_panel';raw=json.loads((out/'results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text());direction=json.loads((out/'direction.json').read_text())
    m=json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    expected={(c['id'],seed,i) for c in m['cases'] for seed in m['refolding']['seeds'] for i in range(5)}
    assert {(r['id'],r['seed'],r['sample']) for r in raw['rows']}==expected and len(raw['rows'])==len(expected)
    selected={(r['id'],r['seed'],r['sample']):r for r in raw['rows'] if r['confidence_selected']}
    pm={(r['id'],r['seed'],r['sample']):r for r in prep['rows']};assert set(pm)==set(selected)
    assert prep['source_results_sha256']==sha(out/'results.json')
    bundle=out/'review_bundle';bundle.mkdir(exist_ok=False);models=[];arrays=[]
    for key,r in selected.items():
        pr=pm[key];dest=bundle/f"{r['id']}_seed{r['seed']}";dest.mkdir()
        src=ROOT/r['raw_scoring_output']/'model.pdb';assert sha(src)==pr['source_sha256']
        shutil.copyfile(src,dest/'raw.pdb');save(dest/'raw_checks.json',r)
        if pr['status']=='ok':
            src=ROOT/pr['dest']/'prepared.pdb';assert sha(src)==pr['prepared_sha256'];shutil.copyfile(src,dest/'prepared.pdb')
        models.append({'id':r['id'],'seed':r['seed'],'sample':r['sample'],'files':{f.name:sha(f) for f in dest.iterdir()}})
    arraydir=bundle/'all_sample_confidence';arraydir.mkdir()
    for c in m['cases']:
        for seed in m['refolding']['seeds']:
            src=run/f'seed-{seed}/fold_out_npz'/f"{c['id']}.npz";dest=arraydir/f"{c['id']}_seed{seed}.npz";shutil.copyfile(src,dest)
            with np.load(dest,allow_pickle=False) as d:
                n=194+len(c['binder_sequence']);assert d['pae'].shape==(5,n,n) and np.isfinite(d['pae']).all()
            arrays.append({'path':str(dest.relative_to(ROOT)),'sha256':sha(dest)})
    count=0
    for r in direction['rows']:
        if r['direction']['status']!='diagnostic':continue
        path=ROOT/r['source'];assert sha(path)==r['source_sha256'];positions=cases[r['id']]['histidine_positions']
        free=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',positions)
        bound=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',positions)
        check=direction_metrics(free,bound,terminal_distance(path,positions))
        for key,value in check.items():assert r['direction'][key]==value
        count+=1
    families={f['id']:f for f in m['prospective_plan']['families']}
    for c in m['cases']:
        if c.get('local_redesign'):assert eligible_sequence(c['binder_sequence'],families[c['parent_id']])
    assert len(direction['followup_ids'])<=2 and all(v['fully_qualified_for_followup'] for v in direction['variants'] if v['id'] in direction['followup_ids'])
    fallback=ROOT/'reports/submission_release_2026-10-25'
    frozen=json.loads((fallback/'manifest.json').read_text())
    for name,digest in frozen['files'].items():assert sha(fallback/name)==digest
    assert sha(ROOT/'reports/submission_draft/submission_draft.csv')==sha(fallback/'submission_draft.csv')
    save(bundle/'manifest.json',{'raw_sha256':sha(out/'results.json'),'preparation_sha256':sha(out/'preparation/results.json'),
        'direction_sha256':sha(out/'direction.json'),'models':models,'all_sample_confidence':arrays})
    save(out/'evidence_audit.json',{'status':'PASS','checked_utc':now(),'predictions':len(expected),'selected_models':len(selected),
        'full_confidence_files':len(arrays),'pka_records_recomputed':count,'fixed_positions_and_target_preserved':True,'frozen_release_and_submission_csv_unchanged':True,
        'bundle_manifest_sha256':sha(bundle/'manifest.json'),'followup_ids':direction['followup_ids']})
    if not direction['followup_ids']:
        save(out/'conditional_followup.json',{'status':'SKIPPED','reason':'No prospectively qualifying improvement','spent_usd':0,'followup_reservation_usd':3})
    print('Archived and audited',len(selected),'selected models,',len(arrays),'complete confidence files,',count,'pKa records.')


if __name__=='__main__':main()
