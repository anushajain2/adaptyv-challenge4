"""Archive selected poses, full per-case confidence arrays and final exclusions."""
import argparse,json,shutil
from pathlib import Path
from novelty_search import ROOT,sha,save,now

def main():
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);run=p.parse_args().run.resolve()
    out=ROOT/'reports/final_his_panel';raw=json.loads((out/'results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text());direction=json.loads((ROOT/'reports/final_his_direction/results.json').read_text())
    pm={(r['id'],r['seed'],r['sample']):r for r in prep['rows']}
    bundle=out/'review_bundle';bundle.mkdir(exist_ok=False);models=[]
    for r in raw['rows']:
        if not r['confidence_selected']:continue
        p=pm[r['id'],r['seed'],r['sample']];dest=bundle/f"{r['id']}_seed{r['seed']}";dest.mkdir()
        source=ROOT/r['raw_scoring_output']/'model.pdb';assert sha(source)==p['source_sha256']
        shutil.copyfile(source,dest/'raw.pdb');save(dest/'raw_checks.json',r);save(dest/'preparation.json',p)
        if p['status']=='ok':
            prepared=ROOT/p['dest']/'prepared.pdb';assert sha(prepared)==p['prepared_sha256'];shutil.copyfile(prepared,dest/'prepared.pdb')
        models.append({'id':r['id'],'seed':r['seed'],'sample':r['sample'],'files':{f.name:sha(f) for f in dest.iterdir()}})
    arrays=bundle/'all_sample_confidence';arrays.mkdir()
    manifest=json.loads((run/'inputs/manifest.json').read_text());confidence=[]
    for seed in manifest['refolding']['seeds']:
        for c in manifest['cases']:
            source=run/f'seed-{seed}/fold_out_npz'/f"{c['id']}.npz"
            dest=arrays/f"{c['id']}_seed{seed}.npz";shutil.copyfile(source,dest)
            confidence.append({'path':str(dest.relative_to(ROOT)),'sha256':sha(dest),'source':str(source.relative_to(ROOT))})
    save(bundle/'manifest.json',{'source_results_sha256':sha(out/'results.json'),'preparation_results_sha256':sha(out/'preparation/results.json'),'direction_results_sha256':sha(ROOT/'reports/final_his_direction/results.json'),'models':models,'all_sample_confidence':confidence})
    save(out/'completion.json',{'created_utc':now(),'predictions':len(raw['rows']),'new_variants':12,'raw_controls_pass':raw['controls_pass'],'raw_variant_passes':sum(v['raw_pass_both_seeds'] for v in direction['variants']),'prepared_controls_pass':direction['prepared_controls_pass'],'direction_passes':sum(v['direction_pass_both_seeds'] for v in direction['variants']),'fully_qualified_passes':sum(v['fully_qualified_direction_pass_both_seeds'] for v in direction['variants']),'tag_selected_ids':direction['tag_selected_ids'],'challenge_submitted':False,'bundle_manifest_sha256':sha(bundle/'manifest.json')})
    print('Archived',len(models),'selected models and',len(confidence),'complete five-sample confidence files.')
if __name__=='__main__':main()
