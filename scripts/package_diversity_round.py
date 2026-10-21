"""Archive all selected confirmation poses and their preparation caveats."""
import json
import shutil
from novelty_search import ROOT, save, sha, now
from build_submission_draft import sequence_audit


def main():
    out=ROOT/'reports/diversity_confirmation'
    raw=json.loads((out/'results.json').read_text())
    prep=json.loads((out/'preparation/results.json').read_text())
    m=json.loads((ROOT/'data/diversity_confirmation/manifest.json').read_text())
    selected=[r for r in raw['rows'] if r['confidence_selected']]
    pm={(r['id'],r['seed']):r for r in prep['rows']}
    assert {(r['id'],r['seed']) for r in selected}==set(pm)
    bundle=out/'review_bundle';bundle.mkdir(exist_ok=False)
    models=[]
    for r in selected:
        p=pm[r['id'],r['seed']];dest=bundle/f"{r['id']}_seed{r['seed']}";dest.mkdir()
        source=ROOT/r['raw_scoring_output']/'model.pdb'
        assert sha(source)==p['source_sha256']
        shutil.copyfile(source,dest/'raw.pdb')
        save(dest/'raw_checks.json',r);save(dest/'preparation.json',p)
        if p['status']=='ok':
            prepared=ROOT/p['dest']/'prepared.pdb';assert sha(prepared)==p['prepared_sha256']
            shutil.copyfile(prepared,dest/'prepared.pdb')
        models.append({'id':r['id'],'seed':r['seed'],'sample':r['sample'],
            'raw_pass':r.get('screen',{}).get('pass',r.get('control_pass')),
            'source':str(source.relative_to(ROOT)),'source_sha256':sha(source),
            'files':{f.name:sha(f) for f in dest.iterdir()}})
    save(bundle/'manifest.json',{'source_results_sha256':sha(out/'results.json'),
        'preparation_results_sha256':sha(out/'preparation/results.json'),'models':models})
    candidates=[]
    for c in raw['candidate_summaries']:
        ps=[pm[c['id'],seed] for seed in m['refolding']['seeds'] if (c['id'],seed) in pm]
        prep_pass=len(ps)==2 and all(p['status']=='ok' and p['displacement']['pass']
                                    and p['prepared_internal_geometry']['pass'] for p in ps)
        candidates.append({'id':c['id'],'sequence':c['sequence'],
            'two_seed_raw_screen_pass':c['all_selected_seeds_full_screen'],
            'passing_samples_of_10':c['passing_samples'],'two_pose_preparation_pass':prep_pass,
            'sequence_audit':sequence_audit(c['sequence']),
            'selected_models':[{'seed':r['seed'],'sample':r['sample'],
                'path':str((bundle/f"{r['id']}_seed{r['seed']}"/'raw.pdb').relative_to(ROOT)),
                'sha256':sha(bundle/f"{r['id']}_seed{r['seed']}"/'raw.pdb'),
                'raw_pass':r['screen']['pass'],'failed_checks':r['screen']['failed']} for r in c['selected_models']],
            'experimental_binding':'untested','pH_release':'unestablished','mouse_cross_reactivity':'unestablished',
            'assay_tag_compatibility':'not modeled for these new sequences',
            'database_review':'pending for any retained new candidate','submission_admitted':False})
    summary={'created_utc':now(),'controls_pass':raw['controls_pass'],'candidates':candidates,
        'raw_two_seed_pass_count':sum(c['two_seed_raw_screen_pass'] for c in candidates),
        'preparation_completed':len(prep['rows']),
        'preparation_execution_ok':sum(p['status']=='ok' for p in prep['rows']),
        'preparation_movement_pass':sum(p.get('displacement',{}).get('pass',False) for p in prep['rows']),
        'submission_modified':False,'source_results_sha256':sha(out/'results.json'),
        'review_bundle_manifest_sha256':sha(bundle/'manifest.json')}
    save(out/'summary.json',summary)
    print('Packaged',len(models),'selected models;',summary['raw_two_seed_pass_count'],'raw two-seed passes.')


if __name__=='__main__':main()
