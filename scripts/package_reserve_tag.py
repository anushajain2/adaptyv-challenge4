"""Preserve complete confidence arrays, selected constructs and tag-contact diagnostics."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from novelty_search import ROOT,sha,save,now
from prepare_controls import read_structure,sequence,save_pdb
from analyze_tag_comparison import core_structure
from analyze_controls import pose_metrics
from audit_tag_contacts import close_pairs
from direction_check import histidines
from analyze_receptor_preparation_v2 import direction_metrics
from summarize_direction import terminal_distance


def main():
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);run=parser.parse_args().run.resolve()
    out=ROOT/'reports/reserve_tag';raw=json.loads((out/'results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text());direction=json.loads((out/'direction.json').read_text())
    m=json.loads((ROOT/'data/reserve_tag/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    expected={(c['id'],seed,i) for c in m['cases'] for seed in m['refolding']['seeds'] for i in range(5)}
    assert len(raw['rows'])==len(expected)==140 and {(r['id'],r['seed'],r['sample']) for r in raw['rows']}==expected
    selected={(r['id'],r['seed'],r['sample']):r for r in raw['rows'] if r['confidence_selected']}
    pm={(p['id'],p['seed'],p['sample']):p for p in prep['rows']};assert set(selected)==set(pm)
    assert prep['source_results_sha256']==sha(out/'results.json')
    bundle=out/'review_bundle';bundle.mkdir(exist_ok=False);arraydir=bundle/'all_sample_confidence';arraydir.mkdir()
    arrays=[];models=[];tags=[]
    pair_lookup={(r['base_id'],r['condition'],r['seed']):r for r in selected.values()}
    for c in m['cases']:
        for seed in m['refolding']['seeds']:
            src=run/f'seed-{seed}/fold_out_npz'/f"{c['id']}.npz";dest=arraydir/f"{c['id']}_seed{seed}.npz";shutil.copyfile(src,dest)
            with np.load(dest,allow_pickle=False) as data:
                n=194+len(c['binder_sequence']);assert data['pae'].shape==(5,n,n) and np.isfinite(data['pae']).all()
                assert data['plddt'].shape==(5,n)
            arrays.append({'path':str(dest.relative_to(ROOT)),'sha256':sha(dest)})
    for key,r in selected.items():
        p=pm[key];c=cases[r['id']];dest=bundle/f"{r['id']}_seed{r['seed']}";dest.mkdir()
        src=ROOT/r['raw_scoring_output']/'model.pdb';assert sha(src)==p['source_sha256']
        shutil.copyfile(src,dest/'raw_complex.pdb');save(dest/'raw_checks.json',r)
        structure=read_structure(src);core=core_structure(structure,c['core_length']);save_pdb(core,dest/'raw_core.pdb')
        assert sequence(structure[0]['B'])==c['binder_sequence']
        if p['status']=='ok':
            prepared=ROOT/p['dest']/'prepared.pdb';assert sha(prepared)==p['prepared_sha256']
            assert sequence(read_structure(ROOT/p['dest']/'free_binder.pdb')[0]['B'])==c['binder_sequence'],'Suffix lost in free-form chemistry'
            shutil.copyfile(prepared,dest/'prepared_complex.pdb');save_pdb(core_structure(read_structure(prepared),c['core_length']),dest/'prepared_core.pdb')
        models.append({'id':r['id'],'seed':r['seed'],'sample':r['sample'],'files':{x.name:sha(x) for x in dest.iterdir()}})
        if c['condition']=='tagged':
            other=pair_lookup.get((c['base_id'],'bare',r['seed']));paired=None
            if other:
                barecore=core_structure(read_structure(ROOT/other['raw_scoring_output']/'model.pdb'),c['core_length'])
                paired=pose_metrics(core,barecore)
            pairs=r['tag']['target_contact_residue_pairs'];pae_value=None
            with np.load(arraydir/f"{r['id']}_seed{r['seed']}.npz",allow_pickle=False) as data:
                pae=data['pae'][r['sample']];plddt=data['plddt'][r['sample']]
                if pairs:
                    aa=np.array([t-1 for t,b in pairs]);bb=np.array([194+b-1 for t,b in pairs])
                    pae_value=float((pae[aa,bb].mean()+pae[bb,aa].mean())/2)
                tag_plddt=float(plddt[194+c['core_length']:].mean()/(100 if plddt.max()>1 else 1))
            tags.append({'id':r['id'],'base_id':r['base_id'],'seed':r['seed'],'core_clashes':close_pairs(core),
                'tag_contacts':pairs,'tag_contact_PAE_A':pae_value,'tag_mean_plddt_0_to_1':tag_plddt,'core_vs_matched_bare_pose':paired})
    overlap=[]
    for name in sorted({r['id'] for r in tags}):
        rs=[r for r in tags if r['id']==name]
        if len(rs)!=2:continue
        a,b=({tuple(pair) for pair in r['tag_contacts']} for r in rs)
        overlap.append({'id':name,'intersection':len(a&b),'union':len(a|b),'jaccard':len(a&b)/len(a|b) if a|b else None})
    save(out/'tag_contact_audit.json',{'models':tags,'across_seed_contacts':overlap,
        'interpretation':'Descriptive only; low-confidence contacts or absence of contacts cannot prove capture-surface compatibility. No gates changed.'})
    count=0
    for r in direction['rows']:
        if r['direction']['status']!='diagnostic':continue
        assert cases[r['id']]['histidine_positions']==[41]
        path=ROOT/r['source'];assert sha(path)==r['source_sha256']
        f=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',[41]);b=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',[41])
        check=direction_metrics(f,b,terminal_distance(path,[41]))
        assert all(r['direction'][key]==value for key,value in check.items());count+=1
    fallback=ROOT/'reports/submission_release_2026-10-25';frozen=json.loads((fallback/'manifest.json').read_text())
    for name,digest in frozen['files'].items():assert sha(fallback/name)==digest
    assert sha(fallback/'submission_draft.csv')==sha(ROOT/'reports/submission_draft/submission_draft.csv')
    save(bundle/'manifest.json',{'raw_sha256':sha(out/'results.json'),'preparation_sha256':sha(out/'preparation/results.json'),
        'direction_sha256':sha(out/'direction.json'),'models':models,'all_sample_confidence':arrays})
    save(out/'evidence_audit.json',{'status':'PASS','checked_utc':now(),'predictions':140,'selected_models':len(selected),
        'full_confidence_arrays':len(arrays),'pka_records_recomputed':count,'full_suffix_chemistry_verified':True,
        'frozen_release_and_submission_csv_unchanged':True,'bundle_manifest_sha256':sha(bundle/'manifest.json')})
    print('Archived',len(models),'models;',len(arrays),'complete confidence arrays;',count,'pKa comparisons.')


if __name__=='__main__':main()
