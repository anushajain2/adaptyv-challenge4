"""Audit all saved predictions and package the selected evidence for review."""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from prepare_controls import ROOT,read_structure,sequence
from analyze_tag_comparison import core_arrays,cross_pae,select_core
from receptor_preparation_v2 import sha


def main():
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);args=p.parse_args();run=args.run.resolve()
    out=ROOT/'reports/tag_comparison';m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text());result=json.loads((out/'results.json').read_text())
    assert len(result['rows'])==140 and len({(r['id'],r['seed'],r['sample']) for r in result['rows']})==140
    assert result['manifest_sha256']==sha(run/'inputs/manifest.json')
    launch=json.loads((out/'run_record/launch_ledger.json').read_text())
    for name,digest in launch['input_sha256'].items():assert sha(run/'inputs'/name)==digest
    for name,digest in launch['analysis_sha256'].items():assert sha(out/'run_record'/name)==digest
    assert sha(ROOT/'scripts/analyze_tag_comparison.py')==launch['analysis_sha256']['analyze_tag_comparison.py']
    bundle=out/'review_bundle';bundle.mkdir(exist_ok=False)
    inventory=[];models=[]
    for seed in m['refolding']['seeds']:
        folder=run/f'seed-{seed}';pre={r['id']:r for r in json.loads((folder/'preflight.json').read_text())}
        assert len(pre)==14
        for case in m['cases']:
            name=case['id'];length=case['core_length'];source=folder/'fold_out_npz'/f'{name}.npz'
            assert pre[name]['binder_tokens']==len(case['binder_sequence']) and pre[name]['target_tokens']==194
            assert pre[name]['binder_template_excluded'] and pre[name]['msa_depth']==1
            rows=[r for r in result['rows'] if r['id']==name and r['seed']==seed]
            selected=[r for r in rows if r['confidence_selected']]
            chosen=select_core(rows);assert len(selected)==int(chosen is not None)
            with np.load(source,allow_pickle=False) as d:
                for r in rows:
                    pae,plddt=core_arrays(d['pae'][r['sample']],d['plddt'][r['sample']],length)
                    assert cross_pae(pae,length)==r['core_mean_cross_PAE_A']
                    s=read_structure(ROOT/r['raw_scoring_output']/'model.pdb')
                    assert sequence(s[0]['A'])==m['target_sequence'] and sequence(s[0]['B'])==case['binder_sequence']
                if selected:
                    r=selected[0];assert r['sample']==chosen
                    dest=bundle/f'{name}_seed{seed}';dest.mkdir()
                    shutil.copyfile(ROOT/r['raw_scoring_output']/'model.pdb',dest/'raw_complex.pdb')
                    np.savez_compressed(dest/'prediction_confidence.npz',pae=d['pae'][chosen],plddt=d['plddt'][chosen],core_length=np.array(length),target_length=np.array(194))
                    (dest/'metrics.json').write_text(json.dumps(r,indent=2))
                    models.append({'id':name,'seed':seed,'sample':chosen,'files':{p.name:sha(p) for p in dest.iterdir()}})
            inventory.append({'id':name,'seed':seed,'source':str(source.relative_to(ROOT)),'sha256':sha(source),'samples':len(rows)})
    record={'status':'pass','predictions_audited':140,'selected_models':len(models),'source_files':inventory,
        'checks':['Input and analysis snapshot hashes match','All target/core/suffix sequences and lengths verified','Binder and suffix excluded from target templates','Core scores recomputed from exported PAE','Selection uses core PAE and target integrity only','All 140 raw models retained without substitute samples'],
        'run_status':json.loads((run/'run_status.json').read_text()),'result_sha256':sha(out/'results.json')}
    (out/'prediction_audit.json').write_text(json.dumps(record,indent=2))
    (bundle/'manifest.json').write_text(json.dumps({'models':models,'scope':'Selected raw complexes and full confidence arrays; controls are not submission candidates. No experimental validation.'},indent=2))
    print('Audited 140 predictions; packaged',len(models),'selected complexes')

if __name__=='__main__':main()
