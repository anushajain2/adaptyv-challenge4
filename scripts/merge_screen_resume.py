"""Join disjoint prediction checkpoints without modifying either raw run."""
import hashlib,json,os,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def merge():
    first=ROOT/'outputs/broad009_complete/broad-009-6wrw'
    resume=ROOT/'outputs/screen_resume010_complete/screen-resume-010-6wrw'
    out=ROOT/'outputs/broad009_merged/broad-009-6wrw'
    checkpoint=json.loads((ROOT/'reports/broad_campaign/run_record/interruption_checkpoint.json').read_text())
    status=json.loads((resume/'run_status.json').read_text())
    assert status['statuses']==[{'seed':307,'returncode':0}]
    cases=json.loads((first/'cases.json').read_text());names={c['id'] for c in cases}
    done=set(checkpoint['completed_ids']);missing=set(checkpoint['missing_ids'])
    assert done.isdisjoint(missing) and done|missing==names
    resumed=json.loads((resume/'inputs/manifest.json').read_text())
    assert {c['id'] for c in resumed['cases']}==missing
    pre=json.loads((resume/'seed-307/preflight.json').read_text())
    assert {r['id'] for r in pre}==missing and all(r['binder_template_excluded'] and r['msa_depth']==1 for r in pre)
    assert len(list((resume/'seed-307/fold_out_npz').glob('*.npz')))==len(missing)
    shutil.copytree(first,out,copy_function=os.link)
    copied=[]
    for name in sorted(missing):
        for sub,ext in (('fold_out_npz','.npz'),('refold_cif','.cif')):
            src=resume/'seed-307'/sub/(name+ext);dst=out/'seed-307'/sub/(name+ext)
            assert not dst.exists(), (name,'would overwrite a preserved prediction')
            os.link(src,dst);copied.append({'path':str(dst.relative_to(out)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
    summary=[]
    for case in cases:
        name=case['id'];n=194+len(case['binder_sequence'])
        with np.load(out/'seed-307/fold_out_npz'/f'{name}.npz',allow_pickle=False) as d:
            assert d['pae'].shape==(1,n,n) and d['plddt'].shape==(1,n)
            assert np.isfinite(d['pae']).all() and np.isfinite(d['plddt']).all()
            metrics={k:d[k].tolist() for k in ('iptm','ptm','design_to_target_iptm','interaction_pae')}
        summary.append({'id':name,'seed':307,'selected_sample':0,**metrics})
    assert not (out/'seed-307/confidence_summary.json').exists()
    (out/'seed-307/confidence_summary.json').write_text(json.dumps(summary,indent=2))
    provenance={'original_run':first.name,'resume_run':resume.name,'completed_original_predictions':len(done),
        'completed_resume_predictions':len(missing),'combined_predictions':len(names),
        'prediction_source_run':{n:first.name if n in done else resume.name for n in sorted(names)},
        'copied_files':copied,'note':'Disjoint outputs joined with no selection or overwriting. Original raw runs remain unchanged. Seed 307 restarted for missing cases; its RNG trajectory differs from uninterrupted execution.'}
    (out/'execution_provenance.json').write_text(json.dumps(provenance,indent=2))
    (ROOT/'reports/broad_campaign/run_record/execution_provenance.json').write_text(json.dumps(provenance,indent=2))
    print(out)
if __name__=='__main__':merge()
