"""Evaluate every fresh-seed sample, selecting by confidence before geometry."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from analyze_generation_comparison import design_pose_metrics,contacts_summary,screen_sample
from analyze_controls import pose_metrics
from audit_control_samples import sample_structures
from prepare_controls import ROOT,read_structure,geometry,sequence
from prepare_target import write_csv
from score_matched_controls import run_ipsae,interface_pae


def analyze(run,output):
    output.mkdir(parents=True,exist_ok=True)
    m=json.loads((run/'inputs/manifest.json').read_text());arms={a['id']:a for a in m['arms']}
    for name,sha in m['input_sha256'].items():
        assert hashlib.sha256((ROOT/'data/broad_validation'/name).read_bytes()).hexdigest()==sha
    full=read_structure(ROOT/'data/raw/11QC.cif')
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    rows=[]
    for seed in m['refolding']['seeds']:
        folder=run/f'seed-{seed}';pre=json.loads((folder/'preflight.json').read_text())
        assert len(pre)==len(m['cases']) and all(p['binder_template_excluded'] and p['msa_depth']==1 for p in pre)
        for case in m['cases']:
            name=case['id'];ref=read_structure(ROOT/'data/broad_validation'/f'{name}.pdb')
            with np.load(folder/'fold_out_npz'/f'{name}.npz',allow_pickle=False) as d:
                paes=d['pae'].copy();plddts=d['plddt'].copy();n=194+len(case['binder_sequence'])
                assert paes.shape==(5,n,n) and plddts.shape==(5,n)
                assert np.all(d['asym_id'][0,:194]==0) and np.all(d['asym_id'][0,194:]==1) and d['token_pad_mask'].all()
            for idx,s,metrics,selected in sample_structures(folder,name,restore_sidechains=True):
                assert sequence(s[0]['A'])==m['target_sequence'] and sequence(s[0]['B'])==case['binder_sequence']
                mean=float((paes[idx,:194,194:].mean()+paes[idx,194:,:194].mean())/2)
                assert abs(mean-metrics['interaction_pae'])<.03
                dest=ROOT/'outputs/broad_validation_scoring'/run.name/f'seed-{seed}'/name/str(idx)
                scores=run_ipsae(s,paes[idx],plddts[idx],dest);scores.update(interface_pae(s,paes[idx]))
                scores['binder_mean_plddt']=float(plddts[idx,194:].mean())
                r={'id':name,'arm':case['arm'],'kind':case['kind'],'seed':seed,'sample':idx,
                   'confidence_selected':selected,'metrics':metrics,'scores':scores,'raw_scoring_output':str(dest.relative_to(ROOT))}
                if case['kind']=='novel_design':
                    r.update(pose=design_pose_metrics(s,ref),geometry=geometry(s,{i:i+189 for i in range(1,195)},full),
                             contacts=contacts_summary(s,arms[case['arm']],conserved))
                    r['screen']=screen_sample(r,m['screen'])
                elif case['kind']=='positive':
                    p=pose_metrics(s,ref);r['pose']=p
                    r['control_pass']=(metrics['design_to_target_iptm']>=.5 and p['binder_CA_fold_rmsd_A']<=2 and
                        p['binder_interface_CA_rmsd_after_target_fit_A'] is not None and p['binder_interface_CA_rmsd_after_target_fit_A']<=2 and p['native_contact_recall']>=.5)
                else:r['control_pass']=metrics['design_to_target_iptm']<.5
                rows.append(r)
            print('Scored',seed,name,flush=True)
    assert len(rows)==len(m['cases'])*10
    checks=[r for r in rows if r['arm']=='control' and r['confidence_selected']]
    assert len(checks)==6
    controls_pass=all(r['control_pass'] for r in checks)
    summaries=[]
    for c in m['cases']:
        if c['kind']!='novel_design':continue
        rs=[r for r in rows if r['id']==c['id'] and r['confidence_selected']];assert len(rs)==2
        candidate_pass=all(r['screen']['pass'] for r in rs)
        passed=controls_pass and candidate_pass
        summaries.append({'id':c['id'],'arm':c['arm'],'protocol':c['protocol'],'sequence':c['binder_sequence'],
            'both_seeds_candidate_screen':candidate_pass,'batch_controls_pass':controls_pass,
            'both_seeds_full_screen':passed,'provisional_computational_parent':passed and c['protocol']=='protein-anything',
            'selected_models':rs,'passing_samples_out_of_10':sum(r['screen']['pass'] for r in rows if r['id']==c['id'])})
    result={'run':run.name,'source_run':m['source_run'],'controls_pass':controls_pass,'control_checks':checks,
        'candidate_summaries':summaries,'rows':rows,
        'interpretation':'Repeated structural self-consistency, not experimental validation. Nanobody confidence remains uncalibrated. No affinity, mouse binding or pH switching established.'}
    (output/'results.json').write_text(json.dumps(result,indent=2))
    write_csv(output/'candidates.csv',[{'id':s['id'],'arm':s['arm'],'sequence':s['sequence'],
        'both_seeds_candidate_screen':s['both_seeds_candidate_screen'],'batch_controls_pass':controls_pass,
        'both_seeds_full_screen':s['both_seeds_full_screen'],'provisional_computational_parent':s['provisional_computational_parent'],
        'passing_samples_out_of_10':s['passing_samples_out_of_10'],
        'selected_models_json':json.dumps(s['selected_models'])} for s in summaries])
    (output/'provisional_parents.fasta').write_text(''.join(f">{s['id']}\n{s['sequence']}\n" for s in summaries if s['provisional_computational_parent']))
    print(json.dumps({'controls_pass':controls_pass,'evaluated':len(summaries),
        'provisional_parents':sum(s['provisional_computational_parent'] for s in summaries)},indent=2))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path,default=ROOT/'reports/broad_validation')
    a=p.parse_args();analyze(a.run,a.output)
