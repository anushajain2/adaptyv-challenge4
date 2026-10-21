"""Class-specific ranking of the one-sample screen; never claim qualified parents."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from analyze_generation_comparison import design_pose_metrics,contacts_summary,screen_sample,aligned_identity
from analyze_controls import pose_metrics
from audit_control_samples import sample_structures
from prepare_controls import ROOT,read_structure,geometry,sequence
from prepare_target import write_csv
from target_integrity_selection import target_ca_rmsd
from read_pilot import KD
from score_matched_controls import run_ipsae,interface_pae


def rank_candidates(rows):
    """Higher equal-weight within-class percentile is better; ties average."""
    for protocol in sorted({r['protocol'] for r in rows}):
        group=[r for r in rows if r['protocol']==protocol]
        metrics=[[-r['metrics']['design_to_target_iptm'] for r in group],
                 [r['scores']['contact_pae_both_directions_A'] if r['scores']['contact_pae_both_directions_A'] is not None else 1e6 for r in group],
                 [-r['hotspot_fraction'] for r in group],
                 [r['pose']['binder_CA_fold_rmsd_A'] for r in group]]
        percentiles=np.array([1-(rankdata(v,method='average')-1)/max(1,len(group)-1) for v in metrics]).mean(axis=0)
        for r,p in zip(group,percentiles):r['within_class_rank_score']=float(p)
    return sorted(rows,key=lambda r:(not r['geometry_eligible'],-r['within_class_rank_score'],r['id']))


def diverse_shortlist(ranked,limit=30):
    """Geometry-eligible only. Round robin strata, avoiding >80% sequence copies."""
    selected=[]
    pools={k:[r for r in ranked if r['arm']==k and r['geometry_eligible'] and not r.get('matches_known_sequence',False)] for k in sorted({r['arm'] for r in ranked})}
    while len(selected)<limit and any(pools.values()):
        for pool in pools.values():
            while pool:
                r=pool.pop(0)
                if all(r['sequence']!=s['sequence'] and (
                    r['protocol']!=s['protocol'] or aligned_identity(
                        r.get('diversity_sequence',r['sequence']),
                        s.get('diversity_sequence',s['sequence']))<.8) for s in selected):
                    selected.append(r);break
            if len(selected)>=limit:break
    return selected


def analyze(run,output):
    output.mkdir(parents=True,exist_ok=True)
    m=json.loads((run/'inputs/manifest.json').read_text());cases=json.loads((run/'cases.json').read_text())
    provenance=json.loads((run/'execution_provenance.json').read_text()) if (run/'execution_provenance.json').exists() else {}
    arms={a['id']:a for a in m['arms']};seed=m['refolding']['seeds'][0];folder=run/f'seed-{seed}'
    assert m['refolding']['diffusion_samples_per_seed']==1
    preflight=json.loads((folder/'preflight.json').read_text())
    expected_designs=sum(a['num_designs'] for a in m['arms'])
    assert len(preflight)==expected_designs+len(m['controls']) and all(p['binder_template_excluded'] and p['msa_depth']==1 for p in preflight)
    feature_audit=[]
    for arm in arms.values():
        paths=sorted((run/arm['id']/'intermediate_designs').glob('*.npz'));assert len(paths)==arm['num_designs']
        for path in paths:
            with np.load(path,allow_pickle=False) as d:
                hot=(np.flatnonzero(d['binding_type'][:194]==1)+1).tolist()
                assert hot==arm['hotspots_crop']
                assert not d['design_mask'][:194].any() and d['design_mask'][194:].any()
            feature_audit.append({'path':str(path.relative_to(run)),'hotspots_crop':hot,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    full=read_structure(ROOT/'data/raw/11QC.cif')
    known_sequences={sequence(read_structure(ROOT/path)[0]['B']) for path in (
        'data/controls_matched/2ds25_5_6wrw.pdb','data/controls_matched/3ds18_6wrv.pdb','configs/7eow.cif')}
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((run/'inputs/spatial_map.csv').open())}
    rows=[];controls=[]
    for case in cases:
        name=case['id'];ref=read_structure(run/'references'/f'{name}.pdb')
        with np.load(folder/'fold_out_npz'/f'{name}.npz',allow_pickle=False) as d:
            paes=d['pae'].copy();plddts=d['plddt'].copy();n=194+len(case['binder_sequence'])
            assert paes.shape==(1,n,n) and plddts.shape==(1,n)
            assert np.all(d['asym_id'][0,:194]==0) and np.all(d['asym_id'][0,194:]==1) and d['token_pad_mask'].all()
        samples=list(sample_structures(folder,name,restore_sidechains=True));assert len(samples)==1
        index,s,metrics,selected=samples[0];assert index==0 and selected
        assert sequence(s[0]['A'])==m['target_sequence'] and sequence(s[0]['B'])==case['binder_sequence']
        avg=float((paes[0,:194,194:].mean()+paes[0,194:,:194].mean())/2)
        assert abs(avg-metrics['interaction_pae'])<.03
        dest=ROOT/'outputs/broad_scoring'/run.name/name
        scores=run_ipsae(s,paes[0],plddts[0],dest);scores.update(interface_pae(s,paes[0]));scores['binder_mean_plddt']=float(plddts[0,194:].mean())
        row={'target_CA_rmsd_to_input_A':target_ca_rmsd(s,ref),'id':name,'arm':case['arm'],'kind':case['kind'],'seed':seed,'sample':0,
             'metrics':metrics,'scores':scores,'raw_scoring_output':str(dest.relative_to(ROOT)),
             'prediction_source_run':provenance.get('prediction_source_run',{}).get(name,run.name)}
        if case['kind']!='novel_design':
            if case['kind']=='positive':
                pose=pose_metrics(s,ref);row['pose']=pose
                row['pass']=(metrics['design_to_target_iptm']>=.5 and pose['binder_CA_fold_rmsd_A']<=2 and
                    pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and pose['binder_interface_CA_rmsd_after_target_fit_A']<=2 and pose['native_contact_recall']>=.5)
            else:row['pass']=metrics['design_to_target_iptm']<.5
            row['target_integrity_pass']=row['target_CA_rmsd_to_input_A']<=m['selection_rule']['max_target_CA_rmsd_to_input_A'];row['pass']=row['pass'] and row['target_integrity_pass'];controls.append(row)
        else:
            arm=arms[case['arm']];seq=case['binder_sequence']
            row.update(protocol=case['protocol'],epitope=case['epitope'],sequence=seq,length=len(seq),
                mean_kyte_doolittle=sum(KD[a] for a in seq)/len(seq),
                geometry=geometry(s,{i:i+189 for i in range(1,195)},full),pose=design_pose_metrics(s,ref),
                contacts=contacts_summary(s,arm,conserved))
            row['hotspot_fraction']=len(row['contacts']['hotspots_contacted_crop'])/len(arm['hotspots_crop'])
            row['screen']=screen_sample(row,m['screen'])
            row['screen']['checks']['target_integrity']=row['target_CA_rmsd_to_input_A']<=m['selection_rule']['max_target_CA_rmsd_to_input_A']
            row['screen']['failed']=[k for k,v in row['screen']['checks'].items() if not v]
            row['screen']['pass']=all(row['screen']['checks'].values())
            row['geometry_eligible']=all(v for k,v in row['screen']['checks'].items() if k!='interface_confidence')
            row['single_sample_full_screen']=row['screen']['pass']
            row['qualified_parent']=False
            row['matches_known_sequence']=seq in known_sequences
            metadata=run/case['generated_source'].replace('intermediate_designs_inverse_folded','intermediate_designs')
            with np.load(metadata.with_suffix('.npz'),allow_pickle=False) as generated:
                designed=(np.flatnonzero(generated['design_mask'][194:])+1).tolist()
            assert designed and max(designed)<=len(seq)
            row['designed_binder_positions']=designed
            row['diversity_sequence']=''.join(seq[i-1] for i in designed) if case['protocol']=='nanobody-anything' else seq
            rows.append(row)
        print('Scored',name,flush=True)
    assert len(rows)==expected_designs and len(controls)==len(m['controls'])==3
    controls_pass=all(r['pass'] for r in controls)
    prior=json.loads((ROOT/m['prior_campaign']).read_text())
    old_sequences={r['sequence'] for r in prior['ranked_candidates']}
    old_leads=[r['sequence'] for r in prior['ranked_candidates'] if r['shortlisted_for_deep_evaluation']]
    for r in rows:
        r['exact_prior_sequence']=r['sequence'] in old_sequences
        r['similar_to_prior_shortlist']=any(aligned_identity(r['sequence'],seq)>=.8 for seq in old_leads)
        r['matches_known_sequence']=r['matches_known_sequence'] or r['exact_prior_sequence'] or r['similar_to_prior_shortlist']
    ranked=rank_candidates(rows)
    shortlist=diverse_shortlist(ranked) if controls_pass else []
    shortlisted={r['id'] for r in shortlist}
    for i,r in enumerate(ranked):r['screening_order']=i+1;r['shortlisted_for_deep_evaluation']=r['id'] in shortlisted
    summaries={a:{'designs':sum(r['arm']==a for r in rows),
        'geometry_eligible':sum(r['arm']==a and r['geometry_eligible'] for r in rows),
        'single_sample_full_screen':sum(r['arm']==a and r['single_sample_full_screen'] for r in rows),
        'max_interface_iptm':max(r['metrics']['design_to_target_iptm'] for r in rows if r['arm']==a),
        'shortlisted':sum(r['arm']==a for r in shortlist)} for a in arms}
    result={'run':run.name,'controls_pass':controls_pass,'controls':controls,'arms':summaries,'ranked_candidates':ranked,
        'execution_provenance':provenance,
        'unique_novel_sequences':len({r['sequence'] for r in rows}),
        'known_sequence_check_scope':'Exact identity to two published positive controls and the starting VHH framework only; not a comprehensive novelty search.',
        'shortlisted_ids':[r['id'] for r in shortlist],'actual_generated_feature_audit':feature_audit,
        'interpretation':'One-sample screening only. No experimentally validated binders or qualified parents. Geometry uses generated backbone references, not experimental poses. Nanobody ipTM threshold is descriptive, not calibrated. Scores are relative within binder class; no cross-class affinity ordering.',
        'manifest_sha256':hashlib.sha256((run/'inputs/manifest.json').read_bytes()).hexdigest()}
    (output/'results.json').write_text(json.dumps(result,indent=2))
    write_csv(output/'ranked_designs.csv',[{'id':r['id'],'arm':r['arm'],'sequence':r['sequence'],'length':r['length'],
        'geometry_eligible':r['geometry_eligible'],'single_sample_full_screen':r['single_sample_full_screen'],
        'shortlisted':r['id'] in shortlisted,'within_class_rank_score':r['within_class_rank_score'],
        'interface_iptm':r['metrics']['design_to_target_iptm'],'contact_pae_A':r['scores']['contact_pae_both_directions_A'],
        'ipsae_min':r['scores']['ipsae_min'],'binder_mean_plddt':r['scores']['binder_mean_plddt'],
        'mean_kyte_doolittle':r['mean_kyte_doolittle'],'hotspot_fraction':r['hotspot_fraction'],
        'binder_fold_rmsd_A':r['pose']['binder_CA_fold_rmsd_A'],'failed_checks':';'.join(r['screen']['failed'])} for r in ranked])
    for name,subset in [('all_designs',ranked),('shortlist_for_evaluation',shortlist)]:
        (output/f'{name}.fasta').write_text(''.join(f">{r['id']}\n{r['sequence']}\n" for r in subset))
    print(json.dumps({'controls_pass':controls_pass,'arms':summaries,'shortlist':len(shortlist)},indent=2))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--output',type=Path,default=ROOT/'reports/apical_expansion')
    a=p.parse_args();analyze(a.run,a.output)
