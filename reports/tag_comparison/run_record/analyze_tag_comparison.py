"""Matched comparison with confidence and interface measures restricted to the core."""
import argparse,copy,csv,json
from pathlib import Path
import numpy as np
from prepare_controls import ROOT,read_structure,sequence,geometry,save_pdb
from analyze_controls import contact_pairs,pose_metrics
from analyze_generation_comparison import design_pose_metrics,contacts_summary,screen_sample
from score_matched_controls import run_ipsae,interface_pae
from target_integrity_selection import target_ca_rmsd
from audit_control_samples import sample_structures
from receptor_preparation_v2 import sha


def core_structure(s,length):
    result=copy.deepcopy(s)
    for r in list(result[0]['B']):
        if r.id[1]>length:result[0]['B'].detach_child(r.id)
    assert len(result[0]['B'])==length
    return result


def core_arrays(pae,plddt,length,n_target=194):
    stop=n_target+length
    assert pae.ndim==2 and pae.shape[0]==pae.shape[1]==len(plddt) and stop<=len(plddt)
    return pae[:stop,:stop].copy(),plddt[:stop].copy()


def cross_pae(pae,length,n_target=194):
    assert pae.shape==(n_target+length,n_target+length)
    return float((pae[:n_target,n_target:].mean()+pae[n_target:,:n_target].mean())/2)


def select_core(rows,limit=2.):
    eligible=[r for r in rows if np.isfinite(r['target_CA_rmsd_to_input_A']) and r['target_CA_rmsd_to_input_A']<=limit and np.isfinite(r['core_mean_cross_PAE_A'])]
    return min(eligible,key=lambda r:(r['core_mean_cross_PAE_A'],r['sample']))['sample'] if eligible else None


def core_screen(row,manifest):
    result=screen_sample(row,manifest['geometry_screen'])
    checks=result['checks'];checks.pop('interface_confidence')
    checks['core_cross_PAE']=row['core_mean_cross_PAE_A']<manifest['core_confidence']['max_mean_cross_PAE_A']
    value=row['scores']['contact_pae_both_directions_A']
    checks['core_contact_PAE']=value is not None and value<=manifest['core_confidence']['max_contact_PAE_A']
    return {'pass':all(checks.values()),'checks':checks,'failed':[k for k,v in checks.items() if not v]}


def main():
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);args=p.parse_args();run=args.run.resolve()
    out=ROOT/'reports/tag_comparison';m=json.loads((run/'inputs/manifest.json').read_text())
    source=ROOT/'data/tag_comparison'
    assert sha(run/'inputs/manifest.json')==sha(source/'manifest.json')
    for name,digest in m['input_sha256'].items():assert sha(source/name)==digest
    full=read_structure(ROOT/'data/raw/11QC.cif');arms={a['id']:a for a in m['arms']}
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    rows=[]
    for seed in m['refolding']['seeds']:
        folder=run/f'seed-{seed}';pre=json.loads((folder/'preflight.json').read_text())
        assert len(pre)==14 and all(x['binder_template_excluded'] and x['msa_depth']==1 for x in pre)
        for c in m['cases']:
            name=c['id'];length=c['core_length'];ref=read_structure(ROOT/c['reference'])
            with np.load(folder/'fold_out_npz'/f'{name}.npz',allow_pickle=False) as d:
                paes=d['pae'].copy();plddts=d['plddt'].copy();n=194+len(c['binder_sequence'])
                assert paes.shape==(5,n,n) and plddts.shape==(5,n)
                assert np.all(d['asym_id'][0,:194]==0) and np.all(d['asym_id'][0,194:]==1) and d['token_pad_mask'].all()
            group=[]
            for idx,s,metrics,upstream_selected in sample_structures(folder,name,restore_sidechains=True):
                assert sequence(s[0]['A'])==m['target_sequence'] and sequence(s[0]['B'])==c['binder_sequence']
                assert abs(cross_pae(paes[idx],len(c['binder_sequence']))-metrics['interaction_pae'])<.03
                core=core_structure(s,length);cp,cl=core_arrays(paes[idx],plddts[idx],length)
                dest=ROOT/'outputs/tag_comparison_scoring'/f'seed-{seed}'/name/str(idx)
                scores=run_ipsae(core,cp,cl,dest/'core');scores.update(interface_pae(core,cp))
                scores['binder_core_mean_plddt']=float(cl[194:].mean())
                for residue,value in zip(list(s[0]['A'])+list(s[0]['B']),plddts[idx]):
                    for atom in residue:atom.bfactor=float(value*100 if max(plddts[idx])<=1 else value)
                save_pdb(s,dest/'model.pdb')
                r={'id':name,'base_id':c['base_id'],'condition':c['condition'],'kind':c['kind'],'arm':c['arm'],
                    'seed':seed,'sample':idx,'confidence_selected':False,'upstream_confidence_selected':upstream_selected,
                    'metrics':metrics,'scores':scores,'core_length':length,'core_mean_cross_PAE_A':cross_pae(cp,length),
                    'target_CA_rmsd_to_input_A':target_ca_rmsd(core,ref),'raw_scoring_output':str(dest.relative_to(ROOT)),
                    'geometry':geometry(core,{i:i+189 for i in range(1,195)},full)}
                if c['kind']=='novel_design':
                    r.update(pose=design_pose_metrics(core,ref),contacts=contacts_summary(core,arms[c['arm']],conserved))
                    r['screen']=core_screen(r,m)
                elif c['kind']=='positive':
                    pose=pose_metrics(core,ref);r['pose']=pose
                    r['control_pass']=bool(r['core_mean_cross_PAE_A']<m['core_confidence']['max_mean_cross_PAE_A'] and pose['binder_CA_fold_rmsd_A']<=2 and pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and pose['binder_interface_CA_rmsd_after_target_fit_A']<=2 and pose['native_contact_recall']>=.5)
                else:r['control_pass']=r['core_mean_cross_PAE_A']>=m['core_confidence']['max_mean_cross_PAE_A']
                if c['condition']=='tagged':
                    tag=copy.deepcopy(s)
                    for res in list(tag[0]['B']):
                        if res.id[1]<=length:tag[0]['B'].detach_child(res.id)
                    pairs=sorted(contact_pairs(tag));g=geometry(tag,{i:i+189 for i in range(1,195)},full)
                    r['tag']={'target_contact_residue_pairs':pairs,'target_contact_pair_count':len(pairs),'geometry':g,
                        'mean_plddt':float(plddts[idx,194+length:].mean()),
                        'pass':len(pairs)==0 and all(g[k]==0 for k in g if 'clash_pairs_lt2A' in k)}
                else:r['tag']={'pass':True,'target_contact_pair_count':0,'note':'No tag in this condition'}
                group.append(r)
            chosen=select_core(group)
            for r in group:
                r['confidence_selected']=r['sample']==chosen;r['target_integrity_pass']=r['target_CA_rmsd_to_input_A']<=2
            rows.extend(group);print('Scored',seed,name,'selected',chosen,flush=True)
    assert len(rows)==140
    controls={};summaries=[]
    for condition in ('bare','tagged'):
        subset=[r for r in rows if r['confidence_selected'] and r['condition']==condition and r['kind']!='novel_design']
        controls[condition]={'pass':len(subset)==6 and all(r['control_pass'] for r in subset),'checks':subset}
    for c in m['cases']:
        if c['kind']!='novel_design':continue
        rs=[r for r in rows if r['id']==c['id'] and r['confidence_selected']]
        passed=controls[c['condition']]['pass'] and len(rs)==2 and all(r['screen']['pass'] for r in rs)
        summaries.append({'id':c['id'],'base_id':c['base_id'],'condition':c['condition'],'selected_models':rs,
            'core_pass_both_seeds':passed,'tag_independent_pass_both_seeds':passed and all(r['tag']['pass'] for r in rs),
            'passing_core_samples':sum(r['target_integrity_pass'] and r['screen']['pass'] for r in rows if r['id']==c['id']),
            'passing_core_and_tag_samples':sum(r['target_integrity_pass'] and r['screen']['pass'] and r['tag']['pass'] for r in rows if r['id']==c['id'])})
    result={'run':run.name,'manifest_sha256':sha(source/'manifest.json'),'controls':controls,'rows':rows,'candidate_summaries':summaries,
        'interpretation':'Core-only confidence/geometry avoids suffix dilution. Tag-contact-free selected poses are computational compatibility hypotheses, not experimental validation or proof of tag independence.'}
    (out/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'controls':{k:v['pass'] for k,v in controls.items()},'summary':[{k:v for k,v in c.items() if k!='selected_models'} for c in summaries]},indent=2))

if __name__=='__main__':main()
