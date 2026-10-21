"""Advisory matched-conformation scores; never assign mouse-binding labels."""
import argparse,json
from pathlib import Path
import numpy as np
from prepare_controls import ROOT,read_structure,sequence
from audit_control_samples import sample_structures
from analyze_generation_comparison import design_pose_metrics
from analyze_controls import pose_metrics
from target_integrity_selection import target_ca_rmsd,select_sample
from score_matched_controls import run_ipsae,interface_pae
from prepare_target import write_csv


def analyze(run):
    inputs=ROOT/'data/species_exploration';out=ROOT/'reports/species_exploration'
    m=json.loads((run/'inputs/manifest.json').read_text());seed=m['refolding']['seeds'][0];folder=run/f'seed-{seed}'
    pre=json.loads((folder/'preflight.json').read_text());assert len(pre)==len(m['cases']) and all(x['binder_template_excluded'] and x['msa_depth']==1 for x in pre)
    rows=[];summaries=[]
    for case in m['cases']:
        name=case['id'];ref=read_structure(inputs/f'{name}.pdb');case_rows=[]
        with np.load(folder/'fold_out_npz'/f'{name}.npz',allow_pickle=False) as d:
            paes=d['pae'].copy();plddts=d['plddt'].copy();n=194+len(case['binder_sequence'])
            assert paes.shape==(3,n,n) and plddts.shape==(3,n) and np.isfinite(paes).all() and np.isfinite(plddts).all()
            assert np.all(d['asym_id'][0,:194]==0) and np.all(d['asym_id'][0,194:]==1) and d['token_pad_mask'].all()
        for idx,s,metrics,upstream_selected in sample_structures(folder,name,restore_sidechains=True):
            assert sequence(s[0]['A'])==case['target_sequence'] and sequence(s[0]['B'])==case['binder_sequence']
            dest=ROOT/'outputs/species_scoring'/name/str(idx)
            scores=run_ipsae(s,paes[idx],plddts[idx],dest);scores.update(interface_pae(s,paes[idx]))
            pose=pose_metrics(s,ref) if case['kind']=='positive' else design_pose_metrics(s,ref)
            row={'id':name,'base_id':case['base_id'],'condition':case['condition'],'kind':case['kind'],'seed':seed,'sample':idx,
                 'target_CA_rmsd_to_input_A':target_ca_rmsd(s,ref),'metrics':metrics,'scores':scores,'pose':pose,
                 'confidence_selected':False,'upstream_confidence_selected':upstream_selected,'raw_scoring_output':str(dest.relative_to(ROOT))}
            if case['condition']!='maf' and case['arm']=='control':
                if case['kind']=='positive':
                    row['human_control_pass']=(metrics['design_to_target_iptm']>=.5 and pose['binder_CA_fold_rmsd_A']<=2 and pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and pose['binder_interface_CA_rmsd_after_target_fit_A']<=2 and pose['native_contact_recall']>=.5)
                else:row['human_control_pass']=metrics['design_to_target_iptm']<.5
            case_rows.append(row)
        selected=select_sample(case_rows,2)
        for row in case_rows:
            row['target_integrity_pass']=bool(row['target_CA_rmsd_to_input_A']<=2)
            row['confidence_selected']=row['sample']==selected
        selected_rows=[r for r in case_rows if r['confidence_selected']]
        summaries.append({'id':name,'base_id':case['base_id'],'condition':case['condition'],'known_role_on_human':case['kind'],'mouse_binding_status':'unknown','selected_model':selected_rows[0] if selected_rows else None})
        rows.extend(case_rows);print('Scored',name,flush=True)
    controls=[r for r in rows if r['confidence_selected'] and 'human_control_pass' in r]
    result={'run':run.name,'human_conformation_controls_pass':len(controls)==6 and all(r['human_control_pass'] for r in controls),'human_control_checks':controls,'candidate_summaries':summaries,'rows':rows,'plan':m['interpretation_rules'],'no_mouse_binding_classification':True,'limitations':['One seed and three samples per condition','No validated mouse positive control','Mouse reference binder pose is transferred from human; not experimental','AlphaFold monomer is not an experimentally validated bound/dimer/glycosylated target','No mouse affinity or cross-reactivity claim; low confidence cannot establish nonbinding']}
    (out/'results.json').write_text(json.dumps(result,indent=2))
    table=[]
    for r in summaries:
        s=r['selected_model']
        table.append({'base_id':r['base_id'],'condition':r['condition'],'selected_sample':s['sample'] if s else None,'interface_iptm':s['metrics']['design_to_target_iptm'] if s else None,'contact_pae_A':s['scores']['contact_pae_both_directions_A'] if s else None,'reference_interface_rmsd_A':s['pose']['binder_interface_CA_rmsd_after_target_fit_A'] if s else None,'reference_contact_recall':s['pose']['native_contact_recall'] if s else None,'mouse_binding_status':'unknown'})
    write_csv(out/'comparison.csv',table)
    print(json.dumps({'predictions':len(rows),'human_controls_pass':result['human_conformation_controls_pass']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();analyze(a.run.resolve())
