"""Prospective matched core/tag experiment; no tag coordinates used as templates."""
import copy,json,re
from datetime import datetime,timezone
import numpy as np
from Bio.PDB import Atom,Residue
from Bio.SeqUtils import seq3
from prepare_controls import ROOT,read_structure,sequence,save_input_cif,save_pdb
from prepare_assay_tag_review import SUFFIX
from receptor_preparation_v2 import sha,PROTOCOL
from repair_mutant_input_names import pinned_target_regex

BASE_IDS=['apical_short_008','apical_short_008_A41H','expansion_apical_short_026','expansion_apical_short_026_E19H','2ds25_5_6wrw','scrambled_6wrw','caplacizumab_6wrw']
SEEDS=[2003,2011]

def append_placeholder(s,suffix):
    result=copy.deepcopy(s);chain=result[0]['B'];n=len(chain)
    origin=np.asarray(chain[n]['C'].coord)+[10.,0,0]
    for i,aa in enumerate(suffix,1):
        r=Residue.Residue((' ',n+i,' '),seq3(aa).upper(),' ')
        center=origin+[3.8*i,0.,0.]
        for name,offset in [('N',[-1.2,0,0]),('CA',[0,0,0]),('C',[1.2,.6,0]),('O',[1.5,1.7,0])]:
            r.add(Atom.Atom(name,np.array(center+offset,dtype=float),0.,1.,' ',f'{name:>4}',0,element=name[0]))
        chain.add(r)
    return result

def main():
    src=ROOT/'data/expansion_his_confirmation';previous=json.loads((src/'manifest.json').read_text())
    calibration_path=ROOT/'reports/expansion_his_confirmation/results.json';calibration=json.loads(calibration_path.read_text())
    controls=[r for r in calibration['rows'] if r['confidence_selected'] and r['arm']=='control']
    positive_max=max(r['metrics']['interaction_pae'] for r in controls if r['kind']=='positive')
    negative_min=min(r['metrics']['interaction_pae'] for r in controls if r['kind']=='negative_hypothesis')
    assert positive_max<negative_min
    threshold=(positive_max+negative_min)/2
    out=ROOT/'data/tag_comparison';out.mkdir(exist_ok=False);cases=[]
    lookup={c['id']:c for c in previous['cases']}
    for base in BASE_IDS:
        old=lookup[base];reference=read_structure(src/f'{base}.pdb')
        for condition in ('bare','tagged'):
            name=base+'_'+condition;assert re.fullmatch(pinned_target_regex(),name)
            s=append_placeholder(reference,SUFFIX) if condition=='tagged' else copy.deepcopy(reference)
            seq=sequence(s[0]['B']);assert seq==old['binder_sequence']+(SUFFIX if condition=='tagged' else '')
            save_pdb(s,out/f'{name}.pdb');save_input_cif(s,out/f'{name}.cif',previous['target_sequence'],seq)
            np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(seq)))
            cases.append({**old,'id':name,'base_id':base,'condition':condition,'core_length':len(old['binder_sequence']),
                          'core_sequence':old['binder_sequence'],'binder_sequence':seq,'target_sequence':previous['target_sequence'],
                          'reference':f'data/expansion_his_confirmation/{base}.pdb'})
    manifest={'scope':'Paired bare-versus-platform-default-tag sensitivity; no submission or new generation',
        'prepared_utc':datetime.now(timezone.utc).isoformat(),'boltzgen_commit':previous['boltzgen_commit'],
        'target_sequence':previous['target_sequence'],'cases':cases,'arms':previous['arms'],
        'refolding':{**previous['refolding'],'seeds':SEEDS,'diffusion_samples_per_seed':5,'export_pae':True,'target_msa':None},
        'selection_rule':{'max_target_CA_rmsd_to_input_A':2.,'ranking':'Lowest symmetric mean PAE across ALL 194 target residues and binder CORE only; tie by sample index. No geometry/contact/shift-based sample selection.'},
        'core_confidence':{'max_mean_cross_PAE_A':threshold,'calibration_source':str(calibration_path.relative_to(ROOT)),
            'calibration_sha256':sha(calibration_path),'previous_selected_positive_max_A':positive_max,'previous_selected_negative_min_A':negative_min,
            'policy':'Midpoint of prior selected positive max and presumed-negative min, fixed before prediction. This is an exploratory replacement for length-confounded whole-chain ipTM, not a validated affinity gate.',
            'max_contact_PAE_A':10.},
        'geometry_screen':previous['screen'],
        'tag_screen':{'max_target_contact_residue_pairs_lt5A':0,'max_full_context_clash_pairs_lt2A':0,
            'interpretation':'Conservative tag-interference flag, not proof that a flexible tag cannot touch the target. Tags never contribute to core hotspot or confidence scores.'},
        'control_policy':'In EACH condition/seed: positive must pass core PAE plus native pose/fold/contact recall; presumed negatives must fail the core PAE threshold. No threshold tuning after results.',
        'preparation':PROTOCOL,
        'direction_policy':'Apply v2 preparation once to ALL selected models, including controls, in both conditions. Only variants with controls and core raw screens passing both seeds enter direction analysis; require prepared core movement/geometry and warning-free full free/bound PROPKA. Tag contacts block tag-independent qualification; report effects as sensitivity. Do not use core-only artificially terminated pKas for tagged constructs.',
        'tag_provenance':json.loads((ROOT/'data/assay_tag_review/manifest.json').read_text())['evidence_status'],
        'suffix':SUFFIX,'placeholder_note':'Tag backbone coordinates merely ensure parser/export completeness. Entire binder plus tag is excluded from structure templates; no tag conformation is supplied to the model.',
        'max_prediction_seconds':2100,'experiment_cap_usd':2,'total_cap_usd':30}
    manifest['input_sha256']={p.name:sha(p) for p in out.iterdir()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    report=ROOT/'reports/tag_comparison';report.mkdir(exist_ok=False)
    (report/'plan.json').write_text(json.dumps(manifest,indent=2))
    print('Prepared',len(cases)*len(SEEDS)*5,'predictions; frozen core PAE threshold',threshold)

if __name__=='__main__':main()
