"""Matched advisory predictions on bound human, AF human and AF mouse targets."""
import copy,hashlib,json,re
import numpy as np
from Bio import SeqIO
from Bio.PDB import Superimposer
from prepare_controls import ROOT,read_structure,sequence,residues,make_complex,save_pdb,save_input_cif
from datetime import datetime,timezone
from repair_mutant_input_names import pinned_target_regex


def main():
    out=ROOT/'data/species_exploration';out.mkdir(exist_ok=False)
    validation=json.loads((ROOT/'reports/expansion_validation/results.json').read_text());assert validation['controls_pass']
    selected=[c['id'] for c in validation['candidate_summaries'] if c['all_selected_seeds_full_screen']]
    assert len(selected)<=6
    original=json.loads((ROOT/'data/expansion_validation/manifest.json').read_text())
    bases=[(c,ROOT/'data/expansion_validation') for c in original['cases'] if c['id'] in selected or c['arm']=='control']
    hist=json.loads((ROOT/'data/his_confirmation/manifest.json').read_text())
    bases += [(c,ROOT/'data/his_confirmation') for c in hist['cases'] if c['id'] in {'apical_short_008','apical_short_008_A41H'}]
    human=read_structure(ROOT/'data/apical_expansion_fixed/target_6wrw.cif')
    afhuman=read_structure(ROOT/'data/mouse_context_audit/human_P02786.cif')
    assert sequence(afhuman[0]['A'])==str(next(SeqIO.parse(ROOT/'data/raw/P02786.fasta','fasta')).seq)
    audit=json.loads((ROOT/'reports/mouse_context_audit/results.json').read_text())
    anchors=[h for h,m in audit['anchors']]
    fit=Superimposer();fit.set_atoms([human[0]['A'][h-189]['CA'] for h in anchors],[afhuman[0]['A'][h]['CA'] for h in anchors]);fit.apply(list(afhuman.get_atoms()))
    afmouse=read_structure(ROOT/'reports/mouse_context_audit/mouse_model_in_human_frame.pdb')
    targets={'hbound':residues(human[0]['A']),'haf':[afhuman[0]['A'][i] for i in range(190,384)],'maf':[afmouse[0]['A'][i] for i in range(192,386)]}
    cases=[]
    for condition,target in targets.items():
        for base,src in bases:
            binder=read_structure(src/f"{base['id']}.pdb")[0]['B'];s=make_complex(target,residues(binder))
            name=condition+'_'+base['id'];assert re.fullmatch(pinned_target_regex(),name)
            t_seq=sequence(s[0]['A']);b_seq=sequence(s[0]['B']);assert len(t_seq)==194 and b_seq==base['binder_sequence']
            save_pdb(s,out/f'{name}.pdb');save_input_cif(s,out/f'{name}.cif',t_seq,b_seq)
            np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(b_seq)))
            cases.append({**base,'id':name,'base_id':base['id'],'condition':condition,'target_sequence':t_seq,'target_native_map':{i:i+(191 if condition=='maf' else 189) for i in range(1,195)},'reference_note':'Binder coordinates excluded from prediction templates; transferred human reference pose used only for descriptive consistency metrics. Mouse binding status unknown.'})
    m={**hist,'cases':cases,'selected_ids':selected+['apical_short_008','apical_short_008_A41H'],
       'scope':'Exploratory species/conformation comparison, not calibrated mouse binding or affinity prediction',
       'refolding':{**hist['refolding'],'seeds':[1601],'diffusion_samples_per_seed':3},
       'max_prediction_seconds':1500,'experiment_cap_usd':1.5,
       'conditions':{'hbound':'Bound human 6WRW native 190-383','haf':'AlphaFold human P02786 native 190-383, same conserved-core alignment','maf':'AlphaFold mouse Q62351 native 192-385, same conserved-core alignment'},
       'human_AF_anchor_rmsd_A':float(fit.rms),'mouse_AF_anchor_rmsd_A':audit['anchor_CA_rmsd_A'],
       'interpretation_rules':{'human_controls':'Record whether positive interface recovery and presumed-negative confidence behavior persist on both human conformations. Failure makes mouse interpretation additionally confounded.',
         'mouse_controls':'No known mouse-binding positive control is available in this panel; human control labels do NOT transfer to mouse.',
         'selection':'Same target-integrity filter and confidence selection; no alternate seeds. Report all three samples.',
         'no_binary_crossreactivity_call':'Neither high nor low predicted scores establish or rule out mouse binding. AF monomer targets have no experimentally validated bound conformation or full glycan/dimer screen here.'},
       'source_results_sha256':hashlib.sha256((ROOT/'reports/expansion_validation/results.json').read_bytes()).hexdigest()}
    # A single global target sequence/map would be false for this mixed-target panel.
    m.pop('target_sequence',None);m.pop('target_native_map',None)
    for key in ('screen','arms','controls','prospective_decision','selection_details','structural_preparation','next_round_allocation','selection_result_sha256'):m.pop(key,None)
    m['source_run']=validation['run'];m['preparation_scope']='Raw prediction comparison only; no mouse geometry qualification, pKa, or post-minimization rescue is performed.'
    m['prepared_utc']=datetime.now(timezone.utc).isoformat()
    m['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(m,indent=2))
    report=ROOT/'reports/species_exploration';report.mkdir(exist_ok=False);(report/'plan.json').write_text(json.dumps(m,indent=2))
    print(json.dumps({'cases':len(cases),'predictions':len(cases)*3,'human_AF_anchor_rmsd_A':float(fit.rms)},indent=2))
if __name__=='__main__':main()
