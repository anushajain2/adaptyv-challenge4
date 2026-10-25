"""Freeze a fresh bare/tagged confirmation of the unqualified MPNN reserve."""
import copy
import json
import re
import numpy as np
from novelty_search import ROOT,sha,save,now
from prepare_controls import read_structure,sequence,save_pdb,save_input_cif
from prepare_tag_comparison import append_placeholder
from prepare_assay_tag_review import SUFFIX
from prepare_tag_models import PROTOCOL
from repair_mutant_input_names import pinned_target_regex

RESERVE='apical_short_008_A41H_mpnn03'
ABLATION='apical_short_008_mpnn03_H41A'
BASES=['apical_short_008','apical_short_008_A41H',RESERVE,ABLATION,
       '2ds25_5_6wrw','scrambled_6wrw','caplacizumab_6wrw']
COHORT='All confidence-selected reserve confirmation bare/tagged models, including failures and controls.'


def main():
    src=ROOT/'data/local_redesign_panel';previous=json.loads((src/'manifest.json').read_text());lookup={c['id']:c for c in previous['cases']}
    oldtagpath=ROOT/'data/tag_comparison/manifest.json';oldtag=json.loads(oldtagpath.read_text())
    assert oldtag['suffix']==SUFFIX and len(SUFFIX)==64
    out=ROOT/'data/reserve_tag';out.mkdir(exist_ok=False);references=out/'references';references.mkdir();cases=[]
    bases={name:copy.deepcopy(lookup[name]) for name in BASES if name!=ABLATION}
    bases[ABLATION]={**copy.deepcopy(lookup[RESERVE]),'id':ABLATION,'histidine_positions':[],
        'parent_id':RESERVE,'local_redesign':False,'ablation_only':True,'mutations':['H41A']}
    for base in BASES:
        c=bases[base];reference=read_structure(src/f'{RESERVE if base==ABLATION else base}.pdb')
        if base==ABLATION:
            residue=reference[0]['B'][41];assert residue.resname=='HIS';residue.resname='ALA'
            for atom in list(residue):
                if atom.name not in {'N','CA','C','O'}:residue.detach_child(atom.id)
            c['binder_sequence']=sequence(reference[0]['B'])
            reserve=bases[RESERVE]['binder_sequence']
            assert [i+1 for i,(a,b) in enumerate(zip(reserve,c['binder_sequence'])) if a!=b]==[41]
        reference_path=references/f'{base}.pdb';save_pdb(reference,reference_path)
        for condition in ['bare','tagged']:
            name=base+'_'+condition;assert re.fullmatch(pinned_target_regex(),name)
            structure=append_placeholder(reference,SUFFIX) if condition=='tagged' else copy.deepcopy(reference)
            seq=sequence(structure[0]['B']);assert seq==c['binder_sequence']+(SUFFIX if condition=='tagged' else '')
            save_pdb(structure,out/f'{name}.pdb');save_input_cif(structure,out/f'{name}.cif',previous['target_sequence'],seq)
            np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(seq)))
            cases.append({**c,'id':name,'base_id':base,'condition':condition,'core_sequence':c['binder_sequence'],
                'core_length':len(c['binder_sequence']),'binder_sequence':seq,'reference':str(reference_path.relative_to(ROOT)),
                'reference_sha256':sha(reference_path),'histidine_positions':[41] if base in [RESERVE,'apical_short_008_A41H'] else []})
    protocol=copy.deepcopy(PROTOCOL);protocol['cohort']=COHORT
    plan={'created_utc':now(),'authorization':'Fresh user-authorized follow-up to an unqualified reserve; the prior automatic follow-up remains skipped.',
        'reserve':RESERVE,'matched_H41A_comparator':ABLATION,'base_ids':BASES,'seeds':[3509,3607],'samples_per_seed':5,
        'selection':'Identical target-integrity filter and minimum binder-core cross PAE to the prior tag experiment; do not choose poses by geometry, pKa or tag contacts.',
        'raw_screen':'Unchanged binder-core structural checks, prior calibrated cross-PAE threshold, contact PAE <=10 A; raw positive/negative controls must separate in each condition and both seeds.',
        'preparation':'Unchanged v2 tag-aware numerical preparation on every selected full construct. Strict target/core movement and internal geometry; full suffix motion also recorded. No selective retries or raw-failure rescue.',
        'primary_readout':'Descriptive held-out replication: report His41 free/bound pKas, bound-minus-free shift and His-only Kd6/Kd7.4 proxy continuously in both conditions and both seeds, with structural eligibility, matched baselines and controls. No new pass threshold replacing the old free-pKa band.',
        'legacy_result':'Also apply the unchanged free-pKa 6-7, shift <=-0.5, ratio >=2, target termini >=10 A and full prepared-control gates. Report that result explicitly; it is not overwritten by the continuous readout.',
        'sensitivity':'Separately record the already-declared interface movement sensitivity on core structures. No automatic promotion from a sensitivity-only result.',
        'tag_policy':'Full-construct chemistry for both free and bound calculations; only core His41 enters His-only ratio. Tag histidines are not counted as designed mechanism. Report all tag-target contacts, their PAE, reproducibility, and full-context clashes. Flexible uncertain contacts are flags, not experimental incompatibility.',
        'ablation_policy':'H41A on redesigned background tests modeled fold/pose preservation. A zero introduced-His score is tautological, not evidence that the ablation lacks experimental pH sensitivity; tag and other titratable residues remain.',
        'decision':'Archive outcomes without changing the frozen eight-design submission. No further paid job, new generation, alternative seeds or submission in this round. Recommend next action from the complete evidence; no automatic promotion.',
        'limitations':['Selected after observing a promising reserve; not an unbiased discovery screen','Same predictor, fresh seeds; no independent model validation',
            'Fixed bound conformation, no unbound ensemble or full proton linkage','Platform-default suffix, not confirmed challenge-specific construct','No experimental binding, pH or mouse claims'],
        'budget':{'round_max_usd':2,'project_gross_cap_usd':30}}
    m={k:copy.deepcopy(oldtag[k]) for k in ['boltzgen_commit','core_confidence','tag_screen','tag_provenance','placeholder_note']}
    m.update(scope='Fresh reserve confirmation and platform-default-tag sensitivity',prepared_utc=now(),cases=cases,
        target_sequence=previous['target_sequence'],arms=previous['arms'],suffix=SUFFIX,
        refolding={**previous['refolding'],'seeds':plan['seeds'],'diffusion_samples_per_seed':5,'export_pae':True,'target_msa':None},
        selection_rule=oldtag['selection_rule'],geometry_screen=oldtag['geometry_screen'],preparation=protocol,
        prospective_plan=plan,prior_tag_manifest_sha256=sha(oldtagpath),source_redesign_manifest_sha256=sha(src/'manifest.json'),
        max_prediction_seconds=2100,experiment_cap_usd=2,total_cap_usd=30)
    m['input_sha256']={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}
    save(out/'manifest.json',m);save(ROOT/'reports/reserve_tag/plan.json',m)
    print('Frozen',len(cases),'constructs;',len(cases)*10,'predictions;',len(SUFFIX),'residue suffix.')


if __name__=='__main__':main()
