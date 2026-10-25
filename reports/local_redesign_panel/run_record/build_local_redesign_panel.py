"""Construct the preregistered sequence-only refolding comparison."""
import copy
import json
import numpy as np
from Bio.SeqUtils import seq3
from novelty_search import ROOT, sha, save, now
from prepare_controls import read_structure, sequence, save_pdb, save_input_cif
from receptor_preparation_v2 import PROTOCOL

COHORT = 'Every confidence-selected local redesign candidate, baseline and control, including failures.'


def main():
    plan_path = ROOT/'data/local_redesign/plan.json'; plan = json.loads(plan_path.read_text())
    sampled_path = ROOT/'reports/local_redesign/mpnn/results.json'; sampled = json.loads(sampled_path.read_text())
    assert sampled['plan_sha256'] == sha(plan_path)
    out = ROOT/'data/local_redesign_panel';out.mkdir(exist_ok=False)
    old = json.loads((ROOT/'data/final_his_panel/manifest.json').read_text())
    m = {k:copy.deepcopy(old[k]) for k in ['target_sequence','screen','selection_rule','refolding','arms']}
    m['arms'] += [a for a in json.loads((ROOT/'data/his_confirmation/manifest.json').read_text())['arms'] if a['id']=='apical_short']
    cases = []

    def copy_case(c, src):
        cases.append(copy.deepcopy(c))
        for ext in ['.pdb','.cif','.npz']:
            (out/(c['id']+ext)).write_bytes((src/(c['id']+ext)).read_bytes())

    for c in old['cases']:
        if c['arm']=='control':copy_case(c,ROOT/'data/final_his_panel')
    for family in plan['families']:
        src = ROOT/family['baseline_input_directory']; base = json.loads((src/'manifest.json').read_text())
        selected = next(f for f in sampled['families'] if f['id']==family['id'])['selected']
        baseline = next(c for c in base['cases'] if c['id']==family['id'])
        assert baseline['binder_sequence'] == family['sequence']
        for name in [family['id'],family['unmodified_parent']]:
            c = copy.deepcopy(next(c for c in base['cases'] if c['id']==name)); c['comparison_only']=True
            c['histidine_positions']=[family['histidine_position']] if name==family['id'] else []
            copy_case(c,src)
        for row in selected:
            structure = read_structure(src/f"{family['id']}.pdb")
            for position,(aa,new) in enumerate(zip(family['sequence'],row['sequence']),1):
                if aa == new:continue
                assert position in family['mutable_positions']
                residue = structure[0]['B'][position];residue.resname=seq3(new).upper()
                for atom in list(residue):
                    if atom.name not in {'N','CA','C','O'}:residue.detach_child(atom.id)
            assert sequence(structure[0]['B']) == row['sequence']
            name=row['id'];save_pdb(structure,out/f'{name}.pdb');save_input_cif(structure,out/f'{name}.cif',m['target_sequence'],row['sequence'])
            np.savez_compressed(out/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(row['sequence'])))
            cases.append({**baseline,'id':name,'binder_sequence':row['sequence'],'parent_id':family['id'],
                          'unmodified_parent_id':family['unmodified_parent'],'histidine_positions':[family['histidine_position']],
                          'mutations':row['mutations'],'local_redesign':True,'mpnn_ensemble_nll':row['ensemble_nll']})
    assert 7 < len(cases) <= 19 and len({c['binder_sequence'] for c in cases})==len(cases)
    protocol=copy.deepcopy(PROTOCOL);protocol['cohort']=COHORT
    m.update(cases=cases,prepared_utc=now(),source_run='local-proteinmpnn-two-backbone-ranking',
             scope='Local sequence redesign around fixed A41H and S90H; matched His and original parents with controls',
             refolding={**m['refolding'],'seeds':plan['refolding']['seeds'],'diffusion_samples_per_seed':5},
             structural_preparation=protocol,prospective_plan=plan,source_mpnn_sha256=sha(sampled_path),
             max_prediction_seconds=1500,experiment_cap_usd=2,total_cap_usd=30)
    m['screen']['replicates']='Both confidence-selected fresh seeds, five samples each; no alternative-sample rescue'
    m['screen']['decision']='Evaluate frozen local redesign panel; no automatic generation expansion'
    m['input_sha256']={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    save(out/'manifest.json',m);save(ROOT/'reports/local_redesign_panel/plan.json',m)
    print(json.dumps({'cases':len(cases),'new_variants':sum(bool(c.get('local_redesign')) for c in cases),'predictions':10*len(cases)},indent=2))


if __name__ == '__main__':main()
