"""Prepare a sequence-only mutant refold with its parent and matched controls."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from prepare_controls import read_structure, save_input_cif, save_pdb, sequence

ROOT = Path(__file__).resolve().parents[1]
PARENT = 'apical_long_015'
MUTANT = 'apical_long_015__D51H_T70H'


def mutate_reference(structure):
    import copy
    s = copy.deepcopy(structure)
    for position, original in ((51, 'ASP'), (70, 'THR')):
        residue = s[0]['B'][position]
        assert residue.resname == original
        residue.resname = 'HIS'
        for atom in list(residue):
            if atom.name not in {'N', 'CA', 'C', 'O'}:
                residue.detach_child(atom.id)
    return s


def main():
    validation = ROOT / 'reports/confirmation/results.json'
    direction = ROOT / 'reports/confirmation_direction/results.json'
    v, d = json.loads(validation.read_text()), json.loads(direction.read_text())
    assert v['controls_pass']
    parent = next(p for p in v['candidate_summaries'] if p['id'] == PARENT)
    hit = next(p for p in d['pairs'] if p['id'] == MUTANT)
    assert parent['provisional_computational_parent'] and hit['direction_criteria_pass_all_four']
    source = ROOT / 'data/confirmation'
    m = json.loads((source / 'manifest.json').read_text())
    cases = [c for c in m['cases'] if c['id'] == PARENT or c['arm'] == 'control']
    out = ROOT / 'data/mutant_validation'
    out.mkdir(exist_ok=False)
    for case in cases:
        for ext in ('.pdb', '.cif', '.npz'):
            name = case['id'] + ext
            assert hashlib.sha256((source / name).read_bytes()).hexdigest() == m['input_sha256'][name]
            (out / name).write_bytes((source / name).read_bytes())
    parent_case = next(c for c in cases if c['id'] == PARENT)
    mutant = mutate_reference(read_structure(out / f'{PARENT}.pdb'))
    seq = sequence(mutant[0]['B'])
    assert seq == hit['sequence'] and len(seq) == len(parent['sequence']) == 98
    assert [i+1 for i, (a, b) in enumerate(zip(seq, parent['sequence'])) if a != b] == [51, 70]
    save_pdb(mutant, out / f'{MUTANT}.pdb')
    save_input_cif(mutant, out / f'{MUTANT}.cif', m['target_sequence'], seq)
    np.savez_compressed(out / f'{MUTANT}.npz', design_mask=np.array([False]*194 + [True]*98))
    cases.insert(1, {**parent_case, 'id': MUTANT, 'binder_sequence': seq,
                     'parent_id': PARENT, 'mutations': ['D51H', 'T70H'],
                     'reference_note': 'Same original generated parent backbone, with only sequence labels changed. Mutant binder coordinates are excluded from prediction templates; reference is used only for scoring.',
                     'generated_source': parent_case['generated_source']})
    m.update(scope='Independent mutant-sequence refolding; no new generation or mutation search',
             prepared_utc=datetime.now(timezone.utc).isoformat(), source_run=v['run'], cases=cases,
             selected_ids=[PARENT, MUTANT],
             selection_result_sha256=hashlib.sha256(validation.read_bytes()).hexdigest(),
             direction_result_sha256=hashlib.sha256(direction.read_bytes()).hexdigest(),
             refolding={**m['refolding'], 'seeds': [809, 907], 'diffusion_samples_per_seed': 5},
             max_prediction_seconds=900, experiment_cap_usd=1.5,
             mutation_decision='Controls, unmodified parent and mutant must all pass in both selected seeds. Retain failures without substitutions or threshold changes. A pass supports structural consistency only, not mutant affinity or pH switching.')
    m['input_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}
    (out / 'manifest.json').write_text(json.dumps(m, indent=2))
    report = ROOT / 'reports/mutant_validation'
    report.mkdir(exist_ok=False)
    (report / 'plan.json').write_text(json.dumps(m, indent=2))
    print(json.dumps({'cases': len(cases), 'predictions': len(cases)*10, 'seeds': [809, 907], 'mutant': MUTANT}, indent=2))


if __name__ == '__main__':
    main()
