"""Check the direction report against saved structures, pKa tables and warnings."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from analyze_controls import pose_metrics
from direction_check import HIS_ATOMS, histidines, read_structure
from prepare_controls import sequence
from prepare_target import ROOT


def main():
    report = ROOT / 'reports/confirmation_direction'
    result = json.loads((report / 'results.json').read_text())
    validation_path = ROOT / 'reports/confirmation/results.json'
    assert hashlib.sha256(validation_path.read_bytes()).hexdigest() == result['plan']['validation_results_sha256']
    validation = json.loads(validation_path.read_text())
    parents = {p['id']: p for p in validation['candidate_summaries']}
    audited = []
    for pair in result['pairs']:
        parent = parents[pair['parent']]
        assert parent['provisional_computational_parent'] and validation['controls_pass']
        for row in pair['models']:
            assert row['status'] == 'ok'
            dest = ROOT / row['evidence_directory']
            bound = read_structure(dest / 'complex.pdb')
            free = read_structure(dest / 'free_binder.pdb')
            assert sequence(bound[0]['B']) == sequence(free[0]['B']) == pair['sequence']
            assert {c.id for c in bound[0]} == {'A', 'B'} and {c.id for c in free[0]} == {'B'}
            positions = list(map(int, row['positions'].split(',')))
            for pos in positions:
                assert bound[0]['B'][pos].resname == 'HIS' and HIS_ATOMS <= set(bound[0]['B'][pos].child_dict)
            source_model = next(m for m in parent['selected_models'] if m['seed'] == row['prediction_seed'])
            source = ROOT / source_model['raw_scoring_output'] / 'model.pdb'
            assert hashlib.sha256(source.read_bytes()).hexdigest() == row['source_model_sha256']
            reference = read_structure(source)
            assert sequence(bound[0]['A']) == sequence(reference[0]['A'])
            values = {}
            warnings = []
            for name in ('complex', 'free_binder'):
                values[name] = histidines(json.loads((dest / f'{name}.pka.json').read_text()), 'B', positions)
                warnings.extend(json.loads((dest / f'{name}.propka_warnings.json').read_text()))
            assert not any('Missing atoms or failed protonation' in w for w in warnings)
            for i in range(2):
                assert math.isclose(values['free_binder'][i], row[f'free_pka_{i+1}'], abs_tol=1e-10)
                assert math.isclose(values['complex'][i], row[f'bound_pka_{i+1}'], abs_tol=1e-10)
                assert math.isclose(values['complex'][i]-values['free_binder'][i], row[f'shift_{i+1}'], abs_tol=1e-10)
            full = read_structure(dest / 'full_complex.pdb')
            xyz = np.array([a.coord for c in full[0] if c.id != 'C' for a in c.get_atoms() if a.element != 'H'])
            binder = np.array([a.coord for a in full[0]['C'].get_atoms() if a.element != 'H'])
            clashes = int(cKDTree(xyz).query_ball_point(binder, np.nextafter(2., 0.), return_length=True).sum())
            assert clashes == row['full_dimer_clashes']
            relaxation = json.loads((dest / 'complex.relaxation.json').read_text())
            assert math.isfinite(relaxation['final_energy_kJ_mol'])
            audited.append({'variant': pair['id'], 'prediction_seed': row['prediction_seed'],
                            'rebuild_seed': row['seed'], 'matches_saved_pka_tables': True,
                            'complete_engineered_histidines': True, 'propka_warnings': warnings,
                            'recomputed_full_context_clashes': clashes,
                            'relaxed_pose_against_parent': pose_metrics(bound, reference),
                            'evidence_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir() if p.is_file()}})
    assert len(audited) == 8
    (report / 'evidence_audit.json').write_text(json.dumps({'models_checked': len(audited), 'all_pass': True, 'rows': audited}, indent=2))
    print(json.dumps({'models_checked': len(audited), 'all_pass': True,
                      'max_binder_backbone_pose_rmsd_after_relaxation_A': max(r['relaxed_pose_against_parent']['binder_CA_pose_rmsd_after_target_fit_A'] for r in audited)}, indent=2))


if __name__ == '__main__':
    main()
