"""Diagnose receptor deformation without changing the frozen validation decision."""
import json
import numpy as np
from Bio.PDB import Superimposer
from analyze_controls import contact_pairs
from prepare_controls import ROOT, read_structure


def main():
    data = json.loads((ROOT / 'reports/broad_validation/results.json').read_text())
    ref = read_structure(ROOT / 'data/broad_validation/2ds25_5_6wrw.pdb')
    native = contact_pairs(ref)
    anchors = sorted({a for a, _ in native})
    binder = sorted({b for _, b in native})
    rows = []
    for row in data['rows']:
        if row['kind'] != 'positive':
            continue
        pred = read_structure(ROOT / row['raw_scoring_output'] / 'model.pdb')
        sup = Superimposer()
        ids = [r.id[1] for r in ref[0]['A']]
        sup.set_atoms([ref[0]['A'][i]['CA'] for i in ids], [pred[0]['A'][i]['CA'] for i in ids])
        rot, trans = sup.rotran
        errors = [{'crop_position': i, 'native_position': i + 189,
                   'CA_displacement_after_global_target_fit_A': float(np.linalg.norm(
                       pred[0]['A'][i]['CA'].coord @ rot + trans - ref[0]['A'][i]['CA'].coord))}
                  for i in ids]
        sup.set_atoms([ref[0]['A'][i]['CA'] for i in anchors], [pred[0]['A'][i]['CA'] for i in anchors])
        rot, trans = sup.rotran
        interface_rms = float(np.sqrt(np.mean([np.sum((pred[0]['B'][i]['CA'].coord @ rot + trans - ref[0]['B'][i]['CA'].coord)**2) for i in binder])))
        rows.append({'seed': row['seed'], 'sample': row['sample'],
                     'confidence_selected': row['confidence_selected'], 'frozen_control_pass': row['control_pass'],
                     'original_pose': row['pose'],
                     'diagnostic_target_interface_anchor_fit_rmsd_A': float(sup.rms),
                     'diagnostic_binder_interface_rmsd_after_target_interface_anchor_fit_A': interface_rms,
                     'target_residues_displaced_over_4A': sum(r['CA_displacement_after_global_target_fit_A'] > 4 for r in errors),
                     'target_displacements': errors})
    result = {'rows': rows, 'reference_target_interface_anchor_positions_crop': anchors,
              'control_samples_passing_original_rule': sum(r['frozen_control_pass'] for r in rows),
              'interpretation': 'Post hoc diagnostic only. Fit uses reference target-interface residues, never binder atoms. The original global-target-frame gate and batch result remain unchanged.'}
    (ROOT / 'reports/broad_validation/control_deformation_audit.json').write_text(json.dumps(result, indent=2))
    for row in rows:
        print(row['seed'], row['sample'], 'selected', row['confidence_selected'], 'pass', row['frozen_control_pass'],
              'global target RMSD', round(row['original_pose']['target_CA_rmsd_to_reference_A'], 3),
              'local interface RMSD', round(row['diagnostic_binder_interface_rmsd_after_target_interface_anchor_fit_A'], 3),
              'target residues >4A', row['target_residues_displaced_over_4A'])


if __name__ == '__main__':
    main()
