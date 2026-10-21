"""Select on target conformation and confidence, never on desired binder geometry."""
import math
import numpy as np
from Bio.PDB import Superimposer


def target_ca_rmsd(predicted, reference):
    ids = [r.id for r in reference[0]['A']]
    assert len(ids) == 194 and [r.id for r in predicted[0]['A']] == ids
    ref_atoms = [reference[0]['A'][i]['CA'] for i in ids]
    pred_atoms = [predicted[0]['A'][i]['CA'] for i in ids]
    assert all(predicted[0]['A'][i].resname == reference[0]['A'][i].resname for i in ids)
    if not np.isfinite([a.coord for a in pred_atoms + ref_atoms]).all():
        return float('inf')
    fit = Superimposer()
    fit.set_atoms(ref_atoms, pred_atoms)
    return float(fit.rms)


def select_sample(rows, max_target_rmsd=2.0):
    eligible = [r for r in rows if math.isfinite(r['target_CA_rmsd_to_input_A'])
                and r['target_CA_rmsd_to_input_A'] <= max_target_rmsd
                and all(math.isfinite(r['metrics'][k]) for k in ('iptm', 'ptm'))]
    if not eligible:
        return None
    return max(eligible, key=lambda r: (.8*r['metrics']['iptm'] + .2*r['metrics']['ptm'], -r['sample']))['sample']
