"""Locate historical movement failures without changing any historical decision."""
import json
import numpy as np
from novelty_search import ROOT, sha, save
from prepare_controls import read_structure

BACKBONE = {'N', 'CA', 'C', 'O'}


def motion_audit(raw, prepared):
    rows = []
    for chain, other in [('A', 'B'), ('B', 'A')]:
        other_xyz = np.array([a.coord for a in raw[0][other].get_atoms() if a.element != 'H'])
        for residue in raw[0][chain]:
            heavy = [a for a in residue if a.element != 'H']
            distance = float(np.linalg.norm(np.array([a.coord for a in heavy])[:, None] - other_xyz, axis=2).min())
            for atom in heavy:
                if atom.name not in BACKBONE:
                    continue
                displacement = float(np.linalg.norm(atom.coord - prepared[0][chain][residue.id][atom.name].coord))
                rows.append({'chain': chain, 'position': residue.id[1], 'residue': residue.resname,
                             'atom': atom.name, 'displacement_A': displacement,
                             'atom_distance_to_other_chain_A': float(np.linalg.norm(atom.coord - other_xyz, axis=1).min()),
                             'residue_distance_to_other_chain_A': distance, 'interface_residue_6A': distance <= 6})
    chains = {}
    for chain in ['A', 'B']:
        subset = [r for r in rows if r['chain'] == chain]
        interface = [r for r in subset if r['interface_residue_6A']]
        chains[chain] = {'backbone_rms_A': float(np.sqrt(np.mean([r['displacement_A']**2 for r in subset]))),
                         'CA_max_A': max(r['displacement_A'] for r in subset if r['atom'] == 'CA'),
                         'interface_backbone_max_A': max((r['displacement_A'] for r in interface), default=None)}
    # Prospective sensitivity criterion only. No reclassification of old records.
    sensitive = all(v['backbone_rms_A'] <= .5 and v['CA_max_A'] <= 1 and
                    v['interface_backbone_max_A'] is not None and v['interface_backbone_max_A'] <= 1
                    for v in chains.values())
    return {'chains': chains, 'largest_movements': sorted(rows, key=lambda r: -r['displacement_A'])[:12],
            'all_backbone_over_1A': [r for r in rows if r['displacement_A'] > 1],
            'interface_sensitivity_pass': sensitive}


def main():
    source = ROOT/'reports/final_his_panel/preparation/results.json'
    result = json.loads(source.read_text()); rows = []
    for row in result['rows']:
        if row['status'] != 'ok':
            continue
        raw = ROOT/row['source']; prepared = ROOT/row['dest']/'prepared.pdb'
        assert sha(raw) == row['source_sha256'] and sha(prepared) == row['prepared_sha256']
        rows.append({'id': row['id'], 'seed': row['seed'], 'historical_movement_pass': row['displacement']['pass'],
                     'raw': str(raw.relative_to(ROOT)), 'raw_sha256': sha(raw),
                     'prepared': str(prepared.relative_to(ROOT)), 'prepared_sha256': sha(prepared),
                     **motion_audit(read_structure(raw), read_structure(prepared))})
    save(ROOT/'reports/local_redesign/movement_audit.json', {
        'source_sha256': sha(source), 'rows': rows,
        'interpretation': 'Retrospective diagnosis, not a rescue of historical failures. A global maximum includes remote carbonyl oxygens. The separate interface sensitivity criterion is exploratory and not experimentally calibrated.'})
    print(json.dumps([{k: r[k] for k in ['id', 'seed', 'all_backbone_over_1A']} for r in rows if not r['historical_movement_pass']], indent=2))


if __name__ == '__main__':
    main()
