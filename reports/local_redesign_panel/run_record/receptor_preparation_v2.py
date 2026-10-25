"""Uniform sensitivity analysis; historical raw decisions are immutable.

Unlike v1, receptor side chains can relax. Both chains' backbones are restrained.
No model-specific retries, residue repairs, alternative rotamers or score selection.
"""
import argparse
import hashlib
import json
import math
import multiprocessing
import random
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from prepare_controls import ROOT, read_structure, sequence, geometry
from prepare_target import save_pdb, ProteinSelect
from direction_check import pka_table

BACKBONE = {'N', 'CA', 'C', 'O'}
PROTOCOL = {
    'version': 'receptor-preparation-v2',
    'interpretation': 'Uniform post-hoc sensitivity analysis, specified before revised outcomes; never replaces raw or v1 decisions.',
    'method': 'Amber14/GBn2, 12 A cutoff, both backbones restrained; all side chains and hydrogens mobile',
    'restraint_expression': 'k*((x-x0)^2+(y-y0)^2+(z-z0)^2)',
    'backbone_k_kJ_mol_nm2': 1000., 'seed': 43, 'hydrogen_placement_pH': 7.4,
    'max_iterations': 500, 'minimization_timeout_seconds': 120, 'hard_timeout_seconds': 180,
    'backbone_rms_limit_A': .5, 'backbone_max_limit_A': 1.,
    'internal_nonbonded_distance_cutoff_A': 1.8,
    'internal_nonbonded_pairs_allowed': 0,
    'protonation': 'Strip explicit H; PROPKA assigns H consistently to free and bound forms of the same prepared conformation. Any incomplete/protonation-failed group invalidates that calculation.',
    'cohort': 'All 28 confidence-selected models from expansion_his_confirmation, irrespective of raw outcome.',
    'policy': 'One preparation per source; no selective retries; source sequence and original atom identities preserved; raw screen failures cannot be rescued. Original pose/contact/context gates and pKa heuristics unchanged.',
    'limitations': ['Fixed conformation free/bound pKa comparison', 'No unbound ensemble or full proton linkage', 'No assay tag modeled', 'No experimental binding claim', 'Finite minimization is not proof of convergence or global stereochemical validity'],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def displacement(raw, prepared):
    result = {}
    for chain in ('A', 'B'):
        assert sequence(raw[0][chain]) == sequence(prepared[0][chain])
        groups = {'backbone': [], 'sidechain': []}
        for r in raw[0][chain]:
            for a in r:
                if a.element == 'H':
                    continue
                assert a.name in prepared[0][chain][r.id], f'Missing source atom {chain}:{r.id}:{a.name}'
                groups['backbone' if a.name in BACKBONE else 'sidechain'].append(
                    float(np.linalg.norm(a.coord - prepared[0][chain][r.id][a.name].coord)))
        result[chain] = {k: {'rms_A': float(np.sqrt(np.mean(np.square(v)))), 'max_A': max(v)} for k, v in groups.items()}
    result['pass'] = all(result[c]['backbone']['rms_A'] <= PROTOCOL['backbone_rms_limit_A'] and
                         result[c]['backbone']['max_A'] <= PROTOCOL['backbone_max_limit_A'] for c in ('A', 'B'))
    return result


def internal_geometry(path):
    """Severe-overlap screen, not a full stereochemical validation."""
    from openmm import app, unit
    pdb = app.PDBFile(str(path))
    xyz = np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    bonds = {tuple(sorted((a.index, b.index))) for a, b in pdb.topology.bonds()}
    records = []
    for chain in pdb.topology.chains():
        atoms = [a for a in chain.atoms() if a.element.symbol != 'H']
        for i, j in sorted(cKDTree(xyz[[a.index for a in atoms]]).query_pairs(np.nextafter(1.8, 0.))):
            a, b = atoms[i], atoms[j]
            if tuple(sorted((a.index, b.index))) in bonds:
                continue
            records.append({'a': f'{chain.id}:{a.residue.id}:{a.residue.name}:{a.name}',
                            'b': f'{chain.id}:{b.residue.id}:{b.residue.name}:{b.name}',
                            'distance_A': float(np.linalg.norm(xyz[a.index] - xyz[b.index]))})
    return {'nonbonded_pairs_lt1_8A': records, 'pass': len(records) == 0}


def relax(source, out):
    from pdbfixer import PDBFixer
    import openmm as mm
    from openmm import app, unit
    random.seed(PROTOCOL['seed']); np.random.seed(PROTOCOL['seed'])
    fixer = PDBFixer(filename=str(source))
    fixer.findMissingResidues(); fixer.missingResidues = {}
    fixer.findMissingAtoms(); fixer.addMissingAtoms(seed=PROTOCOL['seed'])
    xyz = np.array(fixer.positions.value_in_unit(unit.nanometer))
    for r in fixer.topology.residues():
        atoms = {a.name: a.index for a in r.atoms()}
        if {'CA', 'C', 'O', 'OXT'} <= set(atoms):
            center = xyz[atoms['C']]
            u, v = xyz[atoms['CA']] - center, xyz[atoms['O']] - center
            direction = -(u / np.linalg.norm(u) + v / np.linalg.norm(v))
            xyz[atoms['OXT']] = center + .125 * direction / np.linalg.norm(direction)
    fixer.positions = xyz * unit.nanometer
    fixer.addMissingHydrogens(PROTOCOL['hydrogen_placement_pH'])
    ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    system = ff.createSystem(fixer.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.2*unit.nanometer, constraints=None)
    restraint = mm.CustomExternalForce(PROTOCOL['restraint_expression'])
    restraint.addGlobalParameter('k', PROTOCOL['backbone_k_kJ_mol_nm2'])
    for name in ('x0', 'y0', 'z0'): restraint.addPerParticleParameter(name)
    coords = fixer.positions.value_in_unit(unit.nanometer)
    for a in fixer.topology.atoms():
        if a.name in BACKBONE and a.element.symbol != 'H':
            restraint.addParticle(a.index, coords[a.index])
    system.addForce(restraint)
    simulation = app.Simulation(fixer.topology, system, mm.VerletIntegrator(.001*unit.picoseconds),
                                mm.Platform.getPlatformByName('CPU'), {'Threads': '2'})
    simulation.context.setPositions(fixer.positions)
    class Reporter(mm.MinimizationReporter):
        def __init__(self):
            super().__init__(); self.start = time.monotonic(); self.timed_out = False; self.iteration = None
        def report(self, iteration, x, grad, args):
            self.iteration = iteration
            self.timed_out = time.monotonic() - self.start > PROTOCOL['minimization_timeout_seconds']
            return self.timed_out
    reporter = Reporter()
    initial = float(simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    simulation.minimizeEnergy(maxIterations=PROTOCOL['max_iterations'], reporter=reporter)
    assert not reporter.timed_out, 'Minimization timed out'
    state = simulation.context.getState(getPositions=True, getEnergy=True, getForces=True)
    final = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    assert math.isfinite(final) and final <= initial, 'Invalid final energy'
    forces = np.array(state.getForces().value_in_unit(unit.kilojoules_per_mole/unit.nanometer))
    out.with_suffix('.relaxation.json').write_text(json.dumps({
        'initial_energy_kJ_mol': initial, 'final_energy_kJ_mol': final,
        'last_reported_iteration': reporter.iteration, 'rms_force_components_kJ_mol_nm': float(np.sqrt(np.mean(forces**2))),
        'protocol': PROTOCOL}, indent=2))
    with out.open('w') as f: app.PDBFile.writeFile(fixer.topology, state.getPositions(), f, keepIds=True)
    save_pdb(read_structure(out), out)


def worker(source, out, queue):
    try:
        relax(source, out); queue.put(None)
    except Exception as e: queue.put(f'{type(e).__name__}: {e}')


def prepare(item):
    source, folder = ROOT/item['source'], ROOT/item['dest']
    folder.mkdir(parents=True, exist_ok=False)
    row = {**item, 'status': 'failed', 'source_sha256': sha(source)}
    try:
        ctx = multiprocessing.get_context('spawn'); queue = ctx.Queue()
        process = ctx.Process(target=worker, args=(source, folder/'prepared.pdb', queue))
        process.start(); process.join(PROTOCOL['hard_timeout_seconds'])
        if process.is_alive():
            process.terminate(); process.join(5)
            if process.is_alive(): process.kill(); process.join()
            raise RuntimeError('Hard preparation timeout')
        assert process.exitcode == 0, f'Preparation worker exited {process.exitcode}'
        error = queue.get(timeout=5)
        if error: raise RuntimeError(error)
        raw, prepared = read_structure(source), read_structure(folder/'prepared.pdb')
        row.update(status='ok', displacement=displacement(raw, prepared),
                   raw_internal_geometry=internal_geometry(source), prepared_internal_geometry=internal_geometry(folder/'prepared.pdb'),
                   prepared_geometry=geometry(prepared, {i:i+189 for i in range(1,195)}, read_structure(ROOT/'data/raw/11QC.cif')),
                   prepared_sha256=sha(folder/'prepared.pdb'))
        save_pdb(prepared, folder/'free_binder.pdb', ProteinSelect({'B'}))
        row['propka'] = {}
        # Audit every prepared candidate, parent and control even if a movement guard fails.
        for key, path in [('free', folder/'free_binder.pdb'), ('bound', folder/'prepared.pdb')]:
            try:
                table = pka_table(path)
                row['propka'][key] = {'status': 'ok', 'titratable_groups': len(table)}
            except Exception as e: row['propka'][key] = {'status': 'failed', 'error': str(e)}
    except Exception as e: row['error'] = f'{type(e).__name__}: {e}'
    (folder/'preparation.json').write_text(json.dumps(row, indent=2))
    return row


def main():
    p = argparse.ArgumentParser(); p.add_argument('--output', type=Path, default=ROOT/'reports/receptor_preparation_v2')
    a = p.parse_args(); out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    source = ROOT/'reports/expansion_his_confirmation/results.json'
    selected = [r for r in json.loads(source.read_text())['rows'] if r['confidence_selected']]
    assert len(selected) == 28
    plan = [{'id': r['id'], 'seed': r['seed'], 'sample': r['sample'], 'source': r['raw_scoring_output']+'/model.pdb',
             'dest': str((out/'models'/f"{r['id']}_seed{r['seed']}").relative_to(ROOT))} for r in selected]
    record = {'protocol': PROTOCOL, 'source_results_sha256': sha(source), 'script_sha256': sha(Path(__file__)),
              'source_hashes': {r['source']: sha(ROOT/r['source']) for r in plan}, 'plan': plan, 'rows': []}
    (out/'plan.json').write_text(json.dumps(record, indent=2))
    (out/'implementation.py.txt').write_bytes(Path(__file__).read_bytes())
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [pool.submit(prepare, item) for item in plan]
        for future in as_completed(futures):
            r = future.result(); record['rows'].append(r)
            (out/'results.json').write_text(json.dumps(record, indent=2))
            print(len(record['rows']), '/', len(plan), r['id'], r['seed'], r['status'],
                  'movement', r.get('displacement', {}).get('pass'), 'protonation', r.get('propka'), flush=True)
    assert len(record['rows']) == len(plan)

if __name__ == '__main__': main()
