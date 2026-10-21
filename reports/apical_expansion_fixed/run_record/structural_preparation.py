"""Frozen, bounded preparation diagnostic; never overwrite raw model decisions."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from direction_check import bounded_relax
from prepare_controls import ROOT, read_structure, geometry, sequence

PROTOCOL = {
    'version': 1, 'method': 'Amber14/GBn2; fixed receptor; binder backbone restrained at 1000 kJ/mol/nm^2',
    'seed': 43, 'pH_for_hydrogen_placement': 7.4, 'max_iterations': 500,
    'hard_timeout_seconds': 180, 'max_binder_backbone_rms_displacement_A': .5,
    'max_binder_backbone_atom_displacement_A': 1., 'max_target_atom_displacement_A': .01,
    'policy': 'Apply to every selected candidate and control, irrespective of raw outcome. Preserve raw decisions; no alternate sample or repeat relaxation after failure. Prepared evidence is a separate prospective diagnostic, not retrospective qualification.',
}

def displacement(raw, prepared):
    assert all(sequence(raw[0][c]) == sequence(prepared[0][c]) for c in ('A','B'))
    values = {}
    for chain in ('A','B'):
        ds = []
        for r in raw[0][chain]:
            for a in r:
                if a.element == 'H' or (chain == 'B' and a.name not in {'N','CA','C','O'}):
                    continue
                ds.append(float(np.linalg.norm(a.coord-prepared[0][chain][r.id][a.name].coord)))
        values[chain] = {'rms_A': float(np.sqrt(np.mean(np.square(ds)))), 'max_A': max(ds)}
    values['pass'] = (values['A']['max_A'] <= .01 and values['B']['rms_A'] <= .5 and values['B']['max_A'] <= 1.)
    return values

def prepare(source, dest):
    dest.mkdir(parents=True, exist_ok=False)
    raw = read_structure(source)
    result = {'source': str(source.relative_to(ROOT)), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'protocol': PROTOCOL}
    try:
        bounded_relax(source, dest/'prepared.pdb', (), PROTOCOL['seed'])
        prepared = read_structure(dest/'prepared.pdb')
        result.update(status='ok', displacement=displacement(raw,prepared))
        full = read_structure(ROOT/'data/raw/11QC.cif')
        for key, s in [('raw_geometry',raw), ('prepared_geometry',prepared)]:
            result[key] = geometry(s,{i:i+189 for i in range(1,195)},full)
        result['prepared_sha256'] = hashlib.sha256((dest/'prepared.pdb').read_bytes()).hexdigest()
    except Exception as e:
        result.update(status='failed', error=f'{type(e).__name__}: {e}')
    (dest/'preparation.json').write_text(json.dumps(result,indent=2))
    return result

def main():
    out=ROOT/'reports/structural_preparation';out.mkdir(exist_ok=False)
    (out/'protocol.json').write_text(json.dumps(PROTOCOL,indent=2))
    mutant=json.loads((ROOT/'reports/mutant_validation/results.json').read_text())
    confirmation=json.loads((ROOT/'reports/confirmation/results.json').read_text())
    selected=[r for r in mutant['rows'] if r['confidence_selected']]
    selected += [r for r in confirmation['rows'] if r['confidence_selected'] and r['id']=='apical_short_008']
    plan=[{'id':r['id'],'seed':r['seed'],'sample':r['sample'],'source':r['raw_scoring_output']+'/model.pdb'} for r in selected]
    (out/'plan.json').write_text(json.dumps(plan,indent=2))
    results=[]
    for item in plan:
        print('Preparing',item['id'],item['seed'],flush=True)
        r=prepare(ROOT/item['source'],out/'models'/f"{item['id']}_seed{item['seed']}")
        results.append({**item,**r})
        (out/'results.json').write_text(json.dumps({'protocol':PROTOCOL,'rows':results},indent=2))
    print(json.dumps({'completed':len(results),'prepared':sum(r['status']=='ok' for r in results)},indent=2))

if __name__=='__main__':main()
