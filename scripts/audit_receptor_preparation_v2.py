"""Verify revised evidence, source identity and preservation of raw exclusions."""
import importlib.metadata
import json
import platform
from pathlib import Path
from prepare_controls import ROOT,read_structure,sequence
from receptor_preparation_v2 import sha,displacement,internal_geometry
from analyze_receptor_preparation_v2 import direction_metrics
from direction_check import histidines
from summarize_direction import terminal_distance


def main():
    out=ROOT/'reports/receptor_preparation_v2'
    prep=json.loads((out/'results.json').read_text());analysis=json.loads((out/'analysis.json').read_text())
    assert sha(ROOT/'scripts/receptor_preparation_v2.py')==prep['script_sha256']==sha(out/'implementation.py.txt')
    assert sha(ROOT/'reports/expansion_his_confirmation/results.json')==prep['source_results_sha256']
    assert sha(out/'results.json')==analysis['source_preparation_sha256']
    assert sha(ROOT/'reports/expansion_his_direction/results.json')==analysis['source_v1_direction_sha256']
    expected={(r['id'],r['seed'],r['sample']) for r in prep['plan']}
    assert len(expected)==len(prep['rows'])==len(analysis['rows'])==28
    assert expected=={(r['id'],r['seed'],r['sample']) for r in analysis['rows']}
    lookup={(r['id'],r['seed']):r for r in analysis['rows']};direction_count=0
    for r in prep['rows']:
        src=ROOT/r['source'];folder=ROOT/r['dest'];dest=folder/'prepared.pdb'
        assert sha(src)==r['source_sha256']==prep['source_hashes'][r['source']]
        assert r['status']=='ok' and sha(dest)==r['prepared_sha256']
        raw,prepared=read_structure(src),read_structure(dest)
        for c in ('A','B'):assert sequence(raw[0][c])==sequence(prepared[0][c])
        assert displacement(raw,prepared)==r['displacement']
        assert internal_geometry(dest)==r['prepared_internal_geometry']
        for name in ('free_binder','prepared'):
            assert json.loads((folder/f'{name}.propka_warnings.json').read_text())==[]
            assert len(json.loads((folder/f'{name}.pka.json').read_text()))>0
        scored=lookup[(r['id'],r['seed'])]
        expected_eligible=bool(scored['raw_pass'] and r['displacement']['pass'] and r['prepared_internal_geometry']['pass'] and scored['prepared_screen']['pass'] and all(v['status']=='ok' for v in r['propka'].values()))
        assert scored['prepared_eligible']==expected_eligible
        if 'direction' not in scored:continue
        d=scored['direction'];assert d['status']=='ok' and expected_eligible
        cases=json.loads((ROOT/'data/expansion_his_confirmation/manifest.json').read_text())['cases']
        case=next(c for c in cases if c['id']==r['id']);positions=[int(s[1:-1]) for s in case['mutations']]
        free=histidines(json.loads((folder/'free_binder.pka.json').read_text()),'B',positions)
        bound=histidines(json.loads((folder/'prepared.pka.json').read_text()),'B',positions)
        calculated=direction_metrics(free,bound,terminal_distance(dest,positions))
        assert all(d[k]==v for k,v in calculated.items());direction_count+=1
    assert direction_count==10
    assert not analysis['controls_prepared_pass'] and not any(v['direction_pass_both_seeds'] for v in analysis['variants'])
    record={'status':'pass','models_verified':28,'free_bound_tables_verified':56,'direction_calculations_recomputed':direction_count,
        'raw_exclusions_preserved':True,'preparation_implementation_matches_frozen_snapshot':True,
        'source_results_unchanged':True,'python':platform.python_version(),
        'versions':{name:importlib.metadata.version(name) for name in ('openmm','pdbfixer','propka','numpy','scipy','biopython')},
        'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='evidence_audit.json'}}
    (out/'evidence_audit.json').write_text(json.dumps(record,indent=2));print('Verified 28 prepared models, 56 PROPKA tables, 10 direction calculations; original exclusions preserved.')

if __name__=='__main__':main()
