"""Final provenance and calculation checks after all local work finishes."""
import json
from prepare_controls import ROOT,read_structure,sequence
from prepare_tag_models import displacement,internal_geometry,sha
from direction_check import histidines
from analyze_receptor_preparation_v2 import direction_metrics
from summarize_direction import terminal_distance

def main():
    out=ROOT/'reports/tag_comparison';m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    raw=json.loads((out/'results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text());direction=json.loads((out/'direction.json').read_text())
    assert len(raw['rows'])==140 and len(prep['rows'])==28 and len(direction['variants'])==4
    assert prep['script_sha256']==sha(ROOT/'scripts/prepare_tag_models.py')==sha(out/'preparation/implementation.py.txt')
    assert prep['source_results_sha256']==sha(out/'results.json')==direction['source_results_sha256']
    assert direction['preparation_results_sha256']==sha(out/'preparation/results.json')
    prepared={(r['id'],r['seed']):r for r in prep['rows']};checked=0
    for p in prep['rows']:
        source=ROOT/p['source'];assert sha(source)==p['source_sha256']==prep['source_hashes'][p['source']]
        if p['status']!='ok':continue
        path=ROOT/p['dest']/'prepared.pdb';assert sha(path)==p['prepared_sha256']
        a,b=read_structure(source),read_structure(path)
        assert sequence(b[0]['B'])==cases[p['id']]['binder_sequence']
        assert displacement(a,b,p['core_length'])==p['displacement']
        assert internal_geometry(path)==p['prepared_internal_geometry']
    for variant in direction['variants']:
        c=cases[variant['id']];positions=[int(x[1:-1]) for x in c['mutations']]
        for row in variant['models']:
            if row['status']!='ok':
                assert not row['direction_pass'] and not row['tag_independent_direction_pass'];continue
            p=prepared[(row['id'],row['seed'])];folder=ROOT/p['dest'];path=folder/'prepared.pdb'
            assert variant['raw_core_pass_both_seeds'] and p['displacement']['pass'] and p['prepared_internal_geometry']['pass']
            assert all(x['status']=='ok' for x in p['propka'].values())
            free=histidines(json.loads((folder/'free_binder.pka.json').read_text()),'B',positions)
            bound=histidines(json.loads((folder/'prepared.pka.json').read_text()),'B',positions)
            calculated=direction_metrics(free,bound,terminal_distance(path,positions))
            assert all(row[k]==v for k,v in calculated.items());checked+=1
    bundle=json.loads((out/'review_bundle/manifest.json').read_text())
    for model in bundle['models']:
        folder=out/'review_bundle'/f"{model['id']}_seed{model['seed']}"
        for name,digest in model['files'].items():assert sha(folder/name)==digest
    record={'status':'pass','predictions':140,'selected_models':28,'direction_calculations_recomputed':checked,
        'raw_exclusions_preserved':True,'full_suffix_preserved_for_preparation_and_propka':True,
        'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='round_evidence_audit.json'}}
    (out/'round_evidence_audit.json').write_text(json.dumps(record,indent=2));print('Round audit passed;',checked,'admitted direction calculations verified')

if __name__=='__main__':main()
