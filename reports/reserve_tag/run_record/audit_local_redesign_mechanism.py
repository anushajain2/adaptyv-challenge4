"""Decompose existing PROPKA results; never changes frozen eligibility decisions."""
import json
import logging
from propka.run import single
from novelty_search import ROOT,sha,save


def group_terms(path,position,expected_pka):
    warnings=[]
    class Capture(logging.Handler):
        def emit(self,record):warnings.append(record.getMessage())
    handler=Capture(logging.WARNING);logger=logging.getLogger('propka');logger.addHandler(handler)
    try:molecule=single(str(path),write_pka=False)
    finally:logger.removeHandler(handler)
    assert not any('Missing atoms or failed protonation' in x for x in warnings)
    groups=[g for g in molecule.conformations['AVR'].groups if g.type=='HIS' and g.atom.chain_id=='B' and g.atom.res_num==position]
    assert len(groups)==1;group=groups[0]
    assert abs(group.pka_value-expected_pka)<1e-9
    determinants={kind:[{'group':d.group.label,'value_pka_units':float(d.value)} for d in values] for kind,values in group.determinants.items()}
    terms={'volume_desolvation':float(group.energy_volume),'local':float(group.energy_local),
           **{kind:sum(d['value_pka_units'] for d in values) for kind,values in determinants.items()}}
    assert abs(group.model_pka+sum(terms.values())-group.pka_value)<1e-8
    return {'pka':float(group.pka_value),'model_pka':float(group.model_pka),'buried_fraction':float(group.buried),
            'contributions_pka_units':terms,'determinants':determinants,'source':str(path.relative_to(ROOT)),
            'source_sha256':sha(path),'warnings':warnings}


def main():
    folder=ROOT/'reports/local_redesign_panel';source=folder/'direction.json';direction=json.loads(source.read_text())
    cases={c['id']:c for c in json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text())['cases']};rows=[]
    for row in direction['rows']:
        if row['direction']['status']!='diagnostic':continue
        positions=cases[row['id']]['histidine_positions'];assert len(positions)==1
        bound=ROOT/row['source'];free=bound.with_name('free_binder.pdb')
        f=group_terms(free,positions[0],row['direction']['free_pkas'][0]);b=group_terms(bound,positions[0],row['direction']['bound_pkas'][0])
        delta={k:b['contributions_pka_units'][k]-f['contributions_pka_units'][k] for k in f['contributions_pka_units']}
        assert abs(sum(delta.values())-row['direction']['shifts'][0])<1e-8
        rows.append({'id':row['id'],'seed':row['seed'],'strict_eligible':row['strict_eligible'],'free':f,'bound':b,'bound_minus_free_contributions_pka_units':delta})
    save(folder/'mechanism_audit.json',{'source_direction_sha256':sha(source),'rows':rows,
        'interpretation':'Post-result descriptive decomposition of the same PROPKA model, not independent validation or a new selection rule. Terms are pKa offsets, not physical binding energies. Cropped target and one fixed conformation per model remain limitations.'})
    print('Reproduced and decomposed',len(rows),'existing free/bound histidine comparisons.')


if __name__=='__main__':main()
