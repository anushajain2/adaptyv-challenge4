import copy
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT
from run_local_mpnn import eligible_sequence, select_diverse
from audit_preparation_motion import motion_audit
from prepare_controls import read_structure


def test_fixed_histidine_and_target_context_policy():
    plan=json.loads((ROOT/'data/local_redesign/plan.json').read_text())
    for f in plan['families']:
        assert f['histidine_position'] not in f['mutable_positions']
        assert not eligible_sequence(f['sequence'],f)
        seq=list(f['sequence']);seq[f['histidine_position']-1]='A'
        with pytest.raises(AssertionError):eligible_sequence(''.join(seq),f)
        assert all(f['sequence'][p-1] not in 'GPCHDEKR' for p in f['mutable_positions'])


def test_diversity_does_not_fill_with_near_duplicates():
    rows=[{'sequence':'AAAA','ensemble_nll':1.},{'sequence':'AAAV','ensemble_nll':1.1},
          {'sequence':'AAVV','ensemble_nll':1.2},{'sequence':'AAVV','ensemble_nll':2.}]
    assert [r['sequence'] for r in select_diverse(rows)]==['AAAA','AAVV']


def test_remote_oxygen_and_interface_oxygen_are_distinguished():
    audit=json.loads((ROOT/'reports/local_redesign/movement_audit.json').read_text())
    pos=next(r for r in audit['rows'] if r['id']=='2ds25_5_6wrw' and r['seed']==2609)
    mutant=next(r for r in audit['rows'] if r['id']=='diversity_long_037_S72H' and r['seed']==2503)
    assert not pos['historical_movement_pass'] and not mutant['historical_movement_pass']
    assert pos['interface_sensitivity_pass'] and not mutant['interface_sensitivity_pass']
    raw=read_structure(ROOT/pos['raw']);same=motion_audit(raw,copy.deepcopy(raw))
    assert same['interface_sensitivity_pass'] and same['all_backbone_over_1A']==[]


def test_panel_preserves_target_and_fixed_binder_positions():
    path=ROOT/'data/local_redesign_panel/manifest.json'
    if not path.exists():pytest.skip('Panel not constructed yet')
    m=json.loads(path.read_text());plan=m['prospective_plan'];families={f['id']:f for f in plan['families']}
    assert m['refolding']['seeds']==[2909,3001]
    for c in m['cases']:
        if c.get('local_redesign'):assert eligible_sequence(c['binder_sequence'],families[c['parent_id']])
        assert c['target_sequence']==m['target_sequence']
