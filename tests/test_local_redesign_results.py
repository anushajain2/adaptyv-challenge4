"""Evidence consistency checks for completed local-redesign runs."""
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT,sha
from analyze_receptor_preparation_v2 import direction_metrics
from direction_check import histidines
from summarize_direction import terminal_distance


def completed():
    folder=ROOT/'reports/local_redesign_panel'
    if not (folder/'direction.json').exists():pytest.skip('Predictions/preparation not complete yet')
    return folder,json.loads((folder/'direction.json').read_text())


def test_all_predictions_and_selected_poses_accounted_for():
    folder,direction=completed();raw=json.loads((folder/'results.json').read_text())
    m=json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text())
    expected={(c['id'],seed,sample) for c in m['cases'] for seed in [2909,3001] for sample in range(5)}
    assert len(raw['rows'])==len(expected)==130
    assert {(r['id'],r['seed'],r['sample']) for r in raw['rows']}==expected
    selected={(r['id'],r['seed'],r['sample']) for r in raw['rows'] if r['confidence_selected']}
    assert {(r['id'],r['seed'],r['sample']) for r in direction['rows']}==selected
    assert direction['raw_results_sha256']==sha(folder/'results.json')


def test_pka_results_and_promotion_obey_frozen_gate():
    folder,direction=completed();m=json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    for r in direction['rows']:
        if r['direction']['status']!='diagnostic':continue
        path=ROOT/r['source'];assert sha(path)==r['source_sha256'];positions=cases[r['id']]['histidine_positions']
        f=histidines(json.loads(path.with_name('free_binder.pka.json').read_text()),'B',positions)
        b=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',positions)
        check=direction_metrics(f,b,terminal_distance(path,positions))
        assert all(r['direction'][key]==value for key,value in check.items())
    for v in direction['variants']:
        expected=(direction['raw_controls_pass'] and direction['strict_prepared_controls_pass'] and
                  v['strict_structure_and_direction_both'] and v['baseline_strict_eligible_both'] and
                  v['unmodified_parent_strict_eligible_both'] and v['prospective_improvement_heuristic'])
        assert v['fully_qualified_for_followup']==expected
    assert all(v['fully_qualified_for_followup'] for v in direction['variants'] if v['id'] in direction['followup_ids'])
