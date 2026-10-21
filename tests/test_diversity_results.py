import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT, sha


def test_high_passing_sample_count_does_not_rescue_failed_selected_pose():
    r=json.loads((ROOT/'reports/diversity_confirmation/results.json').read_text())
    assert len(r['rows'])==80 and len(r['control_checks'])==6 and r['controls_pass']
    cs={c['id']:c for c in r['candidate_summaries']}
    assert {c['id'] for c in cs.values() if c['all_selected_seeds_full_screen']}=={
        'diversity_short_055','diversity_long_037'}
    for name in ['diversity_short_004','diversity_short_023']:
        c=cs[name]
        assert c['passing_samples']==8 and not c['all_selected_seeds_full_screen']
        assert len(c['selected_models'])==2
        assert any('clash_clearance' in m['screen']['failed'] for m in c['selected_models'])


def test_uniform_preparation_and_bundle_keep_controls_and_failed_candidates():
    out=ROOT/'reports/diversity_confirmation'
    raw=json.loads((out/'results.json').read_text())
    prep=json.loads((out/'preparation/results.json').read_text())
    bundle=json.loads((out/'review_bundle/manifest.json').read_text())
    assert len(prep['rows'])==len(bundle['models'])==16
    assert {(r['id'],r['seed'],r['sample']) for r in prep['rows']}=={
        (r['id'],r['seed'],r['sample']) for r in raw['rows'] if r['confidence_selected']}
    assert all(r['status']=='ok' for r in prep['rows'])
    for model in bundle['models']:
        folder=out/'review_bundle'/f"{model['id']}_seed{model['seed']}"
        for name,digest in model['files'].items():assert sha(folder/name)==digest
    summary=json.loads((out/'summary.json').read_text())
    retained=[c for c in summary['candidates'] if c['two_seed_raw_screen_pass']]
    assert all(c['two_pose_preparation_pass'] and not c['sequence_audit']['His_positions'] for c in retained)
    diversity=json.loads((out/'confirmed_diversity.json').read_text())
    assert len(diversity['comparisons_to_current_draft'])==48
    assert len(diversity['between_new_parents'])==4 and diversity['all_reciprocal_TM_below_0_8']


def test_added_parent_review_and_payload_provenance():
    out=ROOT/'reports/diversity_novelty'
    summary=json.loads((out/'summary.json').read_text())
    assert summary['status']=='COMPLETE_WITH_DATABASE_LIMITS'
    assert all(summary['controls'].values())
    assert set(summary['reviewed_candidate_ids'])=={'diversity_short_055','diversity_long_037'}
    assert all(c['swissprot_close_flags']==0 for c in summary['candidates'])
    manifest=json.loads((out/'manifest.json').read_text())
    for group in ['files','scripts']:
        for path,digest in manifest[group].items():assert sha(ROOT/path)==digest
    update=json.loads((ROOT/'reports/diversity_confirmation/submission_update.json').read_text())
    assert update['original_six_sequences_unchanged'] and update['candidate_count']==8
    for path,digest in update['files'].items():assert sha(ROOT/path)==digest


def test_disulfide_assumption_matches_archived_coordinates():
    from prepare_controls import read_structure
    audit=json.loads((ROOT/'reports/diversity_confirmation/disulfide_audit.json').read_text())
    assert len(audit['rows'])==4
    for row in audit['rows']:
        path=ROOT/row['source']
        assert sha(path)==row['source_sha256']
        chain=read_structure(path)[0]['B']
        distance=float(chain[31]['SG']-chain[42]['SG'])
        assert abs(distance-row['C31_C42_SG_distance_A'])<1e-5
        assert {frozenset(b) for b in row['OpenMM_inferred_SG_SG_bonds']}=={frozenset(['31','42'])}
    assert 'not experimental confirmation' in audit['interpretation']
