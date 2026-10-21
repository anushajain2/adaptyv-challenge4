"""Verify the completed round and record the later draft admission separately."""
import csv
import json
from datetime import datetime, timezone

from novelty_search import ROOT, sha, save


def main():
    out = ROOT / 'reports/diversity_confirmation'
    draft = ROOT / 'reports/submission_draft'
    old = ROOT / 'reports/novelty_review'
    prior = {r['name']: r['sequence'] for r in csv.DictReader((old / 'reviewed_submission_draft.csv').open())}
    current = {r['name']: r['sequence'] for r in csv.DictReader((draft / 'submission_draft.csv').open())}
    added = {'diversity_long_037', 'diversity_short_055'}
    assert len(prior) == 6 and len(current) == 8
    assert set(current) - set(prior) == added
    assert all(current[k] == v for k, v in prior.items())
    audit = json.loads((draft / 'audit.json').read_text())
    assert audit['parent_lineage_count'] == 7
    result = json.loads((out / 'results.json').read_text())
    assert len(result['rows']) == 80 and result['controls_pass']
    assert {c['id'] for c in result['candidate_summaries'] if c['all_selected_seeds_full_screen']} == added
    for folder in ['diversity_generation', 'diversity_confirmation']:
        v = json.loads((ROOT / 'reports' / folder / 'run_record/verification.json').read_text())
        assert v['all_inputs_and_sources_match_launch'] and v['all_cloud_stages_successful']
        for path, digest in v['files'].items():
            assert sha(ROOT / path) == digest, path
    for folder in ['novelty_review', 'diversity_novelty']:
        manifest = json.loads((ROOT / 'reports' / folder / 'manifest.json').read_text())
        for group in ['files', 'scripts']:
            for path, digest in manifest[group].items():
                assert sha(ROOT / path) == digest, path
    billing = json.loads((out / 'run_record/billing_after.json').read_text())
    gross = max(float(billing['metered_cost']), float(billing['metered_cost_breakdown']['ephemeral_apps']))
    assert gross < 30
    apps = json.loads((out / 'run_record/apps_after.json').read_text())
    assert apps and all(a['state'] == 'stopped' for a in apps)
    save(out / 'submission_update.json', {
        'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'added_ids': sorted(added), 'original_six_sequences_unchanged': True,
        'candidate_count': 8, 'parent_lineage_count': 7,
        'files': {str(p.relative_to(ROOT)): sha(p) for p in [
            old / 'reviewed_submission_draft.csv', old / 'reviewed_submission_draft.fasta',
            draft / 'submission_draft.csv', draft / 'submission_draft.fasta',
            out / 'summary.json', ROOT / 'reports/diversity_novelty/summary.json']},
        'challenge_submitted': False,
    })
    paths = [p for folder in ['diversity_generation', 'diversity_confirmation']
             for p in (ROOT / 'reports' / folder).rglob('*')
             if p.is_file() and p.name != 'round_evidence_audit.json']
    save(out / 'round_evidence_audit.json', {
        'status': 'COMPLETE', 'generated': 100, 'unique_generated': 99,
        'initial_predictions_including_controls': 103, 'confirmation_predictions': 80,
        'initial_full_passes': 5, 'confirmed_parents': sorted(added),
        'all_selected_confirmation_controls_pass': True,
        'draft_candidates': 8, 'parent_lineages': 7,
        'recorded_gross_high_water_usd': gross, 'budget_cap_usd': 30,
        'billing_may_lag': True, 'all_apps_stopped_in_archived_snapshot': True,
        'challenge_submitted': False, 'upstream_issue_or_PR_filed': False,
        'files': {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)},
        'review_manifest_sha256': sha(ROOT / 'reports/diversity_novelty/manifest.json'),
        'draft_manifest_sha256': sha(draft / 'manifest.json'),
        'script_sha256': sha(ROOT / 'scripts/finalize_diversity_round.py'),
    })
    print('Verified completed round, preserved original six sequences, and recorded eight-sequence draft admission.')


if __name__ == '__main__':
    main()
