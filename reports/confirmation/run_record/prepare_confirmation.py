"""Freeze a prospective target-conformation filter and an independent 50-sample test."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    previous = ROOT / 'reports/broad_validation/results.json'
    result = json.loads(previous.read_text())
    chosen = [p['id'] for p in result['candidate_summaries'] if p['both_seeds_candidate_screen']]
    assert chosen == ['apical_short_008', 'apical_long_015']
    source = ROOT / 'data/broad_validation'
    manifest = json.loads((source / 'manifest.json').read_text())
    cases = [c for c in manifest['cases'] if c['id'] in chosen or c['arm'] == 'control']
    assert len(cases) == 5
    out = ROOT / 'data/confirmation'
    out.mkdir(exist_ok=False)
    for case in cases:
        for ext in ('.pdb', '.cif', '.npz'):
            name = case['id'] + ext
            assert hashlib.sha256((source / name).read_bytes()).hexdigest() == manifest['input_sha256'][name]
            (out / name).write_bytes((source / name).read_bytes())
    manifest.update(
        scope='Prospective confirmation of two previously held leads; original batch decision unchanged',
        prepared_utc=datetime.now(timezone.utc).isoformat(),
        source_run=result['run'], cases=cases, selected_ids=chosen,
        selection_result_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
        refolding={**manifest['refolding'], 'seeds': [607, 709], 'diffusion_samples_per_seed': 5},
        max_prediction_seconds=900, experiment_cap_usd=1.5,
        selection_rule={
            'max_target_CA_rmsd_to_input_A': 2.0,
            'alignment': 'Rigid least-squares fit using all 194 corresponding target C-alpha atoms; no binder atoms',
            'completeness': 'Exact target sequence and all 194 finite target CA coordinates required',
            'ranking': 'Among eligible samples, maximize 0.8*ipTM + 0.2*pTM; tie -> lowest sample index',
            'no_eligible_samples': 'Case fails; no fallback to a rejected sample',
            'scope': 'Identical rule for positive, negative and novel cases; no binder geometry or epitope used for selection',
            'rationale': 'Target distortion made global-frame pose validation unreliable in a previous control. The 2 A scale matches the existing structural-screen tolerance; this is a new heuristic, not a validated physical boundary.',
        },
    )
    manifest['screen'] = {**manifest['screen'], 'replicates': 'Both fresh seeds; target-conformation-filtered confidence selection',
                          'decision': 'Controls must pass in both seeds; retain all samples and failures. No affinity, mouse-binding or pH-switch claim.'}
    manifest['input_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    report = ROOT / 'reports/confirmation'
    report.mkdir(exist_ok=False)
    (report / 'plan.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'cases': 5, 'predictions': 50, 'selected_ids': chosen, 'seeds': [607, 709]}, indent=2))


if __name__ == '__main__':
    main()
