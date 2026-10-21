"""Reserve one small mutant_validation run; never retry or exceed the approved allocation."""
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from launch_generation_comparison import gross_cost

ROOT = Path(__file__).resolve().parents[1]


def main():
    run_id = 'mutant-013-6wrw'
    folder = ROOT / 'data/mutant_validation_fixed'
    ledger = ROOT / 'outputs/mutant_validation_fixed_launch.json'
    if ledger.exists():
        raise SystemExit('Confirmation already attempted; no automatic retry.')
    manifest = json.loads((folder / 'manifest.json').read_text())
    for name, sha in manifest['input_sha256'].items():
        assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == sha
    assert manifest['selected_ids'] == ['apical_long_015', 'apical_long_015_D51H_T70H']
    assert len(manifest['cases']) == 5 and manifest['refolding']['seeds'] == [809, 907]
    assert manifest['selection_rule']['max_target_CA_rmsd_to_input_A'] == 2.0
    previous = ROOT / 'reports/confirmation/results.json'
    assert hashlib.sha256(previous.read_bytes()).hexdigest() == manifest['selection_result_sha256']
    assert hashlib.sha256((ROOT / 'reports/confirmation_direction/results.json').read_bytes()).hexdigest() == manifest['direction_result_sha256']
    modal = str(ROOT / '.venv/bin/modal')
    apps = json.loads(subprocess.check_output([modal, 'app', 'list', '--json'], text=True, timeout=45))
    if any(a.get('state') != 'stopped' for a in apps):
        raise SystemExit('Outstanding apps must stop first.')
    billing = json.loads(subprocess.check_output([modal, 'billing', 'summary', '--json'], text=True, timeout=45))
    prior = json.loads((ROOT / 'outputs/confirmation_launch.json').read_text())
    # Keep the entire previous reservation to avoid spending lagging billing twice.
    failed = ROOT / 'outputs/mutant012_failed/mutant-012-6wrw'
    failed_status = json.loads((failed / 'run_status.json').read_text())
    assert failed_status['statuses'] == [{'seed': 809, 'returncode': 1}]
    assert not list(failed.glob('seed-*/fold_out_npz/*.npz'))
    rate = .000542 + 4*.0000131 + 32*.00000222
    failed_allowance = (failed_status['elapsed_seconds'] + 600 + 2) * rate + .25
    gross = max(gross_cost(billing), 7.78736345, prior['prior_gross_high_water_usd'] + prior['reserved_usd'] + failed_allowance)
    reserve = 1.5 - failed_allowance
    resource_bound = (750 + 120 + 2) * rate
    assert resource_bound + .25 < reserve and gross + reserve <= 30
    assert 4 + prior['reserved_usd'] + failed_allowance + reserve <= prior['deep_evaluation_allocation_usd'] == 7
    files = {p.name: p.read_bytes() for p in folder.iterdir() if p.suffix in {'.cif', '.npz', '.json', '.a3m'}}
    files['control_refold_worker.py'] = (ROOT / 'scripts/control_refold_worker.py').read_bytes()
    report = ROOT / 'reports/mutant_validation/corrected_run_record'
    report.mkdir(parents=True, exist_ok=True)
    record = {'run_id': run_id, 'input_dir': 'data/mutant_validation_fixed',
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'prior_gross_high_water_usd': gross, 'reserved_usd': reserve,
              'authorized_total_usd': 30, 'deep_evaluation_allocation_usd': 7, 'shared_experiment_cap_usd': 1.5, 'failed_attempt_allowance_usd': failed_allowance,
              'conservative_resource_bound_usd': resource_bound,
              'runtime_bound': {'work_seconds': 650, 'function_timeout_seconds': 750,
                                'startup_timeout_seconds': 120, 'max_containers': 1, 'retries': 0},
              'input_sha256': {k: hashlib.sha256(v).hexdigest() for k, v in files.items()},
              'status': 'launched_or_attempted'}
    with ledger.open('x') as f:
        json.dump(record, f, indent=2)
    for name, value in [('billing_before', billing), ('apps_before', apps), ('launch_ledger', record)]:
        (report / f'{name}.json').write_text(json.dumps(value, indent=2))
    (report / 'manifest.json').write_bytes(files['manifest.json'])
    for name in ('control_refold_worker.py', 'modal_mutant_validation_fixed.py', 'launch_mutant_validation_fixed.py',
                 'prepare_mutant_validation.py', 'analyze_mutant_validation.py', 'target_integrity_selection.py', 'prepare_controls.py', 'repair_mutant_input_names.py'):
        (report / name).write_bytes((ROOT / 'scripts' / name).read_bytes())
    record['analysis_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in report.glob('*.py')}
    ledger.write_text(json.dumps(record, indent=2))
    (report / 'launch_ledger.json').write_text(json.dumps(record, indent=2))
    result = subprocess.run([modal, 'run', '--detach', str(ROOT / 'scripts/modal_mutant_validation_fixed.py'),
                             '--run-id', run_id, '--input-dir', 'data/mutant_validation_fixed'],
                            cwd=ROOT, env=dict(os.environ, TFR1_LAUNCH_RESERVED='1'))
    record['client_exit_code'] = result.returncode
    ledger.write_text(json.dumps(record, indent=2))
    (report / 'launch_ledger.json').write_text(json.dumps(record, indent=2))
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
