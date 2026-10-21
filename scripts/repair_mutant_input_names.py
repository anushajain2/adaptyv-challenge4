"""Preserve failed inputs and create an alias accepted by the pinned filename parser."""
import ast
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OLD = 'apical_long_015__D51H_T70H'
NEW = 'apical_long_015_D51H_T70H'


def pinned_target_regex():
    path = ROOT / 'vendor/boltzgen/src/boltzgen/task/predict/data_from_generated.py'
    tree = ast.parse(path.read_text())
    values = [ast.literal_eval(n.value) for n in ast.walk(tree)
              if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name) and n.target.id == 'target_id_regex']
    assert len(values) == 1
    return values[0]


def main():
    source = ROOT / 'data/mutant_validation'
    out = ROOT / 'data/mutant_validation_fixed'
    out.mkdir(exist_ok=False)
    m = json.loads((source / 'manifest.json').read_text())
    hashes = m['input_sha256']
    for name, sha in hashes.items():
        data = (source / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == sha
        (out / name.replace(OLD, NEW)).write_bytes(data)
    for c in m['cases']:
        if c['id'] == OLD:
            c['variant_label'] = OLD
            c['id'] = NEW
    m['selected_ids'] = [x.replace(OLD, NEW) for x in m['selected_ids']]
    m['id_aliases'] = {NEW: OLD}
    m['max_prediction_seconds'] = 650
    regex = pinned_target_regex()
    assert not re.fullmatch(regex, OLD) and re.fullmatch(regex, NEW)
    assert all(re.fullmatch(regex, c['id']) for c in m['cases'])
    m['input_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}
    (out / 'manifest.json').write_text(json.dumps(m, indent=2))
    record = {'reason': 'Double underscore creates an empty component rejected by the pinned target_id_regex.',
              'filename_regex': regex, 'alias': {NEW: OLD},
              'sequence_coordinates_masks_and_scientific_thresholds_unchanged': True,
              'all_cif_npz_pdb_bytes_identical_under_alias': all(
                  (out / name.replace(OLD, NEW)).read_bytes() == (source / name).read_bytes() for name in hashes),
              'failed_run': 'mutant-012-6wrw', 'corrected_run': 'mutant-013-6wrw',
              'corrected_runtime': {'work_seconds': 650, 'function_timeout_seconds': 750, 'startup_timeout_seconds': 120}}
    (ROOT / 'reports/mutant_validation/input_name_repair.json').write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
