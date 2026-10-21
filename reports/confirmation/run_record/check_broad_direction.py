"""Bounded His-pair checks on parents that passed both fresh prediction seeds."""
import argparse
import hashlib
import json
from pathlib import Path

from direction_check import enumerate_pairs, read_structure, scan
from prepare_target import ROOT, write_csv
from summarize_direction import terminal_distance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation', type=Path, default=ROOT / 'reports/broad_validation/results.json')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/broad_direction')
    parser.add_argument('--report', type=Path, default=ROOT / 'reports/broad_direction')
    args = parser.parse_args()
    source = args.validation
    results = json.loads(source.read_text())
    assert results['controls_pass'], 'Do not engineer parents after failed controls'
    parents = [p for p in results['candidate_summaries'] if p['provisional_computational_parent']]
    parents.sort(key=lambda p: -min(r['metrics']['design_to_target_iptm'] for r in p['selected_models']))
    parents = parents[:3]
    output = args.output
    report = args.report
    prediction_seeds = sorted({r['seed'] for r in results['rows']})
    assert len(prediction_seeds) == 2
    canonical_seed = prediction_seeds[0]
    output.mkdir(exist_ok=False)
    report.mkdir(parents=True, exist_ok=False)
    plan = {
        'validation_results_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'parent_ids': [p['id'] for p in parents],
        'max_pairs_per_parent': 2,
        'prediction_seeds': prediction_seeds,
        'local_rebuild_seeds': [11, 29],
        'pair_selection': f'Enumerate on the selected seed-{canonical_seed} pose, then test identical mutations on both selected prediction poses.',
        'min_distance_to_artificial_target_termini_A': 10.0,
        'decision': 'All four models must satisfy the existing direction criteria and the artificial-terminus distance buffer. This is a direction hypothesis, not mutant fold or affinity validation.',
    }
    (report / 'plan.json').write_text(json.dumps(plan, indent=2))
    rows, pairs = [], []
    for parent in parents:
        models = sorted(parent['selected_models'], key=lambda r: r['seed'])
        assert [r['seed'] for r in models] == prediction_seeds
        assert parent['protocol'] == 'protein-anything' and all(r['screen']['pass'] for r in models)
        paths = {r['seed']: ROOT / r['raw_scoring_output'] / 'model.pdb' for r in models}
        # run_ipsae saves the full predicted structure, including restored sidechains.
        assert all(p.exists() for p in paths.values())
        options = enumerate_pairs(read_structure(paths[canonical_seed]), 2, set(range(1, len(parent['sequence']) + 1)))
        for option in options:
            positions = option['positions']
            label = f"H{positions[0]}_H{positions[1]}"
            group = []
            for model in models:
                seed = model['seed']
                dest = output / parent['id'] / f'seed-{seed}' / label
                scanned = scan(paths[seed], dest, explicit_pair=positions, seeds=(11, 29))
                for row in scanned:
                    row.update(parent=parent['id'], prediction_seed=seed,
                               prediction_sample=model['sample'],
                               source_model_sha256=hashlib.sha256(paths[seed].read_bytes()).hexdigest())
                    pdb = dest / row['variant'] / 'complex.pdb'
                    row['distance_to_artificial_target_termini_A'] = terminal_distance(pdb, positions) if row['status'] == 'ok' else None
                    row['evidence_directory'] = str((dest / row['variant']).relative_to(ROOT))
                group.extend(scanned)
                rows.extend(scanned)
                write_csv(report / 'all_models.csv', rows)
            seq = list(parent['sequence'])
            mutations = [f'{seq[i-1]}{i}H' for i in positions]
            for i in positions:
                seq[i-1] = 'H'
            complete = len(group) == 4 and all(r['status'] == 'ok' for r in group)
            pairs.append({
                'id': parent['id'] + '__' + '_'.join(mutations), 'parent': parent['id'],
                'mutations': mutations, 'sequence': ''.join(seq),
                'models_attempted': len(group), 'models_successful': sum(r['status'] == 'ok' for r in group),
                'both_pKas_lower_in_all_four': complete and all(r['shift_1'] < 0 and r['shift_2'] < 0 for r in group),
                'direction_criteria_pass_all_four': complete and all(
                    r['provisional_candidate'] and r['distance_to_artificial_target_termini_A'] >= 10 for r in group),
                'models': group,
            })
    result = {'plan': plan, 'pairs': pairs, 'parents_evaluated': len(parents),
              'direction_hypotheses': sum(p['direction_criteria_pass_all_four'] for p in pairs),
              'limitations': ['No experimental binding evidence', 'Mutated sequences have not been independently refolded',
                             'PROPKA uses the intact target crop with artificial termini, not the full receptor',
                             'Free binder uses the bound conformation', 'His-only independent-site approximation',
                             'No mouse binding or measured pH-switching evidence']}
    (report / 'results.json').write_text(json.dumps(result, indent=2))
    (report / 'direction_hypotheses.fasta').write_text(''.join(
        f">{p['id']} evidence=COMPUTATIONAL_DIRECTION_ONLY\n{p['sequence']}\n"
        for p in pairs if p['direction_criteria_pass_all_four']))
    print(json.dumps({'parents': len(parents), 'pairs': len(pairs), 'direction_hypotheses': result['direction_hypotheses']}, indent=2))


if __name__ == '__main__':
    main()
