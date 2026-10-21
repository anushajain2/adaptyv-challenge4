"""Show every confirmation sample and the prospective selection filter."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main():
    out = ROOT / 'reports/confirmation'
    data = json.loads((out / 'results.json').read_text())
    cases = ['apical_short_008', 'apical_long_015', '2ds25_5_6wrw', 'scrambled_6wrw', 'caplacizumab_6wrw']
    labels = ['Short lead 008', 'Long lead 015', '2DS25.5 positive', 'Scrambled control', 'Unrelated VHH']
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    for seed, color in [(607, '#2563eb'), (709, '#7c3aed')]:
        rows = [r for r in data['rows'] if r['seed'] == seed]
        for ax, value in [(axes[0], lambda r: r['metrics']['design_to_target_iptm']),
                          (axes[1], lambda r: r['target_CA_rmsd_to_input_A'])]:
            xx = [cases.index(r['id']) + (-.16 if seed == 607 else .16) + .045*(r['sample']-2) for r in rows]
            ax.scatter(xx, [value(r) for r in rows], s=32, color=color, alpha=.75, label=f'Seed {seed}')
            for x, row in zip(xx, rows):
                if row['confidence_selected']:
                    ax.scatter(x, value(row), s=115, facecolors='none', edgecolors='black', linewidth=1.3)
                if not row['target_integrity_pass']:
                    ax.scatter(x, value(row), s=95, marker='x', color='#dc2626', linewidth=2)
    axes[0].axhline(.5, ls='--', color='#64748b', lw=1)
    axes[1].axhline(2, ls='--', color='#64748b', lw=1)
    axes[0].set(ylabel='Interface ipTM', title='Confidence of all 50 samples')
    axes[1].set(ylabel='Target Cα RMSD to supplied structure (Å)', title='Target-conformation filter before selection')
    for ax in axes:
        ax.set_xticks(range(len(cases)), labels, rotation=18, ha='right')
        ax.grid(axis='y', alpha=.15)
    axes[0].legend(loc='upper right', fontsize=8)
    fig.suptitle('TfR1 confirmation · controls pass; one parent advances', fontsize=14)
    fig.supxlabel('Black rings: selected after filtering. Red cross: rejected target conformation. Predictive scores do not establish binding.', fontsize=9)
    fig.savefig(out / 'confirmation.png', dpi=170)
    fig.savefig(out / 'confirmation.pdf')
    plt.close(fig)


if __name__ == '__main__':
    main()
