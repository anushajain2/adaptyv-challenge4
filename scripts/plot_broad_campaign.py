"""Plot the complete one-sample campaign without implying experimental validation."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def plot(folder):
    d=json.loads((folder/'results.json').read_text());rows=d['ranked_candidates'];arms=list(d['arms'])
    colors=dict(zip(arms,['#2563eb','#16a34a','#7c3aed','#d97706','#dc2626']))
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
    for arm in arms:
        rs=[r for r in rows if r['arm']==arm]
        axes[0].scatter([r['metrics']['design_to_target_iptm'] for r in rs],
                        [r['scores']['contact_pae_both_directions_A'] if r['scores']['contact_pae_both_directions_A'] is not None else np.nan for r in rs],
                        color=colors[arm],s=24,alpha=.65,label=arm.replace('_',' '))
        eligible=[r for r in rs if r['geometry_eligible']]
        axes[0].scatter([r['metrics']['design_to_target_iptm'] for r in eligible],
                        [r['scores']['contact_pae_both_directions_A'] for r in eligible],
                        facecolors='none',edgecolors='black',s=65,linewidth=.8)
    for r in d['controls']:
        axes[0].scatter(r['metrics']['design_to_target_iptm'],r['scores']['contact_pae_both_directions_A'],
                        marker='*' if r['kind']=='positive' else 'x',s=110,color='black')
    axes[0].axvline(.5,color='#64748b',ls='--',lw=1)
    axes[0].set(xlabel='Interface ipTM',ylabel='Contact-pair PAE (Å)',title='One refolding sample per sequence')
    axes[0].legend(fontsize=8,loc='best')
    labels=[a.replace('_','\n') for a in arms];x=np.arange(len(arms))
    axes[1].bar(x-.18,[d['arms'][a]['geometry_eligible'] for a in arms],.36,label='Geometry eligible',color='#0284c7')
    axes[1].bar(x+.18,[d['arms'][a]['single_sample_full_screen'] for a in arms],.36,label='Also ipTM ≥0.50*',color='#0f766e')
    axes[1].set(xticks=x,xticklabels=labels,ylabel='Design count (40 per arm)',title='Screening outcomes')
    axes[1].legend(fontsize=8)
    for arm in arms:
        rs=[r for r in rows if r['arm']==arm]
        axes[2].scatter([r['pose']['binder_CA_fold_rmsd_A'] for r in rs],
                        [r['pose']['binder_interface_CA_rmsd_after_target_fit_A'] for r in rs],s=24,alpha=.65,color=colors[arm])
    axes[2].axvline(2,color='#64748b',ls='--',lw=1);axes[2].axhline(2,color='#64748b',ls='--',lw=1)
    axes[2].set(xlabel='Binder self-fit Cα RMSD (Å)',ylabel='Interface Cα RMSD in receptor frame (Å)',title='Generated-reference consistency')
    fig.suptitle('TfR1 · 200-design broad screen · predictions, not validated binding',fontsize=15)
    fig.supxlabel('Circle outlines: geometry eligible. Stars: published positive; crosses: presumed negatives. *Nanobody confidence threshold is descriptive.',fontsize=9)
    fig.savefig(folder/'screening.png',dpi=170);fig.savefig(folder/'screening.pdf');plt.close(fig)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path,default=ROOT/'reports/broad_campaign',nargs='?');plot(p.parse_args().folder)
