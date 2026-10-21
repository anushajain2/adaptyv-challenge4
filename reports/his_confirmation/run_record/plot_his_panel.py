"""Plot all prediction samples and selected representatives, preserving failures."""
import argparse
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tfr1-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    p=argparse.ArgumentParser();p.add_argument('report');a=p.parse_args()
    from pathlib import Path
    out=Path(a.report);d=json.loads((out/'results.json').read_text())
    names=list(dict.fromkeys(r['id'] for r in d['rows']))
    names=sorted(names,key=lambda n:(not any(r['id']==n and r['arm']=='control' for r in d['rows']),n))
    fig,ax=plt.subplots(figsize=(11,max(6,len(names)*.31+1.8)))
    seeds=d['prediction_seeds'];colors=['#007f86','#e58c2b']
    for i,name in enumerate(names):
        rows=[r for r in d['rows'] if r['id']==name]
        for r in rows:
            offset=(seeds.index(r['seed'])-(len(seeds)-1)/2)*.2
            value=r['metrics']['design_to_target_iptm'];color=colors[seeds.index(r['seed'])%2]
            ax.scatter(value,i+offset,s=20,color=color,alpha=.5,marker='o' if r['target_integrity_pass'] else 'x')
            if r['confidence_selected']:
                passed=r.get('control_pass',r.get('screen',{}).get('pass',False))
                ax.scatter(value,i+offset,s=80,facecolors='none',edgecolors='#16804a' if passed else '#c44441',linewidths=1.5)
    ax.axvline(.5,color='#777',linestyle='--',linewidth=.8)
    ax.set_yticks(range(len(names)),names,fontsize=8);ax.invert_yaxis();ax.set_xlim(0,.8)
    ax.set_xlabel('Interface ipTM (confidence score; not affinity)')
    ax.set_title(f"His {d['stage']}: all samples and frozen selections",loc='left',fontweight='bold')
    ax.grid(axis='x',alpha=.15)
    ax.spines[['top','right']].set_visible(False)
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[i%2],markeredgecolor=colors[i%2],label=f'Seed {seed}') for i,seed in enumerate(seeds)]
    handles += [Line2D([0],[0],marker='o',color='none',markeredgecolor=c,markerfacecolor='none',label=label) for c,label in [('#16804a','Selected: raw checks pass'),('#c44441','Selected: raw checks fail')]]
    ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.42,1.005),fontsize=8)
    fig.text(.02,.015,'Crosses: target-integrity rejection. All earlier results remain preserved. Computational evidence only.',fontsize=8)
    fig.tight_layout(rect=(0,.035,1,1));fig.savefig(out/'prediction_panel.png',dpi=170);fig.savefig(out/'prediction_panel.pdf');plt.close(fig)
if __name__=='__main__':main()
