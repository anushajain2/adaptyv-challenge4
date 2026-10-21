"""Show every advisory species/conformation score without mouse-binding labels."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from prepare_controls import ROOT


def main():
    out=ROOT/'reports/species_exploration';d=json.loads((out/'results.json').read_text())
    names=list(dict.fromkeys(r['base_id'] for r in d['rows']))
    conditions=[('hbound','Human crystal','#00868b'),('haf','Human AlphaFold','#d78327'),('maf','Mouse AlphaFold','#7a62a5')]
    fig,ax=plt.subplots(figsize=(12,6.4))
    for j,(condition,label,color) in enumerate(conditions):
        for i,name in enumerate(names):
            rows=[r for r in d['rows'] if r['condition']==condition and r['base_id']==name]
            for k,r in enumerate(rows):
                y=i+(j-1)*.20
                x=r['metrics']['design_to_target_iptm']
                ax.scatter(x,y,color=color,s=30,alpha=.4,marker='o' if r['target_integrity_pass'] else 'x',label=label if i==0 and k==0 else None)
                if r['confidence_selected']:ax.scatter(x,y,facecolors='none',edgecolors=color,s=100,linewidths=1.4)
    ax.set_yticks(range(len(names)),names);ax.invert_yaxis();ax.set_xlim(0,.8)
    ax.set_xlabel('Interface ipTM (confidence score; not affinity)')
    ax.set_title('Species comparison is inconclusive: human AlphaFold positive control fails',loc='left',weight='bold',pad=18)
    ax.grid(axis='x',alpha=.15);ax.legend(loc='lower right',frameon=False)
    ax.spines[['top','right']].set_visible(False)
    fig.text(.02,.02,'Ring: selected representative. Cross: target-integrity rejection. One seed, three samples per condition; no validated mouse positive.',fontsize=9)
    fig.tight_layout(rect=[0,.045,1,1]);fig.savefig(out/'comparison.png',dpi=180);fig.savefig(out/'comparison.pdf');plt.close(fig)

if __name__=='__main__':main()
