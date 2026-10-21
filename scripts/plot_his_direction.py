"""Show direction consistency and the frozen heuristic that A41H misses."""
import json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tfr1-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_controls import ROOT


def main():
    out=ROOT/'reports/his_direction';d=json.loads((out/'results.json').read_text())
    a=next(v for v in d['variants'] if v['id']=='apical_short_008_A41H')
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),gridspec_kw={'width_ratios':[1,1.3]})
    ax=axes[0];ax.axhspan(6,7,color='#e7eef0',label='Frozen free-pKa window')
    for i,m in enumerate(a['models']):
        ax.plot([0,1],[m['free_pkas'][0],m['bound_pkas'][0]],'-o',label=f"Seed {m['seed']}",color=['#007f86','#d98627'][i])
        ax.annotate(f"{m['free_pkas'][0]:.2f}",(-.03,m['free_pkas'][0]),ha='right',fontsize=9)
        ax.annotate(f"{m['bound_pkas'][0]:.2f}",(1.03,m['bound_pkas'][0]),ha='left',fontsize=9)
    ax.set_xlim(-.35,1.4);ax.set_ylim(5,7.8);ax.set_xticks([0,1],['Free','Bound']);ax.set_ylabel('Predicted His41 pKa')
    ax.set_title('A41H: lower pKa in both poses',loc='left',fontsize=11,fontweight='bold');ax.legend(fontsize=8,loc='lower left')
    ax=axes[1];labels=[]
    for i,v in enumerate(d['variants']):
        label='A41H (short parent)' if len(v['mutations'])==1 else 'L7H/L82H (long parent)'
        labels.append(label)
        for j,m in enumerate(v['models']):
            value=m['his_only_Kd6_over_Kd7_4'];y=i+(.08 if j else -.08)
            ax.scatter(value,y,color=['#007f86','#d98627'][j],s=55)
            ax.annotate(f'{value:.2f}',(value,y),xytext=(7,0),textcoords='offset points',va='center',fontsize=9)
    ax.axvline(1,color='#666',linestyle='--',linewidth=1);ax.axvline(2,color='#aaa',linestyle=':',linewidth=1)
    ax.set_xscale('log');ax.set_xlim(.2,10);ax.set_ylim(-.6,1.6);ax.set_yticks(range(2),labels,fontsize=9)
    ax.set_xlabel('His-only Kd(6.0) / Kd(7.4) proxy')
    ax.set_title('Above 1 indicates a release-direction proxy',loc='left',fontsize=10,fontweight='bold')
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Computational direction evidence — no variant passes every strict rule',x=.02,ha='left',fontsize=12,fontweight='bold')
    fig.text(.02,.02,'A41H misses the free-pKa window in one pose. Ratios are fixed-conformation proxies, not measured or calibrated affinity predictions.',fontsize=8)
    fig.tight_layout(rect=(0,.07,1,.93));fig.savefig(out/'direction_evidence.png',dpi=180);fig.savefig(out/'direction_evidence.pdf');plt.close(fig)
if __name__=='__main__':main()
