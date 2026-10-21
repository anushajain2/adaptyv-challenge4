"""Plot all computed directions, with strict-rule and preparation failures explicit."""
import json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tfr1-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from prepare_controls import ROOT


def main():
    out=ROOT/'reports/expansion_his_direction';d=json.loads((out/'results.json').read_text());variants=d['variants']
    if not variants:return
    fig,axes=plt.subplots(1,2,figsize=(13,max(4,len(variants)*.65+1.7)),gridspec_kw={'width_ratios':[1.15,1]})
    labels=[];colors=['#007f86','#e58c2b']
    for i,v in enumerate(variants):
        labels.append(v['id'].replace('expansion_apical_','').replace('apical_short_008_','prior_short_008_'))
        for j,m in enumerate(v['models']):
            y=i+(j-.5)*.22
            if m['status']!='ok':
                axes[0].text(0,y,'Not scored',fontsize=8,va='center');continue
            axes[0].scatter(m['shifts'],[y]*len(m['shifts']),s=40,color=colors[j])
            x=m['his_only_Kd6_over_Kd7_4'];axes[1].scatter(x,y,color=colors[j],s=45)
            axes[1].scatter(x,y,facecolors='none',edgecolors='#16804a' if m['direction_pass'] else '#c44441',s=105,linewidths=1.2)
            axes[1].annotate(f'{x:.2f}',(x,y),xytext=(6,3),textcoords='offset points',fontsize=8)
    axes[0].axvline(0,color='#aaa',linewidth=.8);axes[0].axvline(-.5,color='#666',linestyle='--',linewidth=1)
    axes[0].set_xlabel('Bound minus free His pKa (each introduced site)');axes[0].set_yticks(range(len(labels)),labels,fontsize=8)
    axes[1].set_xscale('log');axes[1].axvline(1,color='#aaa',linewidth=.8);axes[1].axvline(2,color='#666',linestyle='--',linewidth=1)
    axes[1].set_xlabel('His-only Kd(6.0) / Kd(7.4) proxy');axes[1].set_yticks(range(len(labels)),['']*len(labels));axes[1].margins(x=.25);axes[1].set_xticks([.5,1,2,5],['0.5','1','2','5']);axes[1].minorticks_off()
    for ax in axes:ax.set_ylim(len(labels)-.5,-.65);ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.12)
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=c,markeredgecolor=c,label=f'Seed {seed}') for c,seed in zip(colors,[1801,1901])]
    handles +=[Line2D([0],[0],marker='o',color='none',markerfacecolor='none',markeredgecolor=c,label=l) for c,l in [('#16804a','Strict model rule passes'),('#c44441','Strict model rule fails')]]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.62,.93),ncol=2,fontsize=8)
    fig.suptitle('New-parent histidine direction: conditional computational evidence',x=.02,ha='left',fontweight='bold',fontsize=12)
    fig.text(.02,.025,'Ratios are fixed-conformation, independent-His proxies, not measured or calibrated affinity predictions. All failures remain recorded.',fontsize=8)
    fig.tight_layout(rect=(0,.06,1,.86));fig.savefig(out/'direction_evidence.png',dpi=180);fig.savefig(out/'direction_evidence.pdf');plt.close(fig)

if __name__=='__main__':main()
