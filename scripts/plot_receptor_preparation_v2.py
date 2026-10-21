"""Static scientific comparison; these are repeated preparations, not new poses."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_controls import ROOT


def main():
    out=ROOT/'reports/receptor_preparation_v2'
    data=json.loads((out/'analysis.json').read_text())
    ids=['apical_short_008_A41H','expansion_apical_short_026_E19H','expansion_apical_long_092_R3H','expansion_apical_short_026_V60H','expansion_apical_short_018_E57H_K60H']
    fig,axes=plt.subplots(1,5,figsize=(14,4.8),sharey=True)
    for ax,name in zip(axes,ids):
        rows=sorted([r for r in data['rows'] if r['id']==name],key=lambda r:r['seed'])
        for x,row in enumerate(rows):
            d=row['direction'];new=d['shifts'];old=d['original_v1_shifts']
            for i,value in enumerate(new):
                offset=(i-(len(new)-1)/2)*.22
                if old:
                    ax.plot([x+offset-.06,x+offset+.06],[old[i],value],color='#9ba6af',lw=1)
                    ax.scatter(x+offset-.06,old[i],facecolors='none',edgecolors='#758390',s=50,zorder=3)
                ax.scatter(x+offset+.06,value,color='#147d92',s=42,zorder=4)
            if not old:ax.text(x,-1.40,'v1 excluded',ha='center',fontsize=8,color='#65717c')
        ax.axhline(0,color='#bbb',lw=.8);ax.axhline(-.5,color='#bd7530',ls='--',lw=1)
        ax.set_xticks([0,1],['1801','1901']);ax.set_xlabel('Prediction seed')
        ax.set_title(name.replace('expansion_apical_','').replace('apical_','').replace('_','\n'),fontsize=10)
        ax.set_xlim(-.45,1.45);ax.set_ylim(-1.55,.9);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Bound − free histidine pKa')
    fig.suptitle('A41H retains the clearest lowering across both preparations',fontsize=15,y=.98)
    fig.text(.5,.075,'Open markers: original fixed-receptor preparation  •  Filled: uniform revised preparation  •  Dashed: −0.5 heuristic',ha='center',fontsize=10)
    fig.text(.5,.025,'Same predicted poses, not independent replication. No variant passes every strict criterion in both poses; full control-panel preparation gate fails.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.12,1,.92))
    fig.savefig(out/'direction_sensitivity.png',dpi=180);fig.savefig(out/'direction_sensitivity.pdf');plt.close(fig)

if __name__=='__main__':main()
