"""Selected-model core confidence versus length-confounded whole-chain confidence."""
import json,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tfr1-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_controls import ROOT

def main():
    out=ROOT/'reports/tag_comparison';r=json.loads((out/'results.json').read_text());m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())
    names=list(dict.fromkeys(c['base_id'] for c in m['cases']))
    labels=['008 parent','008 A41H','026 parent','026 E19H','Positive','Scrambled','Unrelated VHH']
    fig,axes=plt.subplots(1,2,figsize=(12,5.3))
    for condition,color,offset in [('bare','#176a8a',-.12),('tagged','#d67a25',.12)]:
        for i,name in enumerate(names):
            rows=[x for x in r['rows'] if x['confidence_selected'] and x['base_id']==name and x['condition']==condition]
            for ax,key in zip(axes,['core_mean_cross_PAE_A','iptm']):
                values=[x[key] if key in x else x['metrics'][key] for x in rows]
                if not values:continue
                ax.plot([i+offset]*len(values),values,'o-',color=color,label=condition if i==0 else None,ms=6,alpha=.85)
    axes[0].axhline(m['core_confidence']['max_mean_cross_PAE_A'],ls='--',color='#777',lw=1)
    axes[0].set_ylabel('Target–binder CORE mean PAE (Å)');axes[0].set_title('Core confidence used for selection')
    axes[1].set_ylabel('Whole-complex ipTM');axes[1].set_title('Whole-complex score retained for comparison')
    for ax in axes:
        ax.set_xticks(range(len(names)),labels,rotation=35,ha='right');ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False)
    fig.suptitle('Matched bare and 64-residue-tagged constructs',fontsize=14)
    fig.text(.5,.025,'Two fresh seeds, five samples each; dots show selected models. Core confidence alone does not establish binding or tag compatibility.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.94));fig.savefig(out/'core_confidence.png',dpi=180);fig.savefig(out/'core_confidence.pdf');plt.close(fig)

if __name__=='__main__':main()
