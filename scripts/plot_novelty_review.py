"""Scientific figure: original-coordinate verification of the 2OO2 match."""
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(tempfile.gettempdir())/'tfr1-novelty-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser
from tmtools import tm_align
from tmtools.io import get_residue_data
from novelty_search import ROOT, OUT
from prepare_controls import read_structure


def main():
    verification=json.loads((OUT/'independent_pdb_verification.json').read_text())
    plan=json.loads((OUT/'plan.json').read_text())
    rows=[r for r in verification['rows'] if r['pdb']=='2OO2']
    target=MMCIFParser(QUIET=True).get_structure('2OO2',str(OUT/'independent_pdb/2OO2.cif'))
    y,sy=get_residue_data(target[0]['A'])
    fig=plt.figure(figsize=(10,5.2),facecolor='white')
    for i,row in enumerate(rows):
        query=next(c for c in plan['structures'] if c['id']==row['query'])
        x,sx=get_residue_data(read_structure(ROOT/query['query_file'])[0]['B'])
        r=tm_align(x.astype(float),y.astype(float),sx,sy)
        transformed=x@r.u.T+r.t
        # Verify plotting transform with the alignment's matched pairs.
        qi=ti=0;d=[]
        for a,b,marker in zip(r.seqxA,r.seqyA,r.seqM):
            if a!='-' and b!='-':
                distance=np.linalg.norm(transformed[qi]-y[ti]);d.append(distance)
                # The returned transform optimizes TM score, whereas reported
                # RMSD need not be exactly the RMSD under that transform.
                if marker==':':assert distance<5.001
                if marker=='.':assert distance>=4.999
            qi+=a!='-';ti+=b!='-'
        ax=fig.add_subplot(1,2,i+1,projection='3d')
        ax.plot(*y.T,color='#a3acb8',lw=2.2,label='Known protein AF1782 (2OO2 A)')
        ax.plot(*transformed.T,color='#086b89',lw=2.6,label='short_015 design')
        ax.scatter(*transformed[0],color='#086b89',s=22)
        allxyz=np.vstack([y,transformed]);center=(allxyz.max(0)+allxyz.min(0))/2
        radius=np.ptp(allxyz,axis=0).max()/2+1
        ax.set_xlim(center[0]-radius,center[0]+radius);ax.set_ylim(center[1]-radius,center[1]+radius);ax.set_zlim(center[2]-radius,center[2]+radius)
        ax.set_box_aspect((1,1,1));ax.view_init(elev=16,azim=-55);ax.set_axis_off()
        seed=row['query'].split('seed')[-1]
        ax.set_title(f'Seed {seed} · TM-align RMSD {r.rmsd:.2f} Å\nTM (design / known chain): {r.tm_norm_chain1:.3f} / {r.tm_norm_chain2:.3f}',fontsize=11)
    fig.suptitle('short_015 shares a known backbone arrangement',fontsize=16,fontweight='bold',y=.98)
    handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.08),ncol=2,frameon=False,fontsize=10)
    fig.text(.5,.025,'All 61 design residues align; 61/76 modeled AF1782 residues. Global sequence identity: 27.6%.\nCα traces, binder-only TM-align. Structural similarity does not establish copied sequence or binding function.',ha='center',fontsize=9,color='#414b59')
    fig.subplots_adjust(top=.80,bottom=.15,wspace=.04,left=.02,right=.98)
    fig.savefig(OUT/'short015_known_fold.png',dpi=180)
    fig.savefig(OUT/'short015_known_fold.pdf')
    plt.close(fig)


if __name__=='__main__':main()
