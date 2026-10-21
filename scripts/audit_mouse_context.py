"""Audit mouse AlphaFold correspondence; do not infer cross-species binding."""
import copy,csv,hashlib,json
import numpy as np
from Bio import SeqIO
from Bio.PDB import Superimposer
from prepare_controls import ROOT,read_structure,sequence,save_pdb,make_complex,residues


def main():
    src=ROOT/'data/mouse_context_audit';out=ROOT/'reports/mouse_context_audit';out.mkdir(exist_ok=False)
    mouse=read_structure(src/'mouse_Q62351.cif');human=read_structure(ROOT/'data/apical_expansion_fixed/target_6wrw.cif')
    expected=str(next(SeqIO.parse(ROOT/'data/raw/Q62351.fasta','fasta')).seq)
    assert sequence(mouse[0]['A'])==expected and len(expected)==763
    mapping=list(csv.DictReader((ROOT/'data/prepared/conservation_map.csv').open()))
    anchors=[(int(r['human_uniprot']),int(r['mouse_uniprot'])) for r in mapping if r['mouse_uniprot'] and r['conserved']=='True' and 190<=int(r['human_uniprot'])<=383 and not 198<=int(r['human_uniprot'])<=218 and mouse[0]['A'][int(r['mouse_uniprot'])]['CA'].bfactor>=70]
    fit=Superimposer();fit.set_atoms([human[0]['A'][h-189]['CA'] for h,m in anchors],[mouse[0]['A'][m]['CA'] for h,m in anchors]);fit.apply(list(mouse.get_atoms()))
    loop=[]
    for h in range(204,217):
        distances=sorted((float(np.linalg.norm(human[0]['A'][h-189]['CA'].coord-mouse[0]['A'][m]['CA'].coord)),m) for m in range(206,220))
        row=next(r for r in mapping if int(r['human_uniprot'])==h)
        loop.append({'human_native':h,'human_resname':human[0]['A'][h-189].resname,'sequence_alignment_mouse':row['mouse_uniprot'],'published_correspondence_mouse':row['published_mouse_position'],'alignment_caution':row['alignment_caution'],'nearest_mouse_CA_candidates':[{'mouse_native':m,'resname':mouse[0]['A'][m].resname,'distance_A':d,'pLDDT':mouse[0]['A'][m]['CA'].bfactor} for d,m in distances[:3]]})
    # Archive the source model in the human frame, retaining native mouse numbering.
    save_pdb(mouse,out/'mouse_model_in_human_frame.pdb')
    result={'source':json.loads((src/'source.json').read_text()),'full_mouse_sequence_matches_Q62351':True,'anchors':anchors,'anchor_CA_rmsd_A':float(fit.rms),'excluded_loop_human':[198,218],'loop_audit':loop,'mouse_crop_candidate_native':[192,385],'mouse_crop_length':194,'mouse_crop_sequence':expected[191:385],'interpretation':'This is a predicted monomer and a structural correspondence diagnostic. Nearest C-alpha distance is not a definitive evolutionary alignment or a bound receptor conformation. Do not replace the sequence map or infer mouse binding from this audit.'}
    (out/'results.json').write_text(json.dumps(result,indent=2))
    print('Anchors',len(anchors),'RMSD',fit.rms)
    for r in loop:
        if r['alignment_caution']:print(r['human_native'],r['sequence_alignment_mouse'],r['published_correspondence_mouse'],r['nearest_mouse_CA_candidates'])
if __name__=='__main__':main()
