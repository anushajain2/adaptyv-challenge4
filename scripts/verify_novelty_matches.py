"""Check selected structural matches against original PDB coordinates locally."""
import itertools
import json

import numpy as np
from Bio.PDB import MMCIFParser
from tmtools.io import get_residue_data

from novelty_search import ROOT, OUT, save, sha, now
from analyze_novelty import structural_alignment
from prepare_controls import read_structure
from build_submission_draft import alignment


def main():
    plan=json.loads((OUT/'plan.json').read_text())
    targets={'expansion_apical_short_015':('2OO2','A'),
             'expansion_apical_long_092':('3SSQ','A'),
             'expansion_apical_short_026':('8A53','C')}
    rows=[]
    for c in plan['structures']:
        if c['candidate'] not in targets:continue
        pdb,chain=targets[c['candidate']];path=OUT/'independent_pdb'/f'{pdb}.cif'
        s=MMCIFParser(QUIET=True).get_structure(pdb,str(path))
        x,sx=get_residue_data(read_structure(ROOT/c['query_file'])[0]['B'])
        y,sy=get_residue_data(s[0][chain])
        entry=json.loads((OUT/'independent_pdb'/f'{pdb}.entry.json').read_text())
        rows.append({'query':c['id'],'candidate':c['candidate'],'pdb':pdb,'auth_chain':chain,
                     'pdb_title':entry['struct']['title'],'initial_release_date':entry['rcsb_accession_info']['initial_release_date'],
                     'coordinate_source':'https://files.rcsb.org/download/'+pdb+'.cif','coordinate_sha256':sha(path),
                     'query_length':len(sx),'modeled_subject_length':len(sy),'modeled_subject_sequence':sy,
                     'global_sequence_identity':alignment(sx,sy)[1],
                     **structural_alignment(x,y,sx,sy)})
    save(OUT/'independent_pdb_verification.json',{'rows':rows,'scope':'Original asymmetric-unit coordinates, named author chains. This verifies selected descriptive matches, not the entire database or biological function.'})
    for r in rows:print(r['query'],r['pdb'],'qTM/tTM',round(r['query_normalized_TM'],3),round(r['target_normalized_TM'],3),
                        'coverage',round(r['query_coverage'],3),round(r['target_coverage'],3),'RMSD',round(r['RMSD_A'],3),
                        'global identity',round(r['global_sequence_identity'],3))


if __name__=='__main__':main()
