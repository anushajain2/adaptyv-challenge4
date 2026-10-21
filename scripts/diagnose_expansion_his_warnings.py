"""Read-only geometry audit of PROPKA warning groups; never repair or rescore them."""
import json,re
import numpy as np
from scipy.spatial import cKDTree
from openmm import app,unit
from prepare_controls import ROOT,read_structure

EXPECTED={'ARG':{'N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'},'GLU':{'N','CA','C','O','CB','CG','CD','OE1','OE2'},'ASN':{'N','CA','C','O','CB','CG','OD1','ND2'},'SER':{'N','CA','C','O','CB','OG'}}


def main():
    out=ROOT/'reports/expansion_his_direction';d=json.loads((out/'results.json').read_text());records=[]
    for v in d['variants']:
        for m in v['models']:
            folder=out/'models'/f"{v['id']}_seed{m['seed']}";p=folder/'complex.pdb'
            if not p.exists():continue
            warnings=json.loads((folder/'complex.propka_warnings.json').read_text())
            groups=sorted(set((match[0],int(match[1])) for w in warnings for match in re.findall(r'Missing atoms or failed protonation for (\w+) (\d+) A',w)))
            pdb=app.PDBFile(str(p));atoms=list(pdb.topology.atoms());xyz=np.asarray(pdb.positions.value_in_unit(unit.angstrom));bonds={tuple(sorted((a.index,b.index))) for a,b in pdb.topology.bonds()}
            target=[a for a in atoms if a.residue.chain.id=='A' and a.element.symbol!='H'];pairs=cKDTree(xyz[[a.index for a in target]]).query_pairs(1.8)
            suspicious=[]
            for i,j in pairs:
                a,b=target[i],target[j]
                if tuple(sorted((a.index,b.index))) in bonds:continue
                suspicious.append({'a':f'A:{a.residue.id}:{a.residue.name}:{a.name}','b':f'A:{b.residue.id}:{b.residue.name}:{b.name}','distance_A':float(np.linalg.norm(xyz[a.index]-xyz[b.index])),'warning_residue_involved':any(int(x.residue.id)==pos for x in (a,b) for _,pos in groups)})
            s=read_structure(p);details=[]
            for name,pos in groups:
                r=s[0]['A'][pos];missing=sorted(EXPECTED[name]-set(r.child_dict));ring=[a for pp in m['positions'] for a in s[0]['B'][pp] if a.name in {'CG','ND1','CD2','CE1','NE2'}]
                dist=float(np.linalg.norm(np.array([a.coord for a in r])[:,None]-np.array([a.coord for a in ring])[None],axis=2).min())
                details.append({'crop_position':pos,'native_position':pos+189,'resname':name,'missing_expected_heavy_atoms':missing,'nearest_introduced_His_ring_A':dist})
            records.append({'id':v['id'],'seed':m['seed'],'warning_groups':details,'intra_target_nonbonded_pairs_lt1_8A':sorted(suspicious,key=lambda x:x['distance_A'])})
    result={'models':records,'interpretation':'All findings are descriptive. Bond exclusions use standard PDB/OpenMM topology. Close nonbonded receptor geometry may disrupt distance-based protonation/bond assignment; this is not proof of each warning cause. No coordinates were changed, no pKa rerun, and no failed score was promoted.'}
    (out/'protonation_warning_diagnosis.json').write_text(json.dumps(result,indent=2))
    for r in records:
        if r['warning_groups']:print(r['id'],r['seed'],r['warning_groups'],[x for x in r['intra_target_nonbonded_pairs_lt1_8A'] if x['warning_residue_involved']])

if __name__=='__main__':main()
