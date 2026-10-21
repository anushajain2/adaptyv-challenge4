"""Descriptive post-hoc tag-contact confidence and paired core geometry; no new gates."""
import json
import numpy as np
from scipy.spatial import cKDTree
from prepare_controls import ROOT,read_structure
from analyze_controls import pose_metrics
from analyze_tag_comparison import core_structure


def close_pairs(s):
    a=[x for x in s[0]['A'].get_atoms() if x.element!='H'];b=[x for x in s[0]['B'].get_atoms() if x.element!='H']
    tree=cKDTree([x.coord for x in a]);pairs=[]
    for y in b:
        for i in tree.query_ball_point(y.coord,np.nextafter(2.,0.)):
            x=a[i];pairs.append({'target_native_residue':x.parent.id[1]+189,'target_resname':x.parent.resname,'target_atom':x.name,
                'binder_residue':y.parent.id[1],'binder_resname':y.parent.resname,'binder_atom':y.name,'distance_A':float(np.linalg.norm(x.coord-y.coord))})
    return pairs


def main():
    out=ROOT/'reports/tag_comparison';raw=json.loads((out/'results.json').read_text());selected=[r for r in raw['rows'] if r['confidence_selected']]
    lookup={(r['base_id'],r['condition'],r['seed']):r for r in selected};rows=[]
    for r in selected:
        if r['condition']!='tagged':continue
        folder=out/'review_bundle'/f"{r['id']}_seed{r['seed']}"
        s=read_structure(folder/'raw_complex.pdb');core=core_structure(s,r['core_length'])
        other=lookup.get((r['base_id'],'bare',r['seed']));paired=None
        if other:
            bare=read_structure(ROOT/other['raw_scoring_output']/'model.pdb');paired=pose_metrics(core,bare)
        pairs=r['tag']['target_contact_residue_pairs'];confidence=None
        with np.load(folder/'prediction_confidence.npz') as d:
            if pairs:
                a=np.array([t-1 for t,b in pairs]);b=np.array([194+b-1 for t,b in pairs])
                confidence=float((d['pae'][a,b].mean()+d['pae'][b,a].mean())/2)
        rows.append({'id':r['id'],'seed':r['seed'],'kind':r['kind'],'core_failed_checks':r.get('screen',{}).get('failed'),
            'core_clash_pairs':close_pairs(core),'tag_contact_pairs':pairs,'tag_contact_PAE_A':confidence,
            'tag_mean_plddt_0_to_1':r['tag']['mean_plddt'],'core_vs_matched_bare_pose':paired})
    reproducibility=[]
    for name in sorted({r['id'] for r in rows}):
        rs=[r for r in rows if r['id']==name]
        if len(rs)!=2:continue
        a,b=({tuple(p) for p in r['tag_contact_pairs']} for r in rs)
        reproducibility.append({'id':name,'tag_contact_pair_intersection':len(a&b),'tag_contact_pair_union':len(a|b),'jaccard':len(a&b)/len(a|b) if a|b else None})
    result={'models':rows,'across_seed_contacts':reproducibility,
        'scope':'Post-hoc descriptive audit. Neither low-confidence tag contacts nor their absence prove assay compatibility. All predeclared gates and selections remain unchanged.'}
    (out/'tag_contact_audit.json').write_text(json.dumps(result,indent=2))
    for r in rows:print(r['id'],r['seed'],'tagPAE',r['tag_contact_PAE_A'],'corevsbare',r['core_vs_matched_bare_pose']['binder_interface_CA_rmsd_after_target_fit_A'] if r['core_vs_matched_bare_pose'] else None,'clashes',r['core_clash_pairs'])

if __name__=='__main__':main()
