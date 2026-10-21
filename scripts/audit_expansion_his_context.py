"""Annotate successful direction models with terminal and interface context; no new gate."""
import csv,json
import numpy as np
from prepare_controls import ROOT,read_structure


def distance(a,b):
    return float(np.linalg.norm(np.array([x.coord for x in a])[:,None]-np.array([x.coord for x in b])[None],axis=2).min())


def main():
    out=ROOT/'reports/expansion_his_direction';d=json.loads((out/'results.json').read_text())
    conservation={int(r['native_position']):r for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    rows=[]
    for variant in d['variants']:
        for model in variant['models']:
            path=out/'models'/f"{variant['id']}_seed{model['seed']}"/'complex.pdb'
            if not path.exists():continue
            s=read_structure(path);binder=list(s[0]['B']);sites=[]
            for pos in model['positions']:
                ring=[a for a in s[0]['B'][pos] if a.name in {'CG','ND1','CD2','CE1','NE2'}];assert len(ring)==5
                terminal={'N':distance(ring,[binder[0]['N']]),'C_carboxyl':distance(ring,[a for a in binder[-1] if a.name in {'C','O','OXT'}])}
                contacts=[]
                for r in s[0]['A']:
                    delta=distance(ring,[a for a in r if a.element!='H'])
                    if delta<5:
                        native=r.id[1]+189;contacts.append({'target_native_position':native,'resname':r.resname,'ring_min_distance_A':delta,'human_mouse_conserved':conservation[native]['conserved']=='True'})
                sites.append({'position':pos,'ring_to_binder_termini_A':terminal,'target_contacts_lt5A':contacts,'caution':'Native termini or tag context may affect the calculation; no assay tag modeled.' if min(terminal.values())<10 else 'No ring atom within 10 A of binder termini; tag compatibility remains untested.'})
            rows.append({'id':variant['id'],'seed':model['seed'],'pka_status':model['status'],'sites':sites})
    (out/'context_audit.json').write_text(json.dumps({'models':rows,'scope':'Descriptive audit, not an extra post-hoc pass/fail rule. No energetic, tag-compatibility or mouse-binding claim.'},indent=2))
    print('Audited context for',len(rows),'models')

if __name__=='__main__':main()
