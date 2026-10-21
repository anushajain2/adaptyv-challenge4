"""Check every prepared core, including parents and controls; retain raw exclusions."""
import csv,json
from prepare_controls import ROOT,read_structure,geometry
from analyze_controls import pose_metrics
from analyze_generation_comparison import design_pose_metrics,contacts_summary
from analyze_tag_comparison import core_structure,core_screen
from receptor_preparation_v2 import sha

def main():
    out=ROOT/'reports/tag_comparison';raw=json.loads((out/'results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text())
    m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text());cases={c['id']:c for c in m['cases']};arms={a['id']:a for a in m['arms']}
    selected={(r['id'],r['seed']):r for r in raw['rows'] if r['confidence_selected']}
    assert len(prep['rows'])==len(selected) and {(r['id'],r['seed']) for r in prep['rows']}==set(selected)
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    full=read_structure(ROOT/'data/raw/11QC.cif');rows=[]
    for p in prep['rows']:
        r=selected[(p['id'],p['seed'])];c=cases[p['id']]
        item={'id':r['id'],'condition':r['condition'],'kind':r['kind'],'seed':r['seed'],'sample':r['sample'],
            'raw_pass':r.get('screen',{}).get('pass',r.get('control_pass')),'prepared_eligible':False,'preparation_status':p['status']}
        if p['status']=='ok':
            path=ROOT/p['dest']/'prepared.pdb';assert sha(path)==p['prepared_sha256']
            s=core_structure(read_structure(path),r['core_length']);ref=read_structure(ROOT/c['reference'])
            if c['kind']=='novel_design':
                scored={**r,'pose':design_pose_metrics(s,ref),'geometry':geometry(s,{i:i+189 for i in range(1,195)},full),'contacts':contacts_summary(s,arms[c['arm']],conserved)}
                screen=core_screen(scored,m)
            elif c['kind']=='positive':
                pose=pose_metrics(s,ref);item['pose']=pose
                screen={'pass':r['core_mean_cross_PAE_A']<m['core_confidence']['max_mean_cross_PAE_A'] and pose['binder_CA_fold_rmsd_A']<=2 and pose['binder_interface_CA_rmsd_after_target_fit_A'] is not None and pose['binder_interface_CA_rmsd_after_target_fit_A']<=2 and pose['native_contact_recall']>=.5}
            else:screen={'pass':r['control_pass'],'note':'Original core PAE; not recalculated by minimization'}
            item.update(prepared_screen=screen,movement_pass=p['displacement']['pass'],internal_geometry_pass=p['prepared_internal_geometry']['pass'],protonation_pass=all(v['status']=='ok' for v in p['propka'].values()))
            item['prepared_eligible']=bool(item['raw_pass'] and screen['pass'] and item['movement_pass'] and item['internal_geometry_pass'] and item['protonation_pass'])
        rows.append(item)
    controls={condition:len(rs:= [r for r in rows if r['condition']==condition and r['kind']!='novel_design'])==6 and all(r['prepared_eligible'] for r in rs) for condition in ('bare','tagged')}
    result={'rows':rows,'prepared_controls_pass':controls,'source_preparation_sha256':sha(out/'preparation/results.json'),
        'scope':'Prepared evidence is a separate diagnostic. Raw exclusions remain failures; full-tag movement is reported separately from core guards.'}
    (out/'prepared_core_audit.json').write_text(json.dumps(result,indent=2));print('Prepared controls',controls)

if __name__=='__main__':main()
