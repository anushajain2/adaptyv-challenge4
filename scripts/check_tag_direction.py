"""Retain full suffix chemistry for tagged free/bound pKa sensitivity."""
import csv,json
from prepare_controls import ROOT,read_structure,geometry
from analyze_tag_comparison import core_structure,core_screen
from analyze_generation_comparison import design_pose_metrics,contacts_summary
from analyze_controls import contact_pairs
from direction_check import histidines
from analyze_receptor_preparation_v2 import direction_metrics
from summarize_direction import terminal_distance
from receptor_preparation_v2 import sha


def main():
    root=ROOT/'reports/tag_comparison';raw=json.loads((root/'results.json').read_text())
    m=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text());cases={c['id']:c for c in m['cases']};arms={a['id']:a for a in m['arms']}
    prep=json.loads((root/'preparation/results.json').read_text())
    prepared_audit=json.loads((root/'prepared_core_audit.json').read_text())
    assert prep['source_results_sha256']==sha(root/'results.json')
    selected={(r['id'],r['seed'],r['sample']) for r in raw['rows'] if r['confidence_selected']}
    assert selected=={(r['id'],r['seed'],r['sample']) for r in prep['rows']} and len(selected)==len(prep['rows'])
    prepared={(r['id'],r['seed']):r for r in prep['rows']}
    full=read_structure(ROOT/'data/raw/11QC.cif')
    conserved={int(r['native_position']):r['conserved']=='True' for r in csv.DictReader((ROOT/'data/broad_campaign/spatial_map.csv').open())}
    variants=[]
    for candidate in raw['candidate_summaries']:
        c=cases[candidate['id']]
        if not c.get('mutations'):continue
        positions=[int(x[1:-1]) for x in c['mutations']];models=[]
        for prediction in candidate['selected_models']:
            p=prepared[(c['id'],prediction['seed'])]
            row={'id':c['id'],'seed':prediction['seed'],'status':'not_scored','direction_pass':False,'tag_independent_direction_pass':False}
            try:
                assert candidate['core_pass_both_seeds'],'Raw core screen or controls fail in a selected seed'
                assert p['status']=='ok' and p['displacement']['pass'] and p['prepared_internal_geometry']['pass'],'Prepared core movement/internal geometry failed'
                assert all(v['status']=='ok' for v in p['propka'].values()),'Full-construct PROPKA warning/failure'
                path=ROOT/p['dest']/'prepared.pdb';assert sha(path)==p['prepared_sha256']
                s=read_structure(path);core=core_structure(s,c['core_length']);reference=read_structure(ROOT/c['reference'])
                scored={**prediction,'pose':design_pose_metrics(core,reference),'geometry':geometry(core,{i:i+189 for i in range(1,195)},full),'contacts':contacts_summary(core,arms[c['arm']],conserved)}
                row['prepared_core_screen']=core_screen(scored,m)
                assert row['prepared_core_screen']['pass'],'Prepared core fails unchanged structural checks'
                folder=path.parent
                free=histidines(json.loads((folder/'free_binder.pka.json').read_text()),'B',positions)
                bound=histidines(json.loads(path.with_suffix('.pka.json').read_text()),'B',positions)
                row.update(status='ok',**direction_metrics(free,bound,terminal_distance(path,positions)),source_sha256=sha(path))
                tag_pairs=[list(pair) for pair in sorted(contact_pairs(s)) if pair[1]>c['core_length']]
                row['prepared_tag_target_contacts']=tag_pairs
                g=geometry(s,{i:i+189 for i in range(1,195)},full)
                row['prepared_full_construct_geometry']=g
                tag_clear=not tag_pairs and all(g[k]==0 for k in g if 'clash_pairs_lt2A' in k)
                row['tag_independent_direction_pass']=row['direction_pass'] and prediction['tag']['pass'] and tag_clear
            except Exception as e:row.update(status='excluded',reason=f'{type(e).__name__}: {e}')
            models.append(row)
        parent_rows=[r for r in prepared_audit['rows'] if r['id']==c.get('parent_id','')+'_'+c['condition']]
        variants.append({'id':c['id'],'base_id':c['base_id'],'condition':c['condition'],'parent_id':c.get('parent_id'),
            'raw_core_pass_both_seeds':candidate['core_pass_both_seeds'],'raw_tag_independent_pass_both_seeds':candidate['tag_independent_pass_both_seeds'],
            'prepared_controls_pass':prepared_audit['prepared_controls_pass'][c['condition']],
            'matched_parent_prepared_pass_both_seeds':len(parent_rows)==2 and all(r['prepared_eligible'] for r in parent_rows),
            'models':models,'direction_pass_both_seeds':len(models)==2 and all(r['direction_pass'] for r in models),
            'tag_independent_direction_pass_both_seeds':len(models)==2 and all(r['tag_independent_direction_pass'] for r in models)})
    for v in variants:
        v['fully_qualified_direction_pass_both_seeds']=v['prepared_controls_pass'] and v['tag_independent_direction_pass_both_seeds']
    result={'source_results_sha256':sha(root/'results.json'),'preparation_results_sha256':sha(root/'preparation/results.json'),
        'variants':variants,'scope':'Full binder+suffix free/bound chemistry; only designed mutation sites contribute to reported His-only ratio. No tag histidine counted as a designed mechanism. No experimental binding or capture-surface model.'}
    (root/'direction.json').write_text(json.dumps(result,indent=2))
    for v in variants:print(v['id'],[(r['seed'],r['status'],r.get('shifts'),r.get('reason')) for r in v['models']])

if __name__=='__main__':main()
