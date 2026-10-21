"""Apply frozen raw screens, then select distinct parents for one confirmation."""
import argparse
import json
from tmtools.io import get_residue_data
from analyze_apical_expansion import analyze
from analyze_generation_comparison import aligned_identity
from analyze_novelty import structural_alignment
from prepare_controls import read_structure
from novelty_search import ROOT, sha, save, now


def binder(path):
    return get_residue_data(read_structure(path)[0]['B'])


def reciprocal(a,b):
    x,sx=a;y,sy=b
    r=structural_alignment(x,y,sx,sy)
    return {'query_TM':r['query_normalized_TM'],'target_TM':r['target_normalized_TM'],
            'reciprocal_TM':min(r['query_normalized_TM'],r['target_normalized_TM']),
            'query_coverage':r['query_coverage'],'target_coverage':r['target_coverage']}


def screen(run):
    out=ROOT/'reports/diversity_generation';r=analyze(run,out)
    catalog=json.loads((run/'inputs/prior_catalog.json').read_text())
    oldseq={c['sequence'] for c in catalog['prior_designs']};draft=catalog['current_draft']
    oldmodels=[(c['name'],q,binder(ROOT/q['path'])) for c in draft for q in c['latest_selected_bare_models']]
    newmodels={};comparisons=[]
    for row in r['ranked_candidates']:
        row['exact_prior_sequence']=row['sequence'] in oldseq
        row['similar_to_current_draft']=any(aligned_identity(row['sequence'],c['sequence'])>=.8 for c in draft)
        row['structurally_redundant_to_current_draft']=None
        row['shortlisted_for_deep_evaluation']=False
        if row['single_sample_full_screen']:
            newmodels[row['id']]=binder(ROOT/row['raw_scoring_output']/'model.pdb')
            rs=[]
            for name,q,model in oldmodels:
                comparison={'candidate':row['id'],'reference':name,'reference_path':q['path'],
                            'reference_sha256':q['sha256'],**reciprocal(newmodels[row['id']],model)}
                rs.append(comparison);comparisons.append(comparison)
            row['max_reciprocal_TM_to_current_draft']=max(x['reciprocal_TM'] for x in rs)
            row['structurally_redundant_to_current_draft']=row['max_reciprocal_TM_to_current_draft']>=.8
        row['eligible_for_diversity_confirmation']=r['controls_pass'] and row['single_sample_full_screen'] and not (
            row['exact_prior_sequence'] or row['similar_to_current_draft'] or row['structurally_redundant_to_current_draft'])
    pools={a:[x for x in r['ranked_candidates'] if x['arm']==a and x['eligible_for_diversity_confirmation']]
           for a in ['diversity_short','diversity_long']}
    chosen=[];within=[]
    while len(chosen)<6 and any(pools.values()):
        for pool in pools.values():
            while pool:
                row=pool.pop(0);blocked=False
                for other in chosen:
                    pair={'a':row['id'],'b':other['id'],'sequence_identity':aligned_identity(row['sequence'],other['sequence']),
                          **reciprocal(newmodels[row['id']],newmodels[other['id']])}
                    within.append(pair)
                    blocked=blocked or pair['sequence_identity']>=.8 or pair['reciprocal_TM']>=.8
                if not blocked:
                    row['shortlisted_for_deep_evaluation']=True;chosen.append(row);break
            if len(chosen)==6:break
    r['shortlisted_ids']=[x['id'] for x in chosen]
    for arm in r['arms']:r['arms'][arm]['shortlisted']=sum(x['arm']==arm for x in chosen)
    r['known_sequence_check_scope']='Exact comparison to all 495 prior designs and published controls; current-draft sequence/TM redundancy review. External database search deferred until after confirmation.'
    r['diversity_prior_catalog_sha256']=sha(run/'inputs/prior_catalog.json')
    save(out/'results.json',r)
    save(out/'diversity_audit.json',{'prospective_policy':json.loads((run/'inputs/manifest.json').read_text())['diversity_policy'],
        'prior_design_count':len(catalog['prior_designs']),'comparisons_to_draft':comparisons,'within_new_selection':within,
        'selected_ids':r['shortlisted_ids'],'selection_utc':now(),'controls_pass':r['controls_pass']})
    # Base analysis exports its old shortlist first; replace all shortlist representations.
    import csv
    p=out/'ranked_designs.csv'
    with p.open() as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
    for row in rows:row['shortlisted']=row['id'] in r['shortlisted_ids']
    with p.open('w',newline='') as f:writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    (out/'shortlist_for_evaluation.fasta').write_text(''.join(f">{x['id']}\n{x['sequence']}\n" for x in chosen))
    print('Frozen distinct confirmation candidates:',r['shortlisted_ids'])


def confirmation(run):
    report=ROOT/'reports/diversity_generation/results.json';r=json.loads(report.read_text())
    assert r['run']==run.name and r['controls_pass']
    ids=r['shortlisted_ids']
    if not ids:raise SystemExit('No distinct full-screen passes; no paid confirmation.')
    out=ROOT/'data/diversity_confirmation';out.mkdir(exist_ok=False)
    m=json.loads((run/'inputs/manifest.json').read_text())
    cases=[c for c in json.loads((run/'cases.json').read_text()) if c['id'] in ids or c['arm']=='control']
    for c in cases:
        for ext in ('.pdb','.cif','.npz'):(out/(c['id']+ext)).write_bytes((run/'references'/(c['id']+ext)).read_bytes())
    m.update(cases=cases,selected_ids=ids,source_run=r['run'],source_results_sha256=sha(report),
             prepared_utc=now(),scope='One fresh two-seed confirmation of the frozen diverse full-screen passes',
             refolding={**m['refolding'],'seeds':[2309,2411],'diffusion_samples_per_seed':5},
             max_prediction_seconds=1500,experiment_cap_usd=1.5)
    m['input_sha256']={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    save(out/'manifest.json',m);save(ROOT/'reports/diversity_confirmation/plan.json',m)
    print('Prepared confirmation:',len(ids),'candidates;',len(cases)*10,'predictions.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['screen','confirmation']);p.add_argument('run',type=__import__('pathlib').Path)
    a=p.parse_args();globals()[a.action](a.run.resolve())
