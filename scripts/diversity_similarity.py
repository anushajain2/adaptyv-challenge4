"""Reuse the existing bounded database review for confirmed new parents only."""
import argparse
import copy
import json
import subprocess
from pathlib import Path
from Bio.PDB import PDBIO
import novelty_search as ns
from novelty_search import ROOT, sha, save, now, BinderOnly
from prepare_controls import read_structure, sequence

OUT=ROOT/'reports/diversity_novelty'


def prepare():
    assert not OUT.exists(),'Do not overwrite a frozen review.'
    summary=json.loads((ROOT/'reports/diversity_confirmation/summary.json').read_text())
    candidates=[c for c in summary['candidates'] if c['two_seed_raw_screen_pass']]
    if not candidates:raise SystemExit('No confirmed new parents; no additional search.')
    base=json.loads((ROOT/'reports/novelty_review/plan.json').read_text())
    plan=copy.deepcopy(base);structures=[]
    for c in candidates:
        for model in c['selected_models']:
            structures.append({'id':f"{c['id']}_seed{model['seed']}",'candidate':c['id'],
                'kind':'candidate','sequence':c['sequence'],'source':model['path'],
                'source_sha256':model['sha256'],'original_screen_pass':model['raw_pass']})
    structures.append(copy.deepcopy(next(c for c in base['structures'] if c['kind']=='retrieval_control')))
    for c in structures:
        assert sha(ROOT/c['source'])==c['source_sha256']
        s=read_structure(ROOT/c['source']);assert sequence(s[0]['B'])==c['sequence']
        p=OUT/'queries'/f"{c['id']}.pdb";p.parent.mkdir(parents=True,exist_ok=True)
        io=PDBIO();io.set_structure(s);io.save(str(p),BinderOnly())
        c['query_file']=str(p.relative_to(ROOT));c['query_sha256']=sha(p)
    seqs=[{'id':c['id'],'kind':'candidate','sequence':c['sequence']} for c in candidates]
    seqs+=copy.deepcopy([c for c in base['sequences'] if c['kind']=='retrieval_control'])
    plan.update(prepared_utc=now(),structures=structures,sequences=seqs,
        parent_review_plan_sha256=sha(ROOT/'reports/novelty_review/plan.json'),
        confirmation_summary_sha256=sha(ROOT/'reports/diversity_confirmation/summary.json'))
    plan.pop('draft_audit_sha256',None)
    plan['blast']['query_ids']=[c['id'] for c in seqs if c['id']!='published_2ds25_5_control']
    plan['interpretation']['selection']='Review all raw two-seed passing parents. Preparation caveats remain visible. Shared folds are descriptive; no automatic submission admission.'
    save(OUT/'plan.json',plan)
    save(OUT/'database_inventory.json',json.loads((ROOT/'reports/novelty_review/database_inventory.json').read_text()))
    # Prevent another remote NCBI submission; use the already verified local database.
    save(OUT/'blast/submission.json',{'state':'LOCAL_BLAST_ONLY','reason':'Use verified local Swiss-Prot canonical release 2026_03; do not submit another remote RID.'})
    print('Prepared',len(candidates),'new sequences and',len(structures),'binder-only structure queries.')


def local():
    plan=json.loads((OUT/'plan.json').read_text());dest=OUT/'local_swissprot';dest.mkdir(exist_ok=False)
    old=ROOT/'reports/novelty_review/local_swissprot';oldrun=json.loads((old/'run.json').read_text())
    oldplan=json.loads((old/'plan.json').read_text());work=ROOT/'outputs/local_swissprot'
    blast=work/'ncbi-blast-2.17.0+'/ 'bin/blastp'
    assert sha(blast)==oldrun['blastp_sha256'] and sha(work/'swissprot.fasta')==oldrun['database_fasta_sha256']
    query=''.join(f">{c['id']}\n{c['sequence']}\n" for c in plan['sequences'] if c['id'] in plan['blast']['query_ids'])
    (dest/'queries.fasta').write_text(query)
    save(dest/'plan.json',{**oldplan,'prepared_utc':now(),'query_sha256':sha(dest/'queries.fasta'),
         'reason':'Same local database, executable and parameters as the completed initial review; new confirmed candidates.',
         'original_plan_sha256':sha(OUT/'plan.json'),'prior_run_sha256':sha(old/'run.json')})
    cmd=[str(blast),'-query',str(dest/'queries.fasta'),'-db',str(work/'swissprot'),'-out',str(dest/'results.json')]
    for k,v in oldplan['parameters'].items():cmd+=['-'+k,str(v)]
    record={'started_utc':now(),'state':'RUNNING','command':cmd,'database_release':oldplan['release'],
            'database_sequences':oldrun['database_sequences'],'database_residues':oldrun['database_residues'],
            'database_fasta_sha256':oldrun['database_fasta_sha256'],'blastp_sha256':sha(blast)}
    save(dest/'run.json',record)
    with (dest/'stdout.txt').open('w') as stdout,(dest/'stderr.txt').open('w') as stderr:
        subprocess.run(cmd,stdout=stdout,stderr=stderr,check=True)
    from analyze_novelty import parse_blast,human_tfr1_retrieved
    expected={c['id']:c['sequence'] for c in plan['sequences'] if c['id'] in plan['blast']['query_ids']}
    rows=parse_blast(json.loads((dest/'results.json').read_text()),expected)
    assert human_tfr1_retrieved(next(r for r in rows if r['query']=='human_tfr1_crop_control')['hits'])
    record.update(state='COMPLETE',finished_utc=now(),result_sha256=sha(dest/'results.json'))
    save(dest/'run.json',record)
    print('Local Swiss-Prot complete:',len(rows),'queries; human control pass.')


def public_check():
    # Anonymous reads establish that every exact query is already public before API export.
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    plan=json.loads((OUT/'plan.json').read_text());paths=[c['query_file'] for c in plan['structures']]+[str((OUT/'local_swissprot/queries.fasta').relative_to(ROOT))]
    records=[]
    for path in paths:
        url=f'https://raw.githubusercontent.com/anushajain2/adaptyv-challenge4/{commit}/{path}'
        status,headers,raw=ns.request(url)
        assert status==200 and raw==(ROOT/path).read_bytes()
        records.append({'path':path,'url':url,'sha256':sha(ROOT/path)})
    save(OUT/'public_provenance_check.json',{'checked_utc':now(),'commit':commit,'all_anonymous_bytes_match':True,'files':records})
    print('All exact sequence and structure queries are verified public.')


def remote():
    proof=json.loads((OUT/'public_provenance_check.json').read_text())
    assert proof['all_anonymous_bytes_match']
    for r in proof['files']:assert sha(ROOT/r['path'])==r['sha256']
    ns.OUT=OUT;ns.submit();ns.rcsb()


def collect():
    ns.OUT=OUT;ns.collect()


def analyze():
    import analyze_novelty
    analyze_novelty.OUT=OUT;analyze_novelty.analyze()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','local','public_check','remote','collect','analyze'])
    globals()[p.parse_args().action]()
