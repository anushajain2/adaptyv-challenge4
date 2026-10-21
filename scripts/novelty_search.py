"""Archived public-database similarity searches. No Modal jobs or submissions.

Commands are deliberately separate: prepare freezes queries/interpretation, submit
starts bounded public searches, collect polls once, and rcsb runs sequence searches.
"""
import argparse
import datetime as dt
import hashlib
import gzip
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from Bio.PDB import PDBIO, Select
from prepare_controls import ROOT, read_structure, sequence

OUT = ROOT/'reports/novelty_review'
FOLD = 'https://search.foldseek.com/api'
BLAST = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
RCSB = 'https://search.rcsb.org/rcsbsearch/v2/query'


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2)+'\n')


def request(url, data=None, json_body=False):
    headers={'User-Agent':'TfR1-similarity-review/1.0','Cache-Control':'no-cache'}
    if data is not None:
        headers['Content-Type']='application/json' if json_body else 'application/x-www-form-urlencoded'
        data=json.dumps(data).encode() if json_body else urllib.parse.urlencode(data, doseq=True).encode()
    req=urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(req,timeout=55) as response:
        return response.status, dict(response.headers), response.read()


class BinderOnly(Select):
    def accept_chain(self, chain):
        return chain.id=='B'


def prepare():
    if (OUT/'plan.json').exists():
        raise SystemExit('Frozen plan already exists; use its queries unchanged.')
    draft=json.loads((ROOT/'reports/submission_draft/audit.json').read_text())
    manifest=json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())
    structures=[]
    for c in draft['candidates']:
        for model in c['latest_selected_bare_models']:
            structures.append({'id':f"{c['name']}_seed{model['seed']}", 'candidate':c['name'],
                               'kind':'candidate', 'sequence':c['sequence'], 'source':model['path'],
                               'source_sha256':model['sha256'], 'original_screen_pass':model['selected_screen_pass']})
    ctrl=next(c for c in manifest['cases'] if c['id']=='2ds25_5_6wrw_bare')
    source='data/tag_comparison/2ds25_5_6wrw_bare.pdb'
    structures.append({'id':'published_2ds25_5_control','candidate':None,'kind':'retrieval_control',
                       'sequence':ctrl['core_sequence'],'source':source,'source_sha256':sha(ROOT/source)})
    for c in structures:
        s=read_structure(ROOT/c['source']);assert sequence(s[0]['B'])==c['sequence']
        p=OUT/'queries'/f"{c['id']}.pdb";p.parent.mkdir(parents=True,exist_ok=True)
        io=PDBIO();io.set_structure(s);io.save(str(p),BinderOnly())
        query=read_structure(p)
        assert [ch.id for ch in query[0]]==['B']
        assert sequence(query[0]['B'])==c['sequence']
        assert len([r for r in query[0]['B'] if 'CA' in r])==len(c['sequence'])
        c['query_file']=str(p.relative_to(ROOT));c['query_sha256']=sha(p)
    sequences=[{'id':c['name'],'sequence':c['sequence'],'kind':'candidate'} for c in draft['candidates']]
    sequences.append({'id':'human_tfr1_crop_control','sequence':manifest['target_sequence'],'kind':'retrieval_control'})
    sequences.append({'id':'published_2ds25_5_control','sequence':ctrl['core_sequence'],'kind':'retrieval_control'})
    dbs=json.loads((ROOT/'outputs/novelty_api_docs/foldseek_databases.json').read_text())
    save(OUT/'database_inventory.json',dbs)
    selected=[d for d in dbs['databases'] if d['path'] in ['pdb100','afdb-swissprot']]
    assert len(selected)==2 and all(d['status']=='COMPLETE' for d in selected)
    plan={'prepared_utc':now(),'draft_audit_sha256':sha(ROOT/'reports/submission_draft/audit.json'),
          'structures':structures,'sequences':sequences,
          'foldseek':{'mode':'3diaa','databases':selected,'iterativesearch':False,
                      'scope':'Both selected bare poses per draft candidate; receptor and tags removed. Positive retrieval control is the published 2DS25.5 binder.'},
          'rcsb':{'sequence_type':'protein','identity_cutoff':0,'evalue_cutoff':10,
                  'return_type':'polymer_entity','result_content_type':['experimental'],'max_results':100},
          'blast':{'PROGRAM':'blastp','DATABASE':'swissprot','EXPECT':10,'FILTER':'mL',
                   'COMPOSITION_BASED_STATISTICS':2,'MATRIX':'BLOSUM62','WORD_SIZE':3,'HITLIST_SIZE':100,
                   'query_ids':[r['id'] for r in sequences if r['id']!='published_2ds25_5_control']},
          'interpretation':{
              'close_sequence_review_flag':'Identity >=50%, query coverage >=80%, E-value <=1e-3; report subject coverage separately. A review heuristic, not a competition rule.',
              'near_identical_sequence_review_flag':'Identity >=90%, query coverage >=80%, E-value <=1e-3.',
              'structural_followup':'Inspect the first 10 ranked hits per query/database; obtain structure-aware alignments and query/target-normalized TM scores where available. Common folds do not imply copied sequence or binding function.',
              'limits':'No hit is not proof of novelty. PDB100 server version is 20240101; Swiss-Prot is curated and incomplete, and AlphaFold models are predictions. Search failures must remain missing evidence.',
              'selection':'No candidate is removed solely for a shared fold. Preserve original submission sequences and rank pending interpretation.'},
          'public_apis':[FOLD,BLAST,RCSB],'new_modal_compute_usd':0}
    save(OUT/'plan.json',plan)
    print('Frozen',len(structures),'binder-only structures and',len(sequences),'sequence queries.')


def submit():
    p=json.loads((OUT/'plan.json').read_text())
    # A single multi-query Swiss-Prot request; never retry an uncertain POST.
    bp=OUT/'blast/submission.json'
    if not bp.exists():
        params=dict(p['blast']);ids=params.pop('query_ids')
        params.update(CMD='Put',QUERY=''.join(f">{r['id']}\n{r['sequence']}\n" for r in p['sequences'] if r['id'] in ids))
        save(OUT/'blast/request.json',params)
        save(bp,{'state':'POST_STARTED','utc':now()})
        status,headers,raw=request(BLAST,params)
        (OUT/'blast/submission_response.html').write_bytes(raw)
        rid=re.search(r'RID\s*=\s*([A-Z0-9-]+)',raw.decode())
        wait=re.search(r'RTOE\s*=\s*(\d+)',raw.decode())
        save(bp,{'state':'SUBMITTED' if rid else 'RESPONSE_WITHOUT_RID','utc':now(),
                 'RID':rid.group(1) if rid else None,'RTOE':int(wait.group(1)) if wait else None,
                 'http_status':status,'response_sha256':sha(OUT/'blast/submission_response.html')})
        print('BLAST',json.loads(bp.read_text()),flush=True)
    for c in p['structures']:
        dest=OUT/'foldseek'/c['id'];ticket=dest/'ticket.json'
        if ticket.exists():continue
        assert sha(ROOT/c['query_file'])==c['query_sha256']
        params={'q':(ROOT/c['query_file']).read_text(),'database[]':[d['path'] for d in p['foldseek']['databases']],
                'mode':p['foldseek']['mode'],'iterativesearch':'false'}
        save(ticket,{'state':'POST_STARTED','utc':now()})
        try:
            status,headers,raw=request(FOLD+'/ticket',params)
            response=json.loads(raw);save(ticket,{'utc':now(),'http_status':status,**response})
            print(c['id'],response,flush=True)
        except Exception as e:
            save(ticket,{'state':'SUBMISSION_UNCERTAIN','utc':now(),'error':str(e)});print(c['id'],str(e),flush=True)
        time.sleep(1)


def collect():
    for tp in sorted((OUT/'foldseek').glob('*/ticket.json')):
        if (tp.parent/'result.json').exists() or (tp.parent/'result.json.gz').exists():continue
        t=json.loads(tp.read_text())
        if 'id' not in t:continue
        try:
            _,_,raw=request(FOLD+'/ticket/'+t['id']);s=json.loads(raw)
            save(tp.parent/'status.json',{'checked_utc':now(),**s});print(tp.parent.name,s['status'],flush=True)
            if s['status']=='COMPLETE':
                _,_,raw=request(FOLD+'/result/'+t['id']+'/0')
                json.loads(raw);(tp.parent/'result.json.gz').write_bytes(gzip.compress(raw,mtime=0))
        except Exception as e:save(tp.parent/'collection_error.json',{'utc':now(),'error':str(e)})
    bp=OUT/'blast/submission.json'
    if bp.exists() and not (OUT/'blast/results.json').exists():
        b=json.loads(bp.read_text());rid=b.get('RID')
        if not rid:return
        last=OUT/'blast/last_poll.json'
        if last.exists() and time.time()-json.loads(last.read_text())['epoch']<65:
            print('BLAST polling interval not yet elapsed');return
        save(last,{'epoch':time.time(),'utc':now()})
        status,headers,raw=request(BLAST+'?'+urllib.parse.urlencode({'CMD':'Get','RID':rid,'FORMAT_TYPE':'JSON2_S'}))
        save(OUT/'blast/last_poll_response.json',{'utc':now(),'status':status,'headers':headers,'body_sha256':hashlib.sha256(raw).hexdigest()})
        try:
            result=json.loads(raw)
            if 'BlastOutput2' in result:
                (OUT/'blast/results.json').write_bytes(raw);print('BLAST complete',flush=True)
            else:save(OUT/'blast/pending.json',result)
        except json.JSONDecodeError:
            (OUT/'blast/pending.html').write_bytes(raw);print('BLAST not yet JSON result',flush=True)


def rcsb():
    p=json.loads((OUT/'plan.json').read_text());opts=p['rcsb']
    for c in p['sequences']:
        dest=OUT/'rcsb'/c['id'];result=dest/'response.json'
        if result.exists():continue
        query={'query':{'type':'terminal','service':'sequence','parameters':{'sequence_type':'protein',
                'value':c['sequence'],'identity_cutoff':opts['identity_cutoff'],'evalue_cutoff':opts['evalue_cutoff']}},
               'return_type':'polymer_entity','request_options':{'results_verbosity':'verbose',
                'paginate':{'start':0,'rows':opts['max_results']},'results_content_type':opts['result_content_type']}}
        save(dest/'request.json',query)
        try:
            status,headers,raw=request(RCSB,query,json_body=True)
            data=json.loads(raw) if raw else {}
            save(result,{'retrieved_utc':now(),'http_status':status,'response_headers':headers,'data':data})
            print(c['id'],'hits',data.get('total_count',0),flush=True)
        except urllib.error.HTTPError as e:
            save(dest/'error.json',{'utc':now(),'status':e.code,'body':e.read().decode()});print(c['id'],e,flush=True)
        time.sleep(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','submit','collect','rcsb'])
    globals()[parser.parse_args().action]()
