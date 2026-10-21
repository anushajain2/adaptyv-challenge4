"""Interpret archived database hits; do not equate absent hits with novelty."""
import gzip
import itertools
import json
import math
from pathlib import Path

import numpy as np
import tmtools
from tmtools import tm_align
from tmtools.io import get_residue_data

from prepare_controls import ROOT, read_structure
from novelty_search import OUT, save, sha, now


def load_fold(path):
    if path.exists():return json.loads(path.read_text())
    return json.loads(gzip.decompress(path.with_suffix('.json.gz').read_bytes()))


def flatten_hits(groups):
    return [hit for group in groups for hit in (group if isinstance(group,list) else [group])]


def sequence_flags(identity, query_coverage, evalue):
    return {'close_sequence_review_flag': identity>=.5 and query_coverage>=.8 and evalue<=1e-3,
            'near_identical_sequence_review_flag':identity>=.9 and query_coverage>=.8 and evalue<=1e-3}


def parse_blast(payload, expected_sequences):
    """Reject incomplete/mismatched output; flags apply to individual HSPs."""
    rows_by_query=[];seen=set()
    for item in payload['BlastOutput2']:
        report=item['report'];search=report['results']['search'];rows=[]
        name=search['query_title']
        assert report['program']=='blastp'
        assert name in expected_sequences and name not in seen
        assert search['query_len']==len(expected_sequences[name])
        seen.add(name)
        for h in search.get('hits',[]):
            for a in h['hsps']:
                ident=a['identity']/a['align_len'];qcov=(a['query_to']-a['query_from']+1)/search['query_len']
                tcov=(a['hit_to']-a['hit_from']+1)/h['len']
                rows.append({'descriptions':h['description'],'subject_length':h['len'],**a,
                             'sequence_identity':ident,'query_coverage':qcov,'target_coverage':tcov,
                             **sequence_flags(ident,qcov,a['evalue'])})
        rows_by_query.append({'query':name,'status':'complete','query_length':search['query_len'],
                             'statistics':search.get('stat'),'masking':search.get('query_masking',[]),
                             'message':search.get('message'),'database':report['search_target'],
                             'program':report['program'],'version':report['version'],'parameters':report['params'],
                             'hits':rows})
    assert seen==set(expected_sequences), 'Missing BLAST query results'
    return rows_by_query


def human_tfr1_retrieved(rows):
    """Require human TfR1 accession P02786, not merely a similar control hit."""
    import re
    return any(h['sequence_identity']>=.99 and h['query_coverage']>=.99 and h['evalue']<=1e-10
               and any('P02786' in re.split(r'[|.\s]+',str(d.get('accession',''))+' '+str(d.get('id','')))
                       for d in h['descriptions']) for h in rows)


def structural_alignment(x, y, sx, sy):
    assert x.shape==(len(sx),3) and y.shape==(len(sy),3)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    r=tm_align(np.asarray(x,dtype=np.float64),np.asarray(y,dtype=np.float64),sx,sy)
    pairs=[(a,b) for a,b in zip(r.seqxA,r.seqyA) if a!='-' and b!='-']
    # TM-align's colon marks distances under 5 A; period marks other aligned pairs.
    return {'query_normalized_TM':r.tm_norm_chain1,'target_normalized_TM':r.tm_norm_chain2,
            'RMSD_A':r.rmsd,'aligned_pairs':len(pairs),'query_coverage':len(pairs)/len(sx),
            'target_coverage':len(pairs)/len(sy),
            'aligned_identity':sum(a==b for a,b in pairs)/len(pairs) if pairs else 0,
            'aligned_pairs_within5A':r.seqM.count(':'),
            'query_alignment':r.seqxA,'target_alignment':r.seqyA,'distance_markers':r.seqM}


def analyze():
    plan=json.loads((OUT/'plan.json').read_text());foldrows=[];queries={}
    for c in plan['structures']:
        s=read_structure(ROOT/c['query_file']);x,sx=get_residue_data(s[0]['B']);assert sx==c['sequence']
        queries[c['id']]=(x,sx)
        p=OUT/'foldseek'/c['id']/'result.json'
        if not p.exists() and not p.with_suffix('.json.gz').exists():
            foldrows.append({'query':c['id'],'status':'missing'});continue
        result=load_fold(p)
        assert len(result['queries'])==1 and result['queries'][0]['sequence']==sx
        for db in result['results']:
            hits=flatten_hits(db['alignments']);rows=[]
            for rank,h in enumerate(hits,1):
                qcov=sum(a!='-' for a in h['qAln'])/h['qLen']
                tcov=sum(a!='-' for a in h['dbAln'])/h['dbLen']
                row={k:v for k,v in h.items() if k not in ['tCa','tSeq']}
                row.update(server_rank=rank,query_coverage=qcov,target_coverage=tcov)
                if rank<=10:
                    coords=np.fromstring(h['tCa'],sep=',',dtype=float).reshape(-1,3)
                    assert len(coords)==len(h['tSeq'])==h['dbLen']
                    row['local_TM_align']=structural_alignment(x,coords,sx,h['tSeq'])
                    # Archive only follow-up target coordinates separately; full hits remain compressed.
                    target_file=OUT/'followup_targets'/f"{c['id']}__{db['db']}__rank{rank}.npz"
                    target_file.parent.mkdir(exist_ok=True)
                    np.savez_compressed(target_file,coords=coords,sequence=np.array(h['tSeq']))
                    row['followup_target_file']=str(target_file.relative_to(ROOT))
                    row['followup_target_sha256']=sha(target_file)
                rows.append(row)
            foldrows.append({'query':c['id'],'candidate':c['candidate'],'kind':c['kind'],'database':db['db'],
                             'status':'complete','hit_count':len(rows),'hits':rows,
                             'near_server_result_limit':len(rows)>=950})
        print('Structural alignments complete:',c['id'],flush=True)
    rcsb=[]
    for c in plan['sequences']:
        path=OUT/'rcsb'/c['id']/'response.json'
        if not path.exists():rcsb.append({'query':c['id'],'status':'missing'});continue
        response=json.loads(path.read_text());rows=[]
        for h in response['data'].get('result_set',[]):
            for service in h['services']:
                for node in service['nodes']:
                    for a in node.get('match_context',[]):
                        qcov=(a['query_end']-a['query_beg']+1)/a['query_length']
                        tcov=(a['subject_end']-a['subject_beg']+1)/a['subject_length']
                        rows.append({'target':h['identifier'],**a,'query_coverage':qcov,'target_coverage':tcov,
                                     **sequence_flags(a['sequence_identity'],qcov,a['evalue'])})
        rcsb.append({'query':c['id'],'kind':c['kind'],'status':'complete','http_status':response['http_status'],
                     'total_count':response['data'].get('total_count',0),'returned_hit_count':len(rows),'hits':rows})
    blast=[]
    bp=OUT/'blast/results.json'
    source='NCBI_remote'
    if not bp.exists() and (OUT/'local_swissprot/results.json').exists():
        bp=OUT/'local_swissprot/results.json';source='local_UniProtKB_SwissProt_canonical'
    if bp.exists():
        expected={r['id']:r['sequence'] for r in plan['sequences'] if r['id'] in plan['blast']['query_ids']}
        blast=parse_blast(json.loads(bp.read_text()),expected)
    else:blast=[{'status':'pending_or_missing'}]
    internal=[]
    for a,b in itertools.combinations([c for c in plan['structures'] if c['kind']=='candidate'],2):
        x,sx=queries[a['id']];y,sy=queries[b['id']]
        internal.append({'a':a['id'],'b':b['id'],'candidate_a':a['candidate'],'candidate_b':b['candidate'],
                         **structural_alignment(x,y,sx,sy)})
    control=next(r for r in foldrows if r['query']=='published_2ds25_5_control' and r.get('database')=='pdb100')
    expected=next((h for h in control['hits'][:10] if h['target'].startswith('6wrw-')),None)
    control_pass=expected is not None and expected['local_TM_align']['query_normalized_TM']>=.95 and expected['local_TM_align']['query_coverage']>=.95
    result={'analyzed_utc':now(),'plan_sha256':sha(OUT/'plan.json'),'tmtools_version':tmtools.__version__,
            'foldseek_positive_retrieval_pass':control_pass,'foldseek':foldrows,'rcsb_sequence':rcsb,
            'swissprot_sequence':blast,'swissprot_sequence_source':source if bp.exists() else None,
            'swissprot_result_file':str(bp.relative_to(ROOT)) if bp.exists() else None,
            'internal_structure_alignments':internal,
            'limits':'Database search and TM-align comparisons are descriptive, not proof of de novo provenance, absence of homologs, novelty, binding, or fold-family membership. Local TM-align is restricted to top 10 server-ranked hits per database/query.'}
    save(OUT/'analysis.json',result)
    print('Foldseek retrieval control pass:',control_pass,'BLAST reports:',len(blast))


if __name__=='__main__':analyze()
