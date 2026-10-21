import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import OUT, ROOT, sha
from analyze_novelty import flatten_hits, sequence_flags, structural_alignment, load_fold, parse_blast, human_tfr1_retrieved
from prepare_controls import read_structure, sequence
from tmtools.io import get_residue_data


def test_queries_are_exact_binder_only_and_include_both_failed_representatives():
    plan=json.loads((OUT/'plan.json').read_text())
    assert sha(OUT/'draft_audit_before_review.json')==plan['draft_audit_sha256']
    assert len(plan['structures'])==13 and len(plan['sequences'])==8
    failed=[]
    for q in plan['structures']:
        assert sha(ROOT/q['source'])==q['source_sha256']
        assert sha(ROOT/q['query_file'])==q['query_sha256']
        s=read_structure(ROOT/q['query_file'])
        assert [c.id for c in s[0]]==['B']
        assert sequence(s[0]['B'])==q['sequence']
        if q.get('original_screen_pass') is False:failed.append(q['candidate'])
    assert set(failed)=={'apical_short_008','expansion_apical_short_018'}


def test_sequence_review_flag_requires_significance_and_query_coverage():
    assert sequence_flags(.95,.9,1e-10)['near_identical_sequence_review_flag']
    assert not sequence_flags(.99,.3,1e-10)['close_sequence_review_flag']
    assert not sequence_flags(.65,.9,2)['close_sequence_review_flag']
    assert sequence_flags(.5,.8,1e-3)['close_sequence_review_flag']


def test_server_nested_hits_are_not_mistaken_for_one_hit():
    assert flatten_hits([[{'target':'A'},{'target':'B'}]])==[{'target':'A'},{'target':'B'}]
    p=OUT/'foldseek/published_2ds25_5_control/result.json'
    result=load_fold(p)
    hits=flatten_hits(next(r for r in result['results'] if r['db']=='pdb100')['alignments'])
    assert len(hits)>1 and hits[0]['target'].startswith('6wrw-')


def test_tm_alignment_is_rigid_motion_invariant_and_reports_both_normalizations():
    q=read_structure(OUT/'queries/apical_short_008_A41H_seed2003.pdb')
    coords,seq=get_residue_data(q[0]['B'])
    moved=coords@np.array([[0.,-1,0],[1,0,0],[0,0,1]])+np.array([50.,-20.,7.])
    r=structural_alignment(coords,moved,seq,seq)
    assert r['query_normalized_TM']>.999 and r['target_normalized_TM']>.999
    assert r['query_coverage']==r['target_coverage']==1


def test_known_fold_match_replicates_and_partial_large_protein_match_is_distinguished():
    rows=json.loads((OUT/'independent_pdb_verification.json').read_text())['rows']
    known=[r for r in rows if r['pdb']=='2OO2']
    assert len(known)==2
    for r in known:
        assert r['query_normalized_TM']>.85 and r['target_normalized_TM']>.70
        assert r['query_coverage']==1 and .8<r['target_coverage']<.81
        assert r['global_sequence_identity']<.3
    parts=[r for r in rows if r['pdb']=='8A53']
    assert all(r['query_normalized_TM']>.65 and r['target_normalized_TM']<.2 for r in parts)


def test_completed_local_blast_matches_every_frozen_query_and_recovers_human_accession():
    plan=json.loads((OUT/'plan.json').read_text())
    expected={r['id']:r['sequence'] for r in plan['sequences'] if r['id'] in plan['blast']['query_ids']}
    payload=json.loads((OUT/'local_swissprot/results.json').read_text())
    rows=parse_blast(payload,expected)
    assert len(rows)==7
    control=next(r for r in rows if r['query']=='human_tfr1_crop_control')
    assert human_tfr1_retrieved(control['hits'])
    # Matching geometry/identity from another accession cannot validate this control.
    wrong=[{**h,'descriptions':[{'id':'sp|Q5RDH6|TFR1_PONAB','accession':'Q5RDH6'}]} for h in control['hits']]
    assert not human_tfr1_retrieved(wrong)
    candidates=[r for r in rows if r['query']!='human_tfr1_crop_control']
    assert all(r['masking'] and r['hits'] for r in candidates)
    assert not any(h['close_sequence_review_flag'] for r in candidates for h in r['hits'])
    assert all(r['statistics']['db_num']==575748 for r in rows)


def test_blast_rejects_incomplete_duplicate_or_wrong_length_results():
    import copy
    payload=json.loads((OUT/'local_swissprot/results.json').read_text())
    plan=json.loads((OUT/'plan.json').read_text())
    expected={r['id']:r['sequence'] for r in plan['sequences'] if r['id'] in plan['blast']['query_ids']}
    missing=copy.deepcopy(payload);missing['BlastOutput2'].pop()
    duplicate=copy.deepcopy(payload);duplicate['BlastOutput2'][-1]=duplicate['BlastOutput2'][0]
    wrong=copy.deepcopy(payload);wrong['BlastOutput2'][0]['report']['results']['search']['query_len']+=1
    for invalid in [missing,duplicate,wrong]:
        with pytest.raises(AssertionError):parse_blast(invalid,expected)
