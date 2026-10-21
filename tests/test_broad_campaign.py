import json
import sys
from pathlib import Path
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from launch_broad_campaign import reservation_bound,bundled_files

def test_authorized_counts_epitopes_and_nonoverlapping_hotspots():
    files=bundled_files();m=json.loads(files['manifest.json'])
    assert sum(a['num_designs'] for a in m['arms'])==200
    assert sum(a['num_designs'] for a in m['arms'] if a['protocol']=='nanobody-anything')==40
    for arm in m['arms']:
        spec=yaml.safe_load(files[arm['id']+'.yaml'])
        entries=spec['entities'][0]['file']['binding_types']
        assert all(not {'binding','not_binding'}<=set(e['chain']) for e in entries)
        excluded=set(map(int,entries[0]['chain']['not_binding'].split(',')))
        hot=set(map(int,entries[1]['chain']['binding'].split(',')))
        assert hot==set(arm['hotspots_crop']) and not hot&excluded
        if arm['epitope']=='conserved':
            assert hot=={138,142}
            assert set(arm['reference_contact_region_crop'])==set(range(137,159))
        if arm['protocol']=='protein-anything':
            low,high=arm['binder_length_range']
            assert not set(range(low,high+1))&set(range(73,77))

def test_reserved_resource_cost_fits_first_batch_and_total():
    assert reservation_bound()<6
    assert 3.8935808+16+7+3<30

from analyze_broad_campaign import rank_candidates,diverse_shortlist

def candidate(name,protocol,score,eligible=True):
    return {'id':name,'arm':protocol,'protocol':protocol,'geometry_eligible':eligible,
        'sequence':'ACDEFGHIKLMNPQRSTVWY','metrics':{'design_to_target_iptm':score},
        'scores':{'contact_pae_both_directions_A':20-score*10},'hotspot_fraction':score,
        'pose':{'binder_CA_fold_rmsd_A':2-score}}

def test_ranking_is_within_class_and_geometry_precedes_confidence():
    a=candidate('mini_low','protein-anything',.3)
    b=candidate('mini_high','protein-anything',.9,False)
    c=candidate('nano','nanobody-anything',.2)
    ranked=rank_candidates([a,b,c])
    assert ranked[-1]['id']=='mini_high'
    assert c['within_class_rank_score']==1
    assert a['within_class_rank_score']==0

def test_shortlist_excludes_invalid_geometry_and_near_duplicate_sequences():
    ranked=rank_candidates([candidate('a','protein-anything',.8),candidate('b','protein-anything',.7),
                            candidate('c','nanobody-anything',.9,False)])
    assert len(diverse_shortlist(ranked))==1
    assert diverse_shortlist([candidate('d','protein-anything',.9,False)])==[]

def test_known_control_sequence_cannot_enter_novel_shortlist():
    r=candidate('copied','protein-anything',.9)
    r['matches_known_sequence']=True
    assert diverse_shortlist(rank_candidates([r]))==[]

def test_vhh_diversity_uses_redesigned_loops_and_still_rejects_exact_duplicates():
    a=candidate('vhh_a','nanobody-anything',.8)
    b=candidate('vhh_b','nanobody-anything',.7)
    a['sequence']='A'*100+'Y'*20; b['sequence']='A'*100+'F'*20
    a['diversity_sequence']='Y'*20; b['diversity_sequence']='F'*20
    assert len(diverse_shortlist(rank_candidates([a,b])))==2
    b['sequence']=a['sequence']
    assert len(diverse_shortlist(rank_candidates([a,b])))==1
