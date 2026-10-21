import copy
import csv
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from build_submission_draft import (OUT, ROOT, SUFFIX, alignment, digest,
                                    sequence_audit, validate_submission)
from prepare_controls import read_structure, sequence


def test_draft_sequences_are_bare_unique_provenanced_and_within_track_limit():
    audit = json.loads((OUT/'audit.json').read_text())
    rows = list(csv.DictReader((OUT/'submission_draft.csv').open()))
    sources = json.loads((ROOT/'reports/expansion_his_confirmation/results.json').read_text())
    sources2 = json.loads((ROOT/'reports/expansion_validation/results.json').read_text())
    seqs = {c['id']: c['sequence'] for d in (sources, sources2) for c in d['candidate_summaries']}
    controls = {c['core_sequence'] for c in json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())['cases'] if c['kind']!='novel_design'}
    validate_submission(rows, seqs, controls)
    assert [r['name'] for r in rows] == [c['name'] for c in audit['candidates']]
    assert len(rows)==6 and audit['parent_lineage_count']==5
    assert audit['independent_fold_class_count'] is None
    assert 'E19H' not in ''.join(r['name'] for r in rows)


def test_export_rejects_appended_suffix_controls_and_duplicate_sequences():
    seq='ACDEFGHIKLMNPQRSTVWY'
    row={'name':'test','sequence':seq,'molecule_class':'protein'}
    validate_submission([row], {'test':seq}, set())
    with pytest.raises(AssertionError):
        validate_submission([row], {'test':seq}, {seq})
    tagged=copy.deepcopy(row);tagged['sequence']+=SUFFIX
    with pytest.raises(AssertionError):
        validate_submission([tagged], {'test':tagged['sequence']}, set())
    second={**row,'name':'other'}
    with pytest.raises(AssertionError):
        validate_submission([row,second], {'test':seq,'other':seq}, set())


def test_sequence_flags_include_overlaps_and_exclude_proline_sequon():
    a=sequence_audit('NNSTNPTAAAA')
    assert a['N_X_S_or_T_sequons']==[{'start_1based':1,'sequence':'NNS'},{'start_1based':2,'sequence':'NST'}]
    assert a['identical_residue_runs_ge4']==[{'start_1based':8,'sequence':'AAAA'}]
    assert abs(sequence_audit('EEEEEEEEEE')['sequence_charge_pH7_4'])>9
    with pytest.raises(ValueError):sequence_audit('ACDX')


def test_alignment_uses_longer_sequence_denominator_not_overlap_identity():
    pairs, score=alignment('ACDEFGHIKL','ACDEFGHIKLAAAA')
    assert len(pairs)==10 and score==10/14
    audit=json.loads((OUT/'audit.json').read_text())
    related=[r for r in audit['pairwise_diversity'] if r['same_parent_lineage']]
    assert len(related)==1 and related[0]['global_sequence_identity']==57/58


def test_bundle_preserves_failed_selected_models_and_direction_exclusions():
    audit=json.loads((OUT/'audit.json').read_text())
    failures=[]
    for c in audit['candidates']:
        assert c['experimental_binding']=='unmeasured'
        assert c['tag_compatibility']=='unresolved'
        for m in c['latest_selected_bare_models']:
            assert sequence(read_structure(ROOT/m['path'])[0]['B'])==c['sequence']
            if not m['selected_screen_pass']: failures.append(c['name'])
        for v in c['latest_direction']:
            assert not v['fully_qualified_direction_pass_both_seeds']
    assert set(failures)=={'apical_short_008','expansion_apical_short_018'}
    a41h=next(c for c in audit['candidates'] if c['name']=='apical_short_008_A41H')
    assert any(m['status']=='excluded' for v in a41h['latest_direction'] for m in v['models'])


def test_manifest_hashes_match_sources_models_and_outputs():
    manifest=json.loads((OUT/'manifest.json').read_text())
    for p,h in manifest['inputs'].items():assert digest(ROOT/p)==h
    for p,h in manifest['outputs'].items():assert digest(OUT/p)==h
    for m in manifest['models']:
        assert digest(ROOT/m['path'])==m['sha256']
        if (ROOT/m['source']).exists():
            assert digest(ROOT/m['source'])==m['sha256']
