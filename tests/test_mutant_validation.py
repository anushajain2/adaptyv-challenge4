import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from prepare_controls import ROOT, read_structure, sequence
from prepare_mutant_validation import mutate_reference
from repair_mutant_input_names import pinned_target_regex


def test_mutant_reference_changes_only_requested_sequence_sites_not_backbone():
    parent = read_structure(ROOT / 'data/confirmation/apical_long_015.pdb')
    mutant = mutate_reference(parent)
    assert [i+1 for i, (a, b) in enumerate(zip(sequence(parent[0]['B']), sequence(mutant[0]['B']))) if a != b] == [51, 70]
    assert parent[0]['B'][51].resname == 'ASP' and parent[0]['B'][70].resname == 'THR'
    for chain in ('A', 'B'):
        for residue in parent[0][chain]:
            for atom in ('N', 'CA', 'C', 'O'):
                np.testing.assert_array_equal(residue[atom].coord, mutant[0][chain][residue.id][atom].coord)
    for site in (51, 70):
        assert mutant[0]['B'][site].resname == 'HIS'
        assert set(mutant[0]['B'][site].child_dict) == {'N', 'CA', 'C', 'O'}


def test_run_ids_match_pinned_parser():
    import re
    regex = pinned_target_regex()
    assert not re.fullmatch(regex, 'apical_long_015__D51H_T70H')
    for name in ('apical_long_015_D51H_T70H', 'apical_long_015', '2ds25_5_6wrw', 'scrambled_6wrw', 'caplacizumab_6wrw'):
        assert re.fullmatch(regex, name)
