"""Build a conservative, evidence-linked draft; never submit or launch compute.

CSV is a machine-readable competition upload, not a presentation workbook.
Ranking is an explicit editorial decision; sequence descriptors do not retune gates.
"""
import csv
import hashlib
import itertools
import json
import math
import platform
import re
import shutil
from collections import Counter
from pathlib import Path

import Bio
from Bio.Align import PairwiseAligner
from Bio.PDB import Superimposer
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

from prepare_controls import ROOT, read_structure, sequence
from prepare_assay_tag_review import SUFFIX

OUT = ROOT / 'reports/submission_draft'
ROUNDS = ['his_confirmation', 'expansion_validation', 'expansion_his_confirmation', 'tag_comparison']
SHORTLIST = [
    ('apical_short_008_A41H', 'apical_short_008', 'pH-direction hypothesis',
     'Leading repeated bare-model pKa-lowering hypothesis. Latest bare core screen passes both selected seeds; tagged core fails one. No qualified two-pose direction pass.'),
    ('expansion_apical_long_092', 'expansion_apical_long_092', 'binding hypothesis',
     'Both selected models pass in two fresh confirmation rounds (9/10 then 8/10 samples). Different parent lineage; pH release and tag compatibility unestablished.'),
    ('expansion_apical_short_026', 'expansion_apical_short_026', 'binding hypothesis',
     'Repeated bare structural retention (8/10, 9/10, then 9/10). Tagged core passes both selected seeds, but tag independence remains unresolved.'),
    ('expansion_apical_short_015', 'expansion_apical_short_015', 'diversity backup',
     'Both selected confirmation models pass (6/10 samples); less repeated evidence than the leading parents. Known-backbone similarity to AF1782/2OO2; no fold-novelty claim. No tag or pH validation.'),
    ('expansion_apical_short_018', 'expansion_apical_short_018', 'mixed-geometry diversity backup',
     'Initial fresh confirmation passes both selected seeds (8/10); later confirmation fails one selected clash check (6/10). Retained as an explicitly weaker lineage.'),
    ('apical_short_008', 'apical_short_008', 'matched-parent comparator',
     'Same lineage as A41H; provides a mechanistic comparator rather than another independent scaffold. Latest bare panel fails one selected core clash check (6/10).'),
]
RESERVE = {
    'expansion_apical_short_026_E19H': 'Secondary pH hypothesis: tagged interface pose fails in both seeds; latest bare direction has one excluded pose and one strict heuristic failure. Same lineage as rank 3.',
    'expansion_apical_long_092_R3H': 'Opposite pKa-shift signs across poses; no consistent acidic-release direction.',
    'expansion_apical_long_092_E53H': 'Fresh selected interface-pose failure.',
    'expansion_apical_short_018_K3H': 'Fresh selected clash failure.',
    'expansion_apical_short_018_E57H_K60H': 'Uniform preparation gives no consistent pair direction; terminal context also unresolved.',
    'expansion_apical_short_026_V60H': 'Small pKa shifts and proximity to the C terminus; adds no independent lineage.',
    'apical_long_015': 'Later selected confirmation fails; lower priority than the retained parents.',
    'apical_long_015_L7H_L82H': 'Inconsistent direction; no compelling advantage over the retained shortlist.',
    'apical_long_015_F46H_T70H': 'Selected structural confirmation failure.',
    'apical_short_008_D14H_S47H': 'Selected structural confirmation failure.',
    'apical_short_008_A41H_S47H': 'Selected structural confirmation failure; retain the simpler A41H hypothesis.',
    'expansion_apical_long_020': 'Fresh selected confirmation failure.',
    'expansion_apical_long_089': 'Fresh selected confirmation failure.',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def positions(seq, letters):
    return [i + 1 for i, aa in enumerate(seq) if aa in letters]


def sequence_audit(seq):
    if not seq or set(seq) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('Expected a nonempty canonical protein sequence')
    p = ProteinAnalysis(seq)
    counts = Counter(seq)
    motifs = lambda pattern: [{'start_1based': m.start() + 1, 'sequence': m.group(1)}
                              for m in re.finditer('(?=(' + pattern + '))', seq)]
    runs = [{'start_1based': m.start()+1, 'sequence': m.group()}
            for m in re.finditer(r'(.)\1{3,}', seq)]
    return {
        'length': len(seq), 'canonical_sequence': True,
        'molecular_weight_Da': p.molecular_weight(),
        'estimated_pI': IsoelectricPoint(seq).pi(pH=7, min_=0, max_=14),
        'sequence_charge_pH7_4': p.charge_at_pH(7.4),
        'sequence_charge_pH6_0': p.charge_at_pH(6.0),
        'GRAVY_Kyte_Doolittle': p.gravy(),
        'amino_acid_counts': dict(sorted(counts.items())),
        'composition_entropy_bits': -sum(n/len(seq)*math.log2(n/len(seq)) for n in counts.values()),
        'Cys_positions': positions(seq, 'C'), 'Met_positions': positions(seq, 'M'),
        'Trp_positions': positions(seq, 'W'), 'His_positions': positions(seq, 'H'),
        'N_X_S_or_T_sequons': motifs('N[^P][ST]'),
        'NG_NS_NT_DG_DS_motifs': motifs('N[GST]|D[GS]'),
        'identical_residue_runs_ge4': runs,
        'hydrophobic_runs_ge5_AVILMFWY': [m.group() for m in re.finditer('[AVILMFWY]{5,}', seq)],
        'charged_fraction_DEKR': sum(counts[x] for x in 'DEKR')/len(seq),
        'extinction_coefficient_reduced_M_inv_cm_inv': p.molar_extinction_coefficient()[0],
        'interpretation': 'Descriptive bare-sequence flags only. Motifs are not demonstrated chemical liabilities; pI/charge are not structure-based His pKas or binding predictions. No expression, solubility, aggregation, immunogenicity or stability classifier applied.',
    }


def alignment(a, b):
    obj = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                         open_gap_score=-4, extend_gap_score=-.5).align(a, b)[0]
    pairs = [(int(i), int(j)) for i, j in zip(*obj.indices) if i >= 0 and j >= 0]
    return pairs, sum(a[i] == b[j] for i, j in pairs)/max(len(a), len(b))


def validate_submission(rows, provenance_sequences, controls):
    names, sequences = set(), set()
    for row in rows:
        name, seq = row['name'], row['sequence']
        assert name not in names and seq not in sequences, 'Duplicate name or sequence'
        assert name in provenance_sequences and seq == provenance_sequences[name], 'Sequence provenance mismatch'
        assert seq not in controls and SUFFIX not in seq, 'Control or assay suffix in submission'
        assert 10 <= len(seq) <= 250 and not set(seq)-set('ACDEFGHIKLMNPQRSTVWY')
        assert row['molecule_class'] == 'protein'
        names.add(name); sequences.add(seq)
    assert 1 <= len(rows) <= 20, 'Tracks 2/3 limit'


def summary_evidence(round_name, report, candidate):
    models = candidate['selected_models']
    tag = round_name == 'tag_comparison'
    return {
        'source': f'reports/{round_name}/results.json', 'run': report['run'],
        'condition': candidate.get('condition', 'bare'),
        'selection_rule': report.get('selection_rule', 'Frozen core mean cross-PAE rule; see tag-comparison manifest'),
        'passing_samples': candidate['passing_core_samples' if tag else 'passing_samples'],
        'total_samples': sum(r['id'] == candidate['id'] for r in report['rows']),
        'both_selected_core_or_full_screen': candidate['core_pass_both_seeds' if tag else 'all_selected_seeds_full_screen'],
        'both_selected_tag_independent': candidate.get('tag_independent_pass_both_seeds') if candidate.get('condition') == 'tagged' else None,
        'selected_models': [{k: r[k] for k in ('seed', 'sample', 'metrics', 'scores', 'screen', 'pose', 'geometry')}
                            | {'core_mean_cross_PAE_A': r.get('core_mean_cross_PAE_A')} for r in models],
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prior_manifest = json.loads((OUT/'manifest.json').read_text()) if (OUT/'manifest.json').exists() else {}
    archived = {m['path']: m for m in prior_manifest.get('models', [])}
    sources = {f'reports/{r}/results.json': json.loads((ROOT/f'reports/{r}/results.json').read_text()) for r in ROUNDS}
    tag_manifest = json.loads((ROOT/'data/tag_comparison/manifest.json').read_text())
    direction = json.loads((ROOT/'reports/tag_comparison/direction.json').read_text())
    novelty_path = ROOT/'reports/novelty_review/summary.json'
    novelty = json.loads(novelty_path.read_text()) if novelty_path.exists() else None
    seqs, evidence, latest_bare = {}, {}, {}
    for r in ROUNDS:
        report = sources[f'reports/{r}/results.json']
        for c in report['candidate_summaries']:
            name = c.get('base_id', c['id'])
            if 'sequence' in c:
                assert name not in seqs or seqs[name] == c['sequence']
                seqs[name] = c['sequence']
            evidence.setdefault(name, []).append(summary_evidence(r, report, c))
            if c.get('condition', 'bare') == 'bare':
                latest_bare[name] = c['selected_models']
    for c in tag_manifest['cases']:
        if c['kind'] == 'novel_design':
            assert seqs[c['base_id']] == c['core_sequence']
    controls = {c['id']: c['core_sequence'] for c in tag_manifest['cases'] if c['condition']=='bare' and c['kind']!='novel_design'}
    rows = [{'name': name, 'sequence': seqs[name], 'molecule_class': 'protein'} for name, *_ in SHORTLIST]
    validate_submission(rows, seqs, set(controls.values()))
    with (OUT/'submission_draft.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['name', 'sequence', 'molecule_class'])
        writer.writeheader(); writer.writerows(rows)
    (OUT/'submission_draft.fasta').write_text(''.join(f">{r['name']}\n{r['sequence']}\n" for r in rows))
    candidates, model_index, structures = [], [], {}
    for rank, (name, family, role, rationale) in enumerate(SHORTLIST, 1):
        copied = []
        for model in latest_bare[name]:
            source = ROOT/model['raw_scoring_output']/'model.pdb'
            dest = OUT/'models'/f"{name}_seed{model['seed']}_sample{model['sample']}.pdb"
            dest.parent.mkdir(exist_ok=True)
            if source.exists():
                shutil.copyfile(source, dest)
            else:
                # A fresh clone includes the selected model, but not ignored raw outputs.
                record = archived.get(str(dest.relative_to(ROOT)))
                assert record and record['source'] == str(source.relative_to(ROOT)), 'Missing original or archived model provenance'
                assert dest.exists() and digest(dest) == record['sha256'], 'Archived selected model hash mismatch'
            structure = read_structure(dest)
            assert sequence(structure[0]['B']) == seqs[name]
            assert sequence(structure[0]['A']) == tag_manifest['target_sequence']
            item = {'name': name, 'seed': model['seed'], 'sample': model['sample'],
                    'source': str(source.relative_to(ROOT)), 'path': str(dest.relative_to(ROOT)),
                    'sha256': digest(dest), 'selected_screen_pass': model['screen']['pass'],
                    'selected_screen_failures': model['screen']['failed']}
            copied.append(item); model_index.append(item)
            structures.setdefault(name, structure)
        angles = [(math.degrees(phi), math.degrees(psi)) for pp in PPBuilder().build_peptides(structures[name][0]['B'])
                  for phi, psi in pp.get_phi_psi_list() if phi is not None and psi is not None]
        candidates.append({'rank': rank, 'name': name, 'family': family, 'role': role,
                           'ranking_rationale': rationale, 'sequence': seqs[name],
                           'sequence_sha256': hashlib.sha256(seqs[name].encode()).hexdigest(),
                           'sequence_audit': sequence_audit(seqs[name]), 'evidence_by_round': evidence[name],
                           'latest_direction': [v for v in direction['variants'] if v['base_id'] == name],
                           'latest_selected_bare_models': copied,
                           'backbone_angle_description_first_selected_model': {
                               'angle_count': len(angles),
                               'alpha_bin_fraction': sum(-100<=p<=-30 and -80<=q<=-5 for p,q in angles)/len(angles),
                               'beta_bin_fraction': sum(-180<=p<=-80 and (q>=80 or q<=-120) for p,q in angles)/len(angles)},
                           'experimental_binding': 'unmeasured', 'experimental_pH_selectivity': 'unmeasured',
                           'mouse_cross_reactivity': 'unknown', 'tag_compatibility': 'unresolved'})
    comparisons = []
    for a, b in itertools.combinations(candidates, 2):
        pairs, identity = alignment(a['sequence'], b['sequence'])
        ca = [r['CA'] for r in structures[a['name']][0]['B'] if r.id[0]==' ']
        cb = [r['CA'] for r in structures[b['name']][0]['B'] if r.id[0]==' ']
        sup = Superimposer(); sup.set_atoms([ca[i] for i,j in pairs], [cb[j] for i,j in pairs])
        comparisons.append({'a': a['name'], 'b': b['name'], 'same_parent_lineage': a['family']==b['family'],
                            'global_sequence_identity': identity, 'aligned_residue_pairs': len(pairs),
                            'coverage_a': len(pairs)/len(ca), 'coverage_b': len(pairs)/len(cb),
                            'sequence_aligned_binder_CA_RMSD_A': float(sup.rms)})
    known = [{'candidate': c['name'], 'control': name, 'global_sequence_identity': alignment(c['sequence'], seq)[1]}
             for c in candidates for name, seq in controls.items()]
    result = {
        'status': 'DRAFT_NOT_SUBMITTED', 'as_of_date': '2026-10-21', 'track_scope': 'User confirmed Track 2 or 3',
        'ranking_policy': 'Editorial prioritization: pH hypothesis first, repeated human structural retention next, additional lineages then matched parent. Not a calibrated probability or automatic score.',
        'candidate_count': len(rows), 'parent_lineage_count': len({c['family'] for c in candidates}),
        'independent_fold_class_count': None, 'candidates': candidates,
        'reserve_or_excluded': [{'name': n, 'reason': reason, 'sequence': seqs[n], 'evidence_by_round': evidence[n]} for n, reason in RESERVE.items()],
        'additional_explicit_exclusion': {'name': 'apical_long_015_D51H_T70H', 'reason': 'Independent prediction lost interface; historical fixed-geometry pKa results do not qualify it.', 'source': 'reports/mutant_validation/RESULTS.md'},
        'pairwise_diversity': comparisons, 'local_control_identity_comparison': known,
        'diversity_method': 'Global alignment: match 2, mismatch -1, gap open -4, gap extension -0.5; matches divided by longer full sequence. RMSD is binder-fitted on sequence-aligned CA pairs, not receptor-frame interface RMSD, and is not a structure-search or fold-novelty test. First selected model in latest bare round is used, including failures.',
        'novelty_scope': 'Database review is linked below when available. Five parent lineages are not five established fold classes; no exhaustive similarity search or global novelty proof. All share one BoltzGen apical campaign.',
        'database_similarity_review': {'source':'reports/novelty_review/summary.json', 'sha256':digest(novelty_path), 'status':novelty['status']} if novelty else None,
        'software': {'python': platform.python_version(), 'biopython': Bio.__version__},
        'rules_source': 'https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1',
        'sequence_descriptor_source': 'https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html',
        'new_modal_compute_usd': 0,
    }
    (OUT/'audit.json').write_text(json.dumps(result, indent=2)+'\n')
    dependencies = list(sources) + ['data/tag_comparison/manifest.json', 'reports/tag_comparison/direction.json',
                                   'reports/receptor_preparation_v2/RESULTS.md', 'reports/expansion_his_direction/RESULTS.md',
                                   'reports/tag_comparison/billing_after.json',
                                   'scripts/build_submission_draft.py']
    if novelty: dependencies.append('reports/novelty_review/summary.json')
    manifest = {'status': 'DRAFT_NOT_SUBMITTED', 'inputs': {p: digest(ROOT/p) for p in dependencies},
                'models': model_index, 'outputs': {p.name: digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}}
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Built {len(rows)} draft sequences across {result["parent_lineage_count"]} parent lineages; {len(model_index)} selected models, including failures.')
    for c in candidates:
        a=c['sequence_audit']; print(c['rank'],c['name'],a['length'],'pI',round(a['estimated_pI'],2),'GRAVY',round(a['GRAVY_Kyte_Doolittle'],2),'sequons',a['N_X_S_or_T_sequons'],'runs',a['identical_residue_runs_ge4'])


if __name__ == '__main__':
    main()
