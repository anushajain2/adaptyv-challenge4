"""Freeze a local sequence redesign experiment before sampling new sequences."""
import copy
import json
import numpy as np
from novelty_search import ROOT, sha, save, now
from prepare_controls import read_structure, sequence
from prepare_expansion_his import site_metrics

PIN = '8907e6671bfbfc92303b5f79c4b5e6ce47cdef57'


def choose_sites(models, his_position):
    metrics = [site_metrics(s) for s in models]
    seq = sequence(models[0][0]['B']); rows = []
    assert all(sequence(s[0]['B']) == seq for s in models) and seq[his_position-1] == 'H'
    for position, aa in enumerate(seq, 1):
        distances = []
        for s in models:
            own = np.array([a.coord for a in s[0]['B'][position] if a.element != 'H'])
            his = np.array([a.coord for a in s[0]['B'][his_position] if a.name in {'CG','ND1','CD2','CE1','NE2'}])
            assert len(his) == 5
            distances.append(float(np.linalg.norm(own[:,None]-his, axis=2).min()))
        reasons = []
        if aa in 'GPCHDEKR': reasons.append('preserve charged, His, Gly, Pro and Cys')
        if position <= 5 or position > len(seq)-5: reasons.append('preserve terminal five')
        if any(m[position]['salt_contact'] or m[position]['high_contact_anchor_proxy'] for m in metrics): reasons.append('preserve anchor proxy in either pose')
        if max(distances) > 8: reasons.append('outside 8 A His shell in a pose')
        rows.append({'position': position, 'residue': aa, 'His_ring_distance_A': distances,
                     'excluded_reasons': reasons, 'site_metrics': [m[position] for m in metrics]})
    eligible = sorted([r for r in rows if not r['excluded_reasons']], key=lambda r: (max(r['His_ring_distance_A']), r['position']))[:8]
    assert len(eligible) >= 2, 'Insufficient local sequence freedom; do not sample.'
    return sorted(r['position'] for r in eligible), rows


def main():
    out = ROOT/'data/local_redesign'; out.mkdir(exist_ok=False)
    families = []
    specifications = [
        ('apical_short_008_A41H', 'apical_short_008_A41H_bare', 41, 'tag_comparison', [2003,2011], 'data/his_confirmation', 'apical_short_008'),
        ('diversity_long_037_S90H', 'diversity_long_037_S90H', 90, 'final_his_panel', [2503,2609], 'data/final_his_panel', 'diversity_long_037'),
    ]
    for name, source_id, his, report, seeds, inputs, unmodified in specifications:
        source = ROOT/'reports'/report/'results.json'; data = json.loads(source.read_text())
        selected = sorted([r for r in data['rows'] if r['id'] == source_id and r['confidence_selected']], key=lambda r:r['seed'])
        assert [r['seed'] for r in selected] == seeds and all(r['screen']['pass'] for r in selected)
        models = [read_structure(ROOT/r['raw_scoring_output']/'model.pdb') for r in selected]
        mutable, audit = choose_sites(models, his)
        files = []
        for r in selected:
            src = ROOT/r['raw_scoring_output']/'model.pdb'; dest = out/f'{name}_seed{r["seed"]}.pdb'
            dest.write_bytes(src.read_bytes()); files.append({'path':str(dest.relative_to(ROOT)), 'sha256':sha(dest), 'seed':r['seed'], 'sample':r['sample'], 'source':str(src.relative_to(ROOT))})
        families.append({'id':name,'histidine_position':his,'sequence':sequence(models[0][0]['B']),
                         'unmodified_parent':unmodified,'baseline_input_directory':inputs,'mutable_positions':mutable,
                         'source_results':str(source.relative_to(ROOT)),'source_results_sha256':sha(source),
                         'models':files,'site_audit':audit})
    plan = {'prepared_utc':now(),'status':'PROSPECTIVE_NOT_SUBMITTED','families':families,
        'proteinmpnn':{'repository':'https://github.com/dauparas/ProteinMPNN','commit':PIN,'model':'vanilla_model_weights/v_48_020.pt',
            'weights_sha256':sha(ROOT/'vendor/ProteinMPNN/vanilla_model_weights/v_48_020.pt'),
            'temperatures':[.2,.4],'samples_per_temperature_per_backbone':32,'sampling_seeds':[3109,3203],
            'omit_new_amino_acids':'CHX','backbone_noise':0,'score_decoding_seeds':[3301,3307,3313,3319],
            'ranking':'Mean mutable-position negative log likelihood across both source backbones and four common decoding seeds. This is sequence compatibility, not binding or pH energy.'},
        'selection':{'maximum_per_family':6,'minimum_pairwise_hamming':2,'additional_mutations_range':[1,4],
            'rule':'Deduplicate; retain all fixed positions exactly; reject new N-X-S/T sequons; rank by ensemble NLL, then sequence lexically; greedily require Hamming >=2 among chosen designs. No replenishment after prediction.'},
        'refolding':{'seeds':[2909,3001],'samples_per_seed':5,'target':'6WRW receptor only; empty target MSA; no binder template',
            'selection':'Unchanged target integrity filter then 0.8 ipTM + 0.2 pTM; both selected poses required; unchanged raw geometry/confidence screen.'},
        'preparation':'Unchanged v2 numerical preparation applied once to all selected models, including both His baselines, unmodified parents and three controls. Retain strict historical movement check.',
        'sensitivity':'Separately report chain backbone RMS <=0.5 A, every CA movement <=1 A and backbone movement <=1 A for raw residues within 6 A of the other chain. Do not use this sensitivity result to override strict failures or automatically promote designs.',
        'direction':'Same per-His bound-free pKa <=-0.5, free pKa 6-7, His-only ratio >=2, artificial target termini >=10 A in both poses; valid PROPKA groups, strict geometry and controls required. Report matched His baseline and unmodified parent separately. Report His ring charge/donor/acceptor distances; proximity alone is not an H-bond or binding energy.',
        'followup':'At most two strict two-pose passes ranked by worst-pose His-only ratio then worst-pose ipTM then ID. Require ratio improvement >=20% over matched His baseline in both paired seeds (or baseline direction failure with both baseline structures eligible). Only then run fresh matched bare/tag comparisons with controls. Otherwise preserve fallback and report limitations.',
        'budget':{'local_mpnn_usd':0,'panel_max_usd':2,'conditional_followup_max_usd':3,'unallocated_contingency_usd':1,'total_new_cap_usd':6,'project_gross_cap_usd':30},
        'limitations':['No experimental binding, release or mouse claims','ProteinMPNN is not pH aware and does not place His side chains','Local redesign adds sequence compatibility search, not a new backbone or new lineage','Prediction seeds are not experimental replicates','Frozen eight-design release remains a fallback; no submission authorized']}
    save(out/'plan.json',plan);save(ROOT/'reports/local_redesign/plan.json',plan)
    print(json.dumps([{'id':f['id'],'mutable_positions':f['mutable_positions']} for f in families],indent=2))


if __name__ == '__main__': main()
