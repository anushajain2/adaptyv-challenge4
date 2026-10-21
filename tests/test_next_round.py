import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_controls import ROOT,read_structure,sequence
from prepare_next_round import mutate_reference
from launch_next_round import PHASES,resource_bound
from structural_preparation import displacement


def test_panel_changes_only_requested_sites_and_excludes_binder_template():
    m=json.loads((ROOT/'data/his_panel/manifest.json').read_text())
    assert len(m['cases'])==25 and len({c['binder_sequence'] for c in m['cases']})==25
    for c in m['cases']:
        if 'parent_id' not in c:continue
        parent=read_structure(ROOT/'data/confirmation'/f"{c['parent_id']}.pdb")
        ps=[int(s[1:-1]) for s in c['mutations']]
        s=mutate_reference(parent,ps)
        assert sequence(s[0]['B'])==c['binder_sequence']
        assert [i+1 for i,(a,b) in enumerate(zip(sequence(parent[0]['B']),c['binder_sequence'])) if a!=b]==ps
        for chain in ('A','B'):
            for r in parent[0][chain]:
                for a in ('N','CA','C','O'):np.testing.assert_array_equal(r[a].coord,s[0][chain][r.id][a].coord)
        with np.load(ROOT/'data/his_panel'/f"{c['id']}.npz") as d:
            assert not d['design_mask'][:194].any() and d['design_mask'][194:].all()


def test_reserved_runtime_bounds_and_total_allocation():
    repair=json.loads((ROOT/'reports/apical_expansion_fixed/repair.json').read_text())
    assert abs(sum(v[2] for v in PHASES.values())+repair['failed_attempt_allowance_usd']-12)<1e-9
    assert 7.96348782+12<30
    for phase,v in PHASES.items():assert resource_bound(phase)<v[2]


def test_preparation_rejects_pose_motion_not_removed_by_alignment():
    s=read_structure(ROOT/'data/confirmation/apical_long_015.pdb')
    import copy
    p=copy.deepcopy(s)
    assert displacement(s,p)['pass']
    for a in p[0]['B'].get_atoms():a.coord+=np.array([2.,0,0])
    assert not displacement(s,p)['pass']


def test_expansion_seeds_and_ids_are_distinct_and_apical_only():
    old=json.loads((ROOT/'data/broad_campaign/manifest.json').read_text())
    new=json.loads((ROOT/'data/apical_expansion/manifest.json').read_text())
    assert sum(a['num_designs'] for a in new['arms'])==200
    assert all(a['epitope']=='apical' and a['protocol']=='protein-anything' for a in new['arms'])
    for field in ('generation_seed','inverse_fold_seed','id'):
        assert not {a[field] for a in old['arms']} & {a[field] for a in new['arms']}


def test_modal_decorators_match_reserved_resource_bounds():
    import ast
    for phase,(_,folder,_,_,timeout,startup) in PHASES.items():
        tree=ast.parse((ROOT/'scripts'/f'modal_{folder}.py').read_text())
        calls=[d for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) for d in n.decorator_list if isinstance(d,ast.Call) and isinstance(d.func,ast.Attribute) and d.func.attr=='function']
        assert len(calls)==1
        opts={k.arg:ast.literal_eval(k.value) for k in calls[0].keywords if k.arg in {'timeout','startup_timeout','retries','max_containers','cpu','memory','gpu','scaledown_window'}}
        assert opts==dict(timeout=timeout,startup_timeout=startup,retries=0,max_containers=1,cpu=(2,4),memory=(16384,32768),gpu='L40S',scaledown_window=2)


def test_generation_count_includes_diffusion_samples():
    from apical_expansion_worker import design_count
    assert design_count(10,10)==100
    assert design_count(40,1)==40


def test_confirmation_shortlist_is_bounded_and_never_uses_failed_variants():
    from prepare_his_confirmation import choose
    m={'cases':[]};rows=[]
    for parent in ('apical_long_015','apical_short_008'):
        for i in range(6):
            name=f'{parent}_v{i}'
            m['cases'].append({'id':name,'parent_id':parent})
            rows.append({'id':name,'all_selected_seeds_candidate_screen':i!=5,'selected_models':[{'metrics':{'design_to_target_iptm':i/10}}]})
    result={'controls_pass':True,'candidate_summaries':rows}
    chosen=choose(result,m)
    assert len(chosen)==6 and all(not x.endswith(('v0','v1','v5')) for x in chosen)
    result['controls_pass']=False
    assert choose(result,m) is None


def test_species_inputs_use_actual_target_sequences_and_same_binders():
    from Bio import SeqIO
    m=json.loads((ROOT/'data/species_exploration/manifest.json').read_text())
    assert len(m['cases'])==27 and m['refolding']['diffusion_samples_per_seed']==3
    mouse=str(next(SeqIO.parse(ROOT/'data/raw/Q62351.fasta','fasta')).seq)[191:385]
    human=str(next(SeqIO.parse(ROOT/'data/raw/P02786.fasta','fasta')).seq)[189:383]
    assert len(mouse)==len(human)==194 and mouse!=human
    groups={}
    for c in m['cases']:
        s=read_structure(ROOT/'data/species_exploration'/f"{c['id']}.pdb")
        assert sequence(s[0]['A'])==c['target_sequence']==(mouse if c['condition']=='maf' else human)
        assert sequence(s[0]['B'])==c['binder_sequence']
        groups.setdefault(c['base_id'],set()).add(c['binder_sequence'])
        with np.load(ROOT/'data/species_exploration'/f"{c['id']}.npz") as d:
            assert not d['design_mask'][:194].any() and d['design_mask'][194:].all()
    assert len(groups)==9 and all(len(seqs)==1 for seqs in groups.values())
