import ast,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_controls import ROOT,read_structure,sequence
from prepare_expansion_his import PARENTS,choose_confirmation
from launch_expansion_his import PHASES,BOUND


def test_new_panel_changes_only_chosen_sites_and_keeps_target_and_backbone():
    m=json.loads((ROOT/'data/expansion_his_panel/manifest.json').read_text());lookup={c['id']:c for c in m['cases']}
    assert len(lookup)==32 and sum(c.get('parent_id') in PARENTS for c in m['cases'])==24
    assert len({c['binder_sequence'] for c in m['cases']})==32
    for c in m['cases']:
        s=read_structure(ROOT/'data/expansion_his_panel'/f"{c['id']}.pdb")
        assert sequence(s[0]['A'])==m['target_sequence'] and sequence(s[0]['B'])==c['binder_sequence']
        with np.load(ROOT/'data/expansion_his_panel'/f"{c['id']}.npz",allow_pickle=False) as d:
            assert not d['design_mask'][:194].any() and d['design_mask'][194:].all()
        if c.get('parent_id') not in PARENTS:continue
        parent=read_structure(ROOT/'data/expansion_his_panel'/f"{c['parent_id']}.pdb")
        expected=sorted(int(x[1:-1]) for x in c['mutations'])
        assert [i+1 for i,(a,b) in enumerate(zip(lookup[c['parent_id']]['binder_sequence'],c['binder_sequence'])) if a!=b]==expected
        for chain in ('A','B'):
            for r in parent[0][chain]:
                for name in ('N','CA','C','O'):np.testing.assert_array_equal(r[name].coord,s[0][chain][r.id][name].coord)


def test_chosen_sites_and_pairs_obey_two_pose_anchor_exclusions():
    audits=json.loads((ROOT/'data/expansion_his_panel/contact_audit.json').read_text())
    for audit in audits:
        chosen=set(audit['chosen_sites']);assert len(chosen)==5
        for row in audit['site_audit']:
            if row['position'] in chosen:
                assert not row['excluded_reasons']
                for x in row['models']:
                    assert not x['salt_contact'] and not x['high_contact_anchor_proxy']
                    assert 2<=x['buried_sasa_A2']<=110 and x['free_sasa_A2']>=25 and x['min_target_distance_A']<=5
        assert len(audit['chosen_pairs'])==3
        for p in audit['chosen_pairs']:assert set(p['positions'])<=chosen and all(5<=x<=16 for x in p['CA_distances_A'])


def test_confirmation_never_selects_failed_or_comparison_variants():
    manifest={'cases':[]};rows=[]
    for p in PARENTS+['apical_short_008']:
        for i in range(5):
            name=f'{p}_H{i}';manifest['cases'].append({'id':name,'parent_id':p})
            rows.append({'id':name,'all_selected_seeds_candidate_screen':i!=4,'selected_models':[{'metrics':{'design_to_target_iptm':.5+i/20}}]})
    r={'controls_pass':True,'candidate_summaries':rows};chosen=choose_confirmation(r,manifest)
    assert len(chosen)==6 and all(x.endswith(('H2','H3')) for x in chosen) and not any(x.startswith('apical_short_008') for x in chosen)
    r['controls_pass']=False;assert choose_confirmation(r,manifest) is None


def test_cloud_resource_limits_fit_frozen_reservations():
    assert BOUND<2 and 7.96348782+12+4<30
    for run,folder in PHASES.values():
        tree=ast.parse((ROOT/'scripts'/f'modal_{folder}.py').read_text())
        calls=[d for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) for d in n.decorator_list if isinstance(d,ast.Call) and isinstance(d.func,ast.Attribute) and d.func.attr=='function']
        assert len(calls)==1
        opts={k.arg:ast.literal_eval(k.value) for k in calls[0].keywords if k.arg in {'timeout','startup_timeout','retries','max_containers','cpu','memory','gpu','scaledown_window'}}
        assert opts==dict(timeout=2250,startup_timeout=120,retries=0,max_containers=1,cpu=(2,4),memory=(16384,32768),gpu='L40S',scaledown_window=2)
