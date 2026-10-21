import ast
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from novelty_search import ROOT, sha
from launch_diversity_round import PHASES, PRIOR_ALLOWANCE, resource_bound


def test_new_round_preserves_audited_target_hotspots_and_raw_gates():
    old=json.loads((ROOT/'data/apical_expansion_fixed/manifest.json').read_text())
    new=json.loads((ROOT/'data/diversity_generation/manifest.json').read_text())
    assert [a['num_designs'] for a in new['arms']]==[60,40]
    assert new['target_sequence']==old['target_sequence']
    assert new['selection_rule']==old['selection_rule']
    for key in old['screen']:
        if key!='decision':assert old['screen'][key]==new['screen'][key]
    for a,b in zip(old['arms'],new['arms']):
        assert a['hotspots_crop']==b['hotspots_crop']==[19,21,22]
        assert a['generation_seed']!=b['generation_seed'] and a['inverse_fold_seed']!=b['inverse_fold_seed']
        assert (ROOT/'data/apical_expansion_fixed'/f"{a['id']}.yaml").read_bytes()==(ROOT/'data/diversity_generation'/f"{b['id']}.yaml").read_bytes()
    for name,digest in new['input_sha256'].items():assert sha(ROOT/'data/diversity_generation'/name)==digest
    catalog=json.loads((ROOT/'data/diversity_generation/prior_catalog.json').read_text())
    assert len(catalog['prior_designs'])==495 and len(catalog['current_draft'])==6


def test_round_and_modal_hard_limits_fit_remaining_conservative_budget():
    assert sum(p[2] for p in PHASES.values())==4
    assert PRIOR_ALLOWANCE+4<30
    for phase,(_,folder,reserve,work,timeout,_) in PHASES.items():
        assert resource_bound(phase)<reserve and work<timeout
        tree=ast.parse((ROOT/'scripts'/f'modal_{folder}.py').read_text())
        calls=[d for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) for d in n.decorator_list
               if isinstance(d,ast.Call) and isinstance(d.func,ast.Attribute) and d.func.attr=='function']
        assert len(calls)==1
        opts={k.arg:ast.literal_eval(k.value) for k in calls[0].keywords if k.arg in {
            'timeout','startup_timeout','retries','max_containers','cpu','memory','gpu','scaledown_window'}}
        assert opts==dict(timeout=timeout,startup_timeout=120,retries=0,max_containers=1,
                         cpu=(2,4),memory=(16384,32768),gpu='L40S',scaledown_window=2)
