"""Verify independently predicted His evidence against saved atom and pKa tables."""
import hashlib
import json
import numpy as np
from prepare_controls import ROOT,read_structure,sequence
from direction_check import HIS_ATOMS,histidines,release_ratio


def main():
    out=ROOT/'reports/expansion_his_direction';r=json.loads((out/'results.json').read_text());checks=[]
    m=json.loads((ROOT/'data/expansion_his_confirmation/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    for v in r['variants']:
        for model in v['models']:
            if model['status']!='ok':
                assert not model['direction_pass']
                checks.append({'id':v['id'],'seed':model['seed'],'status':'not_scored','failure_preserved':True})
                continue
            folder=out/'models'/f"{v['id']}_seed{model['seed']}"
            s=read_structure(folder/'complex.pdb');free=read_structure(folder/'free_binder.pdb')
            assert sequence(s[0]['A'])==m['target_sequence']
            assert sequence(s[0]['B'])==sequence(free[0]['B'])==cases[v['id']]['binder_sequence']
            assert hashlib.sha256((folder/'complex.pdb').read_bytes()).hexdigest()==model['source_sha256']
            for residue in s[0]['B']:
                for atom in residue:
                    if atom.element!='H':np.testing.assert_array_equal(atom.coord,free[0]['B'][residue.id][atom.name].coord)
            for pos in model['positions']:assert s[0]['B'][pos].resname=='HIS' and HIS_ATOMS<=set(s[0]['B'][pos].child_dict)
            values={}
            for key,filename in [('free','free_binder'),('bound','complex')]:
                assert json.loads((folder/f'{filename}.propka_warnings.json').read_text())==[]
                values[key]=histidines(json.loads((folder/f'{filename}.pka.json').read_text()),'B',model['positions'])
                np.testing.assert_allclose(values[key],model[f'{key}_pkas'],rtol=0,atol=1e-12)
            shifts=np.array(values['bound'])-values['free'];np.testing.assert_allclose(shifts,model['shifts'],atol=1e-12)
            ratio=release_ratio(values['free'],values['bound']);assert abs(ratio-model['his_only_Kd6_over_Kd7_4'])<1e-12
            decision=all(d<=-.5 for d in shifts) and any(6<=x<=7 for x in values['free']) and ratio>=2 and model['artificial_terminus_distance_A']>=10
            assert decision==model['direction_pass']
            checks.append({'id':v['id'],'seed':model['seed'],'saved_tables_and_atoms_match':True,'strict_rule_recomputed':decision})
        assert v['direction_pass_both_seeds']==all(x['direction_pass'] for x in v['models'])
    assert len(checks)==2*len(r['variants'])
    (out/'evidence_audit.json').write_text(json.dumps({'models':checks,'result':'All scored saved models and pKa tables match; no warnings in scored models; strict and preparation failures preserved'},indent=2))
    print('Verified',len(checks),'direction-model records')
if __name__=='__main__':main()
