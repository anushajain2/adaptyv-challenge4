"""Audit complete sample inventories, frozen selection, preparation and direction evidence."""
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from prepare_controls import ROOT
from prepare_expansion_his import choose_confirmation


def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    checks=[]
    for phase in ['panel','confirmation']:
        folder='expansion_his_'+phase;out=ROOT/'reports'/folder;m=load(ROOT/'data'/folder/'manifest.json');d=load(out/'results.json');record=out/'run_record';launch=load(record/'launch_ledger.json')
        assert load(record/'input_hashes.json')==launch['input_sha256']
        for name,value in m['input_sha256'].items():assert sha(ROOT/'data'/folder/name)==value,(folder,name)
        for name,value in launch['analysis_sha256'].items():assert sha(record/name)==value,(folder,name)
        assert sha(ROOT/'scripts/analyze_his_panel.py')==launch['analysis_sha256']['analyze_his_panel.py']
        expected={(c['id'],seed,i) for c in m['cases'] for seed in m['refolding']['seeds'] for i in range(m['refolding']['diffusion_samples_per_seed'])}
        actual={(r['id'],r['seed'],r['sample']) for r in d['rows']};assert expected==actual and len(actual)==len(d['rows'])
        prep=load(out/'preparation/results.json');selected={(r['id'],r['seed'],r['sample']) for r in d['rows'] if r['confidence_selected']}
        assert selected=={(r['id'],r['seed'],r['sample']) for r in prep['rows']} and len(prep['rows'])==len(selected)
        assert prep['source_results_sha256']==sha(out/'results.json')
        for r in prep['rows']:
            assert sha(ROOT/r['source'])==r['source_sha256']
            if r['status']=='ok':assert sha(ROOT/r['dest']/'prepared.pdb')==r['prepared_sha256']
        checks.append({'phase':phase,'predictions':len(actual),'selected_models':len(selected),'prepared':len(prep['rows']),'preparation_execution_ok':sum(r['status']=='ok' for r in prep['rows']),'preparation_movement_pass':sum(r.get('displacement',{}).get('pass',False) for r in prep['rows']),'controls_pass':d['controls_pass'],'all_sample_and_source_checks_pass':True})
    panel=load(ROOT/'reports/expansion_his_panel/results.json');pm=load(ROOT/'data/expansion_his_panel/manifest.json');cm=load(ROOT/'data/expansion_his_confirmation/manifest.json')
    assert choose_confirmation(panel,pm)==cm['selection_details']['chosen_variants']
    out=ROOT/'reports/expansion_his_direction';direction=load(out/'results.json');assert direction['source_results_sha256']==sha(ROOT/'reports/expansion_his_confirmation/results.json')
    launch=load(ROOT/'reports/expansion_his_confirmation/run_record/launch_ledger.json')
    assert sha(ROOT/'scripts/check_expansion_his_direction.py')==launch['analysis_sha256']['check_expansion_his_direction.py']
    audit=load(out/'evidence_audit.json');assert len(audit['models'])==sum(len(v['models']) for v in direction['variants'])
    # Also verify files from failed calculations; retaining a failure is not a reason to skip provenance.
    preparation=load(ROOT/'reports/expansion_his_confirmation/preparation/results.json')
    prep={(r['id'],r['seed']):r for r in preparation['rows']}
    for v in direction['variants']:
        for model in v['models']:
            folder=out/'models'/f"{v['id']}_seed{model['seed']}"
            assert sha(folder/'complex.pdb')==prep[(v['id'],model['seed'])]['prepared_sha256']
            if model['status']=='failed':
                warnings=load(folder/'complex.propka_warnings.json')
                assert any('Missing atoms or failed protonation' in w for w in warnings)
                assert not model['direction_pass']

    bundle=load(out/'review_bundle/manifest.json')
    for v in bundle['variants']+bundle['parents']:
        for m in v['models']:
            for name,value in m['files'].items():assert sha(ROOT/m['path']/name)==value
    record={'audited_utc':datetime.now(timezone.utc).isoformat(),'phases':checks,'frozen_confirmation_selection_matches':True,'direction_implementation_matches_confirmation_snapshot':True,'strict_direction_new_variant_passes':sum(v['direction_pass_both_seeds'] and not v['comparison_only'] for v in direction['variants']),'comparison_direction_passes':sum(v['direction_pass_both_seeds'] and v['comparison_only'] for v in direction['variants']),'review_bundle_hashes_match':True}
    (out/'round_evidence_audit.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))

if __name__=='__main__':main()
