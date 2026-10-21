"""Verify saved launch hashes, selected review models and complete preparation records."""
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
from prepare_controls import ROOT
from launch_generation_comparison import gross_cost


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def main():
    reports=ROOT/'reports';out=reports/'next_round';checks={}
    for folder in ['his_panel','his_confirmation','apical_expansion_fixed','expansion_validation','species_exploration']:
        record=reports/folder/'run_record';launch=read(record/'launch_ledger.json');manifest=read(record/'manifest.json')
        assert read(record/'input_hashes.json')==launch['input_sha256'],folder
        for file,sha in manifest['input_sha256'].items():assert digest(ROOT/'data'/folder/file)==sha,(folder,file)
        for file,sha in launch['analysis_sha256'].items():assert digest(record/file)==sha,(folder,file)
        status=read(record/'run_status.json');assert all(s['returncode']==0 for s in status['statuses']),folder
        checks[folder]={'cloud_launch_hashes_match':True,'local_input_hashes_match':True,'archived_analysis_hashes_match':True,'successful_worker_status':True}
    prep=[]
    for folder,expected in [('structural_preparation',12),('his_panel',25),('his_confirmation',20),('apical_expansion_fixed',203),('expansion_validation',18)]:
        path=reports/folder/('results.json' if folder=='structural_preparation' else 'preparation/results.json');d=read(path)
        assert len(d['rows'])==expected,(folder,len(d['rows']),expected)
        if folder!='structural_preparation':assert d['source_results_sha256']==digest(reports/folder/'results.json')
        exclusions=[]
        for r in d['rows']:
            assert digest(ROOT/r['source'])==r['source_sha256']
            if r['status']=='ok':
                dest=ROOT/r['dest'] if 'dest' in r else reports/folder/'models'/f"{r['id']}_seed{r['seed']}"
                assert digest(dest/'prepared.pdb')==r['prepared_sha256']
            if r['status']!='ok' or not r.get('displacement',{}).get('pass'):
                exclusions.append({'id':r['id'],'seed':r['seed'],'status':r['status'],'error':r.get('error'),'displacement':r.get('displacement')})
        prep.append({'experiment':folder,'expected':expected,'completed':len(d['rows']),'execution_ok':sum(r['status']=='ok' for r in d['rows']),'movement_pass':sum(r.get('displacement',{}).get('pass',False) for r in d['rows']),'excluded_prepared_models':exclusions})
    bundle=read(out/'review_bundle/manifest.json');assert len(bundle['models'])==10
    for model in bundle['models']:
        for name,sha in model['hashes'].items():assert digest(ROOT/model['path']/name)==sha
    assert len(read(reports/'species_exploration/results.json')['rows'])==81
    apps=read(out/'apps_after.json');assert all(a['state']=='stopped' for a in apps)
    billing=read(out/'billing_after.json');gross=max(gross_cost(billing),gross_cost(read(out/'billing_first_after.json')));assert gross<30
    record={'audited_utc':datetime.now(timezone.utc).isoformat(),'run_checks':checks,'preparation':prep,'review_models_verified':10,'all_cloud_apps_stopped':True,'gross_meter_high_water_usd':gross,'latest_gross_meter_usd':gross_cost(billing),'authorized_gross_cap_usd':30,'billing_may_lag':True}
    (out/'evidence_audit.json').write_text(json.dumps(record,indent=2))
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
