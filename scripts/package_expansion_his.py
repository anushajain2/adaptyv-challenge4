"""Package passing structural variants and all direction outcomes without promoting failures."""
import hashlib,json,shutil
from pathlib import Path
from prepare_controls import ROOT
from prepare_expansion_his import PARENTS


def main():
    out=ROOT/'reports/expansion_his_direction';d=json.loads((out/'results.json').read_text());confirmation=json.loads((ROOT/'reports/expansion_his_confirmation/results.json').read_text());m=json.loads((ROOT/'data/expansion_his_confirmation/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    bundle=out/'review_bundle';bundle.mkdir(exist_ok=False);records=[]
    for c in confirmation['candidate_summaries']:
        case=cases[c['id']]
        if not case.get('mutations') or not c['all_selected_seeds_full_screen']:continue
        evidence=next(x for x in d['variants'] if x['id']==c['id']);models=[]
        for r in c['selected_models']:
            folder=bundle/f"{c['id']}_seed{r['seed']}";folder.mkdir();files={}
            raw=ROOT/r['raw_scoring_output']/'model.pdb';shutil.copyfile(raw,folder/'raw_complex.pdb');files['raw_complex.pdb']=hashlib.sha256(raw.read_bytes()).hexdigest()
            evidence_dir=out/'models'/folder.name
            if evidence_dir.exists():
                for f in evidence_dir.iterdir():
                    if f.is_file():shutil.copyfile(f,folder/f.name);files[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
            (folder/'raw_checks.json').write_text(json.dumps(r,indent=2));models.append({'seed':r['seed'],'files':files,'path':str(folder.relative_to(ROOT))})
        records.append({'id':c['id'],'parent_id':case['parent_id'],'sequence':c['sequence'],'new_parent_variant':case['parent_id'] in PARENTS,'comparison_only':case.get('comparison_only',False),'matched_parent_passes_both_seeds':evidence['matched_parent_passes_both_seeds'],'strict_direction_pass':evidence['direction_pass_both_seeds'],'models':models})
    parents=[]
    parent_ids={r['parent_id'] for r in records}
    for c in confirmation['candidate_summaries']:
        if c['id'] not in parent_ids:continue
        models=[]
        for r in c['selected_models']:
            folder=bundle/f"{c['id']}_seed{r['seed']}";folder.mkdir();raw=ROOT/r['raw_scoring_output']/'model.pdb';shutil.copyfile(raw,folder/'raw_complex.pdb')
            files={'raw_complex.pdb':hashlib.sha256(raw.read_bytes()).hexdigest()}
            prepared=ROOT/'reports/expansion_his_confirmation/preparation/models'/folder.name
            for name in ['prepared.pdb','preparation.json']:
                if (prepared/name).exists():shutil.copyfile(prepared/name,folder/name);files[name]=hashlib.sha256((prepared/name).read_bytes()).hexdigest()
            (folder/'raw_checks.json').write_text(json.dumps(r,indent=2));models.append({'seed':r['seed'],'path':str(folder.relative_to(ROOT)),'files':files})
        parents.append({'id':c['id'],'sequence':c['sequence'],'both_selected_raw_checks_pass':c['all_selected_seeds_full_screen'],'models':models})
    (bundle/'manifest.json').write_text(json.dumps({'variants':records,'parents':parents,'experimental_binding_and_switching_unvalidated':True,'all_failed_direction_evidence_preserved':True},indent=2))
    (bundle/'structurally_retained_variants.fasta').write_text(''.join(f">{r['id']} COMPUTATIONAL_ONLY strict_direction_pass={r['strict_direction_pass']}\n{r['sequence']}\n" for r in records))
    print('Packaged',len(records),'structurally retained variants with explicit direction outcomes')

if __name__=='__main__':main()
