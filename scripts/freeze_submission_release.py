"""Create a portable, hash-verified review package from a committed draft."""
import hashlib,json,re,shutil,subprocess,zipfile
from datetime import datetime,timezone
from pathlib import Path
from novelty_search import ROOT,sha,save


def portable_markdown(source,commit):
    def replace(match):
        label,target=match.groups()
        if target.startswith(('https://','http://','#','mailto:')):return match.group(0)
        path=(source.parent/target).resolve()
        assert path.is_relative_to(ROOT) and path.exists(),target
        kind='tree' if path.is_dir() else 'blob'
        return f'[{label}](https://github.com/anushajain2/adaptyv-challenge4/{kind}/{commit}/{path.relative_to(ROOT)})'
    return re.sub(r'\[([^\]]+)\]\(([^)]+)\)',replace,source.read_text())


def main():
    draft=ROOT/'reports/submission_draft';out=ROOT/'reports/submission_release_2026-10-25'
    assert not out.exists(),'Release already frozen; do not silently overwrite.'
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert not subprocess.check_output(['git','status','--porcelain','--','reports/submission_draft','reports/final_his_panel','reports/final_his_direction','scripts/freeze_submission_release.py'],cwd=ROOT,text=True).strip(),'Commit reviewed sources first.'
    audit=json.loads((draft/'audit.json').read_text());manifest=json.loads((draft/'manifest.json').read_text())
    for p,h in manifest['inputs'].items():assert sha(ROOT/p)==h,p
    for p,h in manifest['outputs'].items():assert sha(draft/p)==h,p
    out.mkdir();models=out/'models';models.mkdir()
    for name in ['submission_draft.csv','submission_draft.fasta','audit.json']:
        shutil.copyfile(draft/name,out/name)
    for name in ['README.md','METHODS.md']:
        (out/('REVIEW.md' if name=='README.md' else name)).write_text(portable_markdown(draft/name,commit))
    for m in manifest['models']:
        source=ROOT/m['path'];assert sha(source)==m['sha256'];shutil.copyfile(source,models/source.name)
    panel=ROOT/'reports/final_his_panel/RESULTS.md'
    (out/'FINAL_EXPERIMENT.md').write_text(portable_markdown(panel,commit))
    review=f'''# TfR1 submission package — frozen for participant review

This package contains {audit['candidate_count']} ranked bare sequences across {audit['parent_lineage_count']} parent lineages. It has **not been submitted**. Computational consistency does not establish binding, pH switching or mouse cross-reactivity.

Upload `submission_draft.csv` after reviewing the sequences, ranking and limitations in `REVIEW.md`. `METHODS.md` and `FINAL_EXPERIMENT.md` describe the methods and final bounded experiment. `audit.json` preserves per-round evidence, including failed selected models; `models/` contains the {len(manifest['models'])} currently selected bare models. The accompanying FASTA is for review, not the required CSV upload format.

The CSV meets the shared Tracks 2/3 size/format limits. Track 3 is recommended for participants outside the accepted Track 2 program; support eligibility is not inferred from these files. No Claude usage is claimed for this Codex-driven work.

The challenge closes November 1, 2026, 23:59 AoE (November 2, 17:29 IST). Review the [official challenge](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1) and [track requirements](https://proteinbase.com/competitions/anthropic-adaptyv-2026) before submitting. The submitting participant must review the designs; this package does not accept terms or submit on their behalf.

Source snapshot: [GitHub commit {commit[:7]}](https://github.com/anushajain2/adaptyv-challenge4/tree/{commit}). Detailed report links are pinned to this commit so this folder can be moved without breaking the evidence references. Original paths in `audit.json` are repository-relative.
'''
    (out/'README.md').write_text(review)
    save(out/'manifest.json',{'status':'FROZEN_FOR_REVIEW_NOT_SUBMITTED','created_utc':datetime.now(timezone.utc).isoformat(),'source_commit':commit,'source_draft_manifest_sha256':sha(draft/'manifest.json'),'candidate_count':audit['candidate_count'],'parent_lineages':audit['parent_lineage_count'],'files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}})
    archive=ROOT/'outputs/submission_release_2026-10-25.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file():z.write(p,str(p.relative_to(out.parent)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for rel,digest in json.loads((out/'manifest.json').read_text())['files'].items():
            assert hashlib.sha256(z.read(out.name+'/'+rel)).hexdigest()==digest
    save(ROOT/'reports/submission_release_2026-10-25_archive.json',{'archive':str(archive.relative_to(ROOT)),'sha256':sha(archive),'bytes':archive.stat().st_size,'source_commit':commit,'verified_all_member_hashes':True})
    print('Frozen review package:',out,'; ZIP:',archive)

if __name__=='__main__':main()
