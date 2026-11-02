"""Adapt the prepared package to the fields shown in the user's submission screenshot."""
import csv
import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'reports/submission_2026-11-02_recommended'
OUT=ROOT/'reports/submission_upload_2026-11-02'
NODE=Path('/Users/anusha/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,data):path.write_text(json.dumps(data,indent=2)+'\n')

METHODOLOGY='''This submission contains eight ranked de novo miniproteins across six generated parent lineages targeting the human TfR1 apical domain. The leading designs explore a His41-based acidic-release mechanism; their original generated parent is included as a matched comparator.

I used BoltzGen with Adaptyv-derived Modal wrappers against a 6WRW-bound receptor crop (native residues 190–383; hotspots 208/210/211). Refolding used receptor-only templates, excluding binder and tag coordinates. Known TfR1 binders, scrambled sequences and an unrelated VHH calibrated the workflow. Per-round selection rules assessed interface confidence, fold retention, receptor-fitted interface placement, hotspot contacts and full-dimer/glycan-context clashes. Failed selected models remain documented; favorable alternate samples did not replace them. Published binders were controls, not submitted sequences or starting binder templates.

Histidine substitutions and local ProteinMPNN redesigns were applied to our own de novo parents. Uniform restrained preparation and matched free/bound PROPKA calculations tested pKa direction. The A41H/A34K/A38R candidate retained the core screen in 9/10 bare and 10/10 tagged predictions; fresh bare pKa shifts were −1.12 and −0.96. These are computational hypotheses, not measured affinities. Free-state pKa heuristics, some control-preparation failures and uncertain tag placement prevent a fully qualified pH-switch claim. Unbound conformations, mouse cross-reactivity, expression and experimental pH release remain unvalidated. The submission omits assay tags; a documented platform-default suffix was modeled only as a sensitivity test.

Ranking balances the proposed pH mechanism, structural consistency and lineage diversity. The close-known-backbone short_015 backup is excluded. Parent sequence/structure database reviews and their limits are disclosed; later mutants do not inherit a new exhaustive novelty search. OpenAI Codex assisted implementation, documentation research and analysis. Scientific tools included BoltzGen, ProteinMPNN, PROPKA, OpenMM, ipSAE, Foldseek, TM-align and BLAST+. No Claude usage or wet-lab validation is claimed.

The attached ZIP contains detailed methods, candidate-level limitations, evidence provenance and 16 selected complexes whose binder sequences match the CSV. Repository: https://github.com/anushajain2/adaptyv-challenge4
Final upload package: https://github.com/anushajain2/adaptyv-challenge4/tree/main/reports/submission_upload_2026-11-02
'''

def main():
    original_manifest=json.loads((SOURCE/'manifest.json').read_text())
    for name,digest in original_manifest['files'].items():assert sha(SOURCE/name)==digest
    assert not OUT.exists(),'Do not overwrite an existing prepared upload package'
    shutil.copytree(SOURCE,OUT)
    subprocess.run([str(NODE),'scripts/submission_csv_20261102.mjs','revised',str(OUT),'single_chain'],cwd=ROOT,check=True)
    old=list(csv.DictReader((SOURCE/'submission.csv').open()));new=list(csv.DictReader((OUT/'submission.csv').open()))
    assert [(r['name'],r['sequence']) for r in old]==[(r['name'],r['sequence']) for r in new]
    assert len(new)==8 and all(r['molecule_class']=='single_chain' for r in new)
    assert (OUT/'submission.csv').stat().st_size<5_000_000
    evidence=json.loads((OUT/'evidence.json').read_text())
    evidence['upload_schema_note']='The live form screenshot lists single_chain; older scientific metadata/FAQ used protein. These describe the same single-chain miniprotein format. Sequences and ranks are unchanged.'
    for c in evidence['candidates']:c['molecule_class']='single_chain'
    save(OUT/'evidence.json',evidence)
    (OUT/'methodology.txt').write_text(METHODOLOGY)
    (OUT/'form_text.txt').write_text('DESCRIBE YOUR METHODOLOGY — paste the following text\n\n'+METHODOLOGY)
    (OUT/'README.md').write_text('''# TfR1 Track 3 upload package

**Prepared for the participant to upload; not submitted.** This version matches the fields shown in the November 2 submission-form screenshot. It contains the recommended eight sequences across six parent lineages, unchanged from the previous recommended package.

## Upload designs

Use [submission.csv](submission.csv), then click **Check designs**. Columns are `name`, `sequence`, and `molecule_class`; all class values are `single_chain`, as listed by the form. The public FAQ previously used `protein`. Only this class label changed: sequence identities and ranking are identical. No assay suffix is appended. [submission.fasta](submission.fasta) is an alternative sequence format; if used, set each design's class to single chain in the form.

The form automatically checks a novelty score of at least 3/4. This server-side score has **not** been run or verified by the local pipeline. Review any returned failures before submitting; do not assume our database searches guarantee acceptance. The screenshot also limits submissions to one every 24 hours.

## Describe your methodology

Paste the entire contents of [methodology.txt](methodology.txt) into the text box. Attach the accompanying ZIP under **Attach supporting files**. It contains [METHODS.md](METHODS.md), [evidence.json](evidence.json), the latest experiment reports, and 16 sequence-matched structures in `models/`. The form accepts `.md` or `.zip`, up to five files of 50 MB each; the package is below 1 MB. Uploading the ZIP alone supplies all supporting files.

## Link design methods

This field is optional. Leave it blank unless an appropriate existing Proteinbase design-method entry is available. The methodology already identifies the scientific tools and includes the public repository link. Do not invent a registered method entry.

## Confirm and submit

Review the eight accepted designs, their order, the methodology and the competition terms. Complete the terms checkbox and final submission yourself. The screenshot states that methodology and supporting files may become public if designs are selected. No terms have been accepted or submission made by this packaging workflow.

The published deadline is **November 2, 2026 at 17:29 IST**. The current list replaces short_015 with the A41H/A34K/A38R reserve as an explicitly unvalidated pH hypothesis. Scientific gate failures remain in the evidence; neither tag compatibility nor experimental pH release or mouse binding is established. Earlier packages remain preserved separately.
''')
    validation=json.loads((OUT/'validation.json').read_text())
    validation.update(upload_molecule_class='single_chain',sequences_and_ranking_unchanged_from_recommended=True,
        csv_below_5_MB=True,server_novelty_check='NOT_RUN',terms_accepted=False,submitted=False,
        schema_source='User-provided screenshot, Screenshot 2026-11-02 at 4.44.19 PM.png',checked_utc=datetime.now(timezone.utc).isoformat())
    save(OUT/'validation.json',validation)
    save(OUT/'manifest.json',{'status':'PREPARED_NOT_SUBMITTED','track':3,'created_utc':datetime.now(timezone.utc).isoformat(),
        'source_package':str(SOURCE.relative_to(ROOT)),'source_package_manifest_sha256':sha(SOURCE/'manifest.json'),
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'files':{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='manifest.json'}})
    archive=ROOT/'outputs'/f'{OUT.name}.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file():z.write(p,str(p.relative_to(OUT.parent)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for name,digest in json.loads((OUT/'manifest.json').read_text())['files'].items():assert hashlib.sha256(z.read(OUT.name+'/'+name)).hexdigest()==digest
    assert archive.stat().st_size<1_000_000
    save(ROOT/'reports'/f'{OUT.name}_archive.json',{'archive':str(archive.relative_to(ROOT)),'sha256':sha(archive),'bytes':archive.stat().st_size,'all_member_hashes_verified':True})
    print('Verified upload package:',archive,'bytes:',archive.stat().st_size,'methodology words:',len(METHODOLOGY.split()))

if __name__=='__main__':main()
