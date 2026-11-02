"""Package the selected Track 3 CSV and current evidence without changing old releases."""
import csv
import argparse
import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from novelty_search import ROOT, sha, save
from prepare_controls import read_structure, sequence
from build_submission_draft import sequence_audit
from freeze_submission_release import portable_markdown

OUT=ROOT/'reports/submission_2026-11-02'
RESERVE='apical_short_008_A41H_mpnn03'
NOTES={
 RESERVE:'Replicated pKa-direction challenger; latest core passes 9/10 bare and 10/10 tagged. Free-pKa band and full control/tag qualification remain unmet. Same lineage as A41H and short_008; no new exhaustive novelty search.',
 'apical_short_008_A41H':'Repeated bare pKa-direction hypothesis; latest bare core passes 9/10. Latest tagged preparation fails in one seed. Free-pKa band and full qualification remain unmet.',
 'expansion_apical_long_092':'Repeated parent-interface retention; pH release, mouse binding and tag independence unestablished. AAAA run is a descriptive sequence flag.',
 'expansion_apical_short_026':'Repeated bare retention; both earlier selected tagged core poses retained. Tag independence unresolved; NNT sequon is a motif flag, not proven glycosylation.',
 'diversity_long_037':'Latest raw confirmation retains both selected poses (8/10 samples). One latest prepared pose fails movement. NTS sequon; no engineered His or established pH mechanism.',
 'diversity_short_055':'Confirmed parent with 9/10 sample passes. Models assume C31–C42 disulfide formation; lower confidence and no demonstrated pH mechanism.',
 'expansion_apical_short_015':'Known-backbone backup: close AF1782/2OO2 match across all 61 residues. Structural novelty/eligibility concern; no novel-fold claim.',
 'expansion_apical_short_018':'Weaker independent lineage: later confirmation has one selected raw clash failure. Explicit long-shot binding hypothesis.',
 'apical_short_008':'Matched parent of A41H and reserve, not an independent scaffold. Latest bare core passes 10/10; no established pH-release mechanism.'}


def main():
    global OUT
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=OUT)
    OUT=parser.parse_args().output.resolve()
    assert OUT.is_relative_to(ROOT/'reports')
    assert not (OUT/'manifest.json').exists(),'Create a new release rather than overwrite a completed package'
    selected=json.loads((ROOT/'outputs/submission-authoring/selection.json').read_text())
    rows=list(csv.DictReader((OUT/'submission.csv').open()))
    assert [(r['name'],r['sequence'],r['molecule_class']) for r in rows]==[(r['name'],r['sequence'],r['molecule_class']) for r in selected['rows']]
    assert len(rows)==8 and len({r['sequence'] for r in rows})==8
    old=ROOT/'reports/submission_release_2026-10-25'
    frozen=json.loads((old/'manifest.json').read_text())
    for name,digest in frozen['files'].items():assert sha(old/name)==digest
    draft=json.loads((ROOT/'reports/submission_draft/manifest.json').read_text())
    historical=json.loads((ROOT/'reports/submission_draft/audit.json').read_text())
    history={c['name']:c for c in historical['candidates']}
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    base=f'https://github.com/anushajain2/adaptyv-challenge4/blob/{commit}'
    reserve=json.loads((ROOT/'reports/reserve_tag/results.json').read_text())
    reserve_direction=json.loads((ROOT/'reports/reserve_tag/direction.json').read_text())
    local=json.loads((ROOT/'reports/local_redesign_panel/results.json').read_text())
    local_direction=json.loads((ROOT/'reports/local_redesign_panel/direction.json').read_text())
    modeldir=OUT/'models';modeldir.mkdir()
    candidates=[];model_records=[]
    for rank,row in enumerate(rows,1):
        name=row['name'];seq=row['sequence'];assert row['molecule_class']=='protein' and 10<=len(seq)<=250
        c={**row,'rank':rank,'family':'apical_short_008' if name.startswith('apical_short_008') else name,
           'interpretation':NOTES[name],'sequence_audit':sequence_audit(seq),'experimental_validation':False}
        if name in history:c['historical_evidence']=history[name]
        if name.startswith('apical_short_008'):
            raw=[r for r in reserve['rows'] if r['base_id']==name and r['condition']=='bare' and r['confidence_selected']]
            c['latest_raw_report']='reports/reserve_tag/results.json'
            c['latest_selected_raw']=raw
            c['latest_preparation_and_direction']=[r for r in reserve_direction['rows'] if r['base_id']==name]
            sources=[(ROOT/r['raw_scoring_output']/'model.pdb',r['seed'],r['sample'],r['screen']['pass'],r['screen']['failed']) for r in raw]
        elif name=='diversity_long_037':
            raw=[r for r in local['rows'] if r['id']==name and r['confidence_selected']]
            c['latest_raw_report']='reports/local_redesign_panel/results.json';c['latest_selected_raw']=raw
            c['latest_preparation_and_direction']=[r for r in local_direction['rows'] if r['id']==name]
            sources=[(ROOT/r['raw_scoring_output']/'model.pdb',r['seed'],r['sample'],r['screen']['pass'],r['screen']['failed']) for r in raw]
        else:
            sources=[]
            for m in draft['models']:
                if m['name']==name:
                    src=ROOT/m['path'];assert sha(src)==m['sha256']
                    sources.append((src,m['seed'],m['sample'],m['selected_screen_pass'],m['selected_screen_failures']))
        assert len(sources)==2,name
        for src,seed,sample,passed,failed in sources:
            assert sequence(read_structure(src)[0]['B'])==seq,name
            dest=modeldir/f'{name}_seed{seed}_sample{sample}.pdb';shutil.copyfile(src,dest)
            model_records.append({'name':name,'seed':seed,'sample':sample,'source':str(src.relative_to(ROOT)),
                'path':str(dest.relative_to(OUT)),'sha256':sha(dest),'raw_screen_pass':passed,'failed_checks':failed})
        candidates.append(c)
    families=len({c['family'] for c in candidates})
    assert families==(6 if selected['choice']=='revised' else 7)
    assert all('WSHPQFEK' not in r['sequence'] for r in rows)
    source_paths=['reports/submission_draft/audit.json','reports/submission_draft/manifest.json',
        'reports/reserve_tag/results.json','reports/reserve_tag/direction.json','reports/reserve_tag/evidence_audit.json',
        'reports/local_redesign_panel/results.json','reports/local_redesign_panel/direction.json',
        'reports/novelty_review/summary.json','reports/diversity_novelty/summary.json',
        'reports/submission_release_2026-10-25/manifest.json','data/reserve_tag/manifest.json']
    save(OUT/'evidence.json',{'status':'PREPARED_NOT_SUBMITTED','track':3,'selection':selected['choice'],
        'candidate_count':8,'parent_lineage_count':families,'candidates':candidates,'models':model_records,
        'historical_evidence_note':'Historical nested ranks/rationales are preserved as dated provenance; current rank and interpretation are top-level fields.',
        'source_reports':{p:{'sha256':sha(ROOT/p),'url':f'{base}/{p}'} for p in source_paths}})
    for src,dest in [('reports/reserve_tag/RESULTS.md','LATEST_RESERVE_RESULTS.md'),('reports/local_redesign_panel/RESULTS.md','LOCAL_REDESIGN_RESULTS.md')]:
        (OUT/dest).write_text(portable_markdown(ROOT/src,commit))
    table='\n'.join(f"| {c['rank']} | {c['name']} | {len(c['sequence'])} | {c['interpretation']} |" for c in candidates)
    choice_note=('This final editorial list replaces short_015 with the A41H/A34K/A38R reserve, explicitly as an unvalidated mechanism hypothesis. It does not reclassify the reserve as passing the frozen scientific gate. The original release remains unchanged.' if selected['choice']=='revised' else 'Sequences and rank order are unchanged from the original frozen shortlist. short_015 retains a close known-backbone match, which creates an explicit structural-novelty/eligibility concern. The reserve is supporting evidence only, not a submitted sequence.')
    (OUT/'README.md').write_text(f'''# TfR1 Track 3 submission prepared for review

**Prepared, not submitted.** Upload [submission.csv](submission.csv) to the [official form](https://proteinbase.com/competitions/anthropic-adaptyv-2026/rounds/4/submit) after participant review. The user confirmed Track 3. The package contains **8 bare sequences across {families} generated parent lineages**. None has experimentally established binding, pH release or mouse cross-reactivity.

The published deadline is **November 2, 2026 at 17:29 IST**, equivalent to November 1 at 23:59 AoE. [Official requirements](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1) were checked on November 2. Track 3 permits up to 20 designs; the required CSV columns are `name`, `sequence`, `molecule_class`. All entries here are `protein`, with canonical sequences of 10–250 residues. No assay suffix or known control binder is included.

{choice_note}

## Ranked candidates

| Rank | Design | Residues | Evidence and key limitation |
| --- | --- | ---: | --- |
{table}

## Files for submission

Use `submission.csv` for sequence upload and [form_text.txt](form_text.txt) for the title, workflow and AI-disclosure fields. [METHODS.md](METHODS.md) describes the final scientific scope. The [evidence record](evidence.json) links original reports, preserves failures and supplies sequence reviews; `models/` contains both selected bare complexes for every candidate, with chain A receptor and chain B binder. The latest available rounds differ among candidates; no samples were selected because of favorable pKa results. [FASTA](submission.fasta) is a sequence-review companion, not a substitute for the required CSV.

The [latest reserve comparison](LATEST_RESERVE_RESULTS.md) and [local redesign comparison](LOCAL_REDESIGN_RESULTS.md) preserve post-shortlist results. Three related short_008 sequences, if included, still count as one lineage. Model consistency and PROPKA ratios are computational hypotheses, not laboratory measurements. Flexible-tag placement, unbound conformations, disulfide formation where relevant and mouse binding remain unresolved.

## Provenance and status

Scientific source snapshot: [commit {commit[:7]}](https://github.com/anushajain2/adaptyv-challenge4/tree/{commit}). The manifest records sequence, model and supporting-file hashes. Older frozen releases remain unchanged. Packaging adds no Modal compute; the last conservative campaign accounting was $13.47 of the authorized $30, including credits.

Participant review and any required form declarations remain the participant's responsibility. This package does not accept terms or assert that an online submission has been made. No upstream issue or pull request was filed.
''')
    save(OUT/'validation.json',{'status':'PASS','track':3,'candidate_count':8,'parent_lineages':families,
        'required_columns':True,'canonical_sequences_and_length':True,'unique_names_and_sequences':True,
        'bare_sequences_only':True,'known_controls_excluded':True,'all_16_model_chain_B_sequences_match_CSV':True,
        'frozen_release_unchanged':True,'new_modal_compute_usd':0,'checked_utc':datetime.now(timezone.utc).isoformat()})
    save(OUT/'manifest.json',{'status':'PREPARED_NOT_SUBMITTED','track':3,'created_utc':datetime.now(timezone.utc).isoformat(),
        'source_commit':commit,'files':{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file()}})
    archive=ROOT/'outputs'/f'{OUT.name}.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file():z.write(p,str(p.relative_to(OUT.parent)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for rel,digest in json.loads((OUT/'manifest.json').read_text())['files'].items():assert hashlib.sha256(z.read(OUT.name+'/'+rel)).hexdigest()==digest
    save(ROOT/'reports'/f'{OUT.name}_archive.json',{'archive':str(archive.relative_to(ROOT)),'sha256':sha(archive),'bytes':archive.stat().st_size,'all_member_hashes_verified':True})
    print('Prepared and verified:',archive)


if __name__=='__main__':main()
