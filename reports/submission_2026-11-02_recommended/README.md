# TfR1 Track 3 submission prepared for review

**Prepared, not submitted.** Upload [submission.csv](submission.csv) to the [official form](https://proteinbase.com/competitions/anthropic-adaptyv-2026/rounds/4/submit) after participant review. The user confirmed Track 3. The package contains **8 bare sequences across 6 generated parent lineages**. None has experimentally established binding, pH release or mouse cross-reactivity.

The published deadline is **November 2, 2026 at 17:29 IST**, equivalent to November 1 at 23:59 AoE. [Official requirements](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1) were checked on November 2. Track 3 permits up to 20 designs; the required CSV columns are `name`, `sequence`, `molecule_class`. All entries here are `protein`, with canonical sequences of 10–250 residues. No assay suffix or known control binder is included.

This final editorial list replaces short_015 with the A41H/A34K/A38R reserve, explicitly as an unvalidated mechanism hypothesis. It does not reclassify the reserve as passing the frozen scientific gate. The original release remains unchanged.

## Ranked candidates

| Rank | Design | Residues | Evidence and key limitation |
| --- | --- | ---: | --- |
| 1 | apical_short_008_A41H_mpnn03 | 58 | Replicated pKa-direction challenger; latest core passes 9/10 bare and 10/10 tagged. Free-pKa band and full control/tag qualification remain unmet. Same lineage as A41H and short_008; no new exhaustive novelty search. |
| 2 | apical_short_008_A41H | 58 | Repeated bare pKa-direction hypothesis; latest bare core passes 9/10. Latest tagged preparation fails in one seed. Free-pKa band and full qualification remain unmet. |
| 3 | expansion_apical_long_092 | 89 | Repeated parent-interface retention; pH release, mouse binding and tag independence unestablished. AAAA run is a descriptive sequence flag. |
| 4 | expansion_apical_short_026 | 62 | Repeated bare retention; both earlier selected tagged core poses retained. Tag independence unresolved; NNT sequon is a motif flag, not proven glycosylation. |
| 5 | diversity_long_037 | 95 | Latest raw confirmation retains both selected poses (8/10 samples). One latest prepared pose fails movement. NTS sequon; no engineered His or established pH mechanism. |
| 6 | diversity_short_055 | 64 | Confirmed parent with 9/10 sample passes. Models assume C31–C42 disulfide formation; lower confidence and no demonstrated pH mechanism. |
| 7 | expansion_apical_short_018 | 61 | Weaker independent lineage: later confirmation has one selected raw clash failure. Explicit long-shot binding hypothesis. |
| 8 | apical_short_008 | 58 | Matched parent of A41H and reserve, not an independent scaffold. Latest bare core passes 10/10; no established pH-release mechanism. |

## Files for submission

Use `submission.csv` for sequence upload and [form_text.txt](form_text.txt) for the title, workflow and AI-disclosure fields. [METHODS.md](METHODS.md) describes the final scientific scope. The [evidence record](evidence.json) links original reports, preserves failures and supplies sequence reviews; `models/` contains both selected bare complexes for every candidate, with chain A receptor and chain B binder. The latest available rounds differ among candidates; no samples were selected because of favorable pKa results. [FASTA](submission.fasta) is a sequence-review companion, not a substitute for the required CSV.

The [latest reserve comparison](LATEST_RESERVE_RESULTS.md) and [local redesign comparison](LOCAL_REDESIGN_RESULTS.md) preserve post-shortlist results. Three related short_008 sequences, if included, still count as one lineage. Model consistency and PROPKA ratios are computational hypotheses, not laboratory measurements. Flexible-tag placement, unbound conformations, disulfide formation where relevant and mouse binding remain unresolved.

## Provenance and status

Scientific source snapshot: [commit c04fa3f](https://github.com/anushajain2/adaptyv-challenge4/tree/c04fa3fc05216bbe8c08fa7127a783778793f540). The manifest records sequence, model and supporting-file hashes. Older frozen releases remain unchanged. Packaging adds no Modal compute; the last conservative campaign accounting was $13.47 of the authorized $30, including credits.

Participant review and any required form declarations remain the participant's responsibility. This package does not accept terms or assert that an online submission has been made. No upstream issue or pull request was filed.
