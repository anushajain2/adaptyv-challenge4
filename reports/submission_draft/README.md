# TfR1 submission draft — six sequences, five parent lineages

**Draft for Tracks 2/3; not submitted.** A41H is the leading pH-release hypothesis. The other rows preserve stronger binding hypotheses, lineage diversity and a matched parent. None is experimentally validated, and none has established tag-compatible pH switching or mouse cross-reactivity.

The user confirmed Track 2 or 3. The [official requirements](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1) allow up to 20 designs for these tracks, require ranked CSV rows with `name`, `sequence`, `molecule_class`, and accept proteins of 10–250 residues. Six is an intentional shortlist, not a completed 20-design campaign. This package would not meet Track 1's 20-design minimum. Deadline: November 1, 2026, 23:59 AoE.

## Ranked draft

| Rank | Design | Length | Role and decisive limitation |
|---:|---|---:|---|
| 1 | apical_short_008_A41H | 58 | Best repeated bare pKa-direction hypothesis; latest tagged core fails one selected clash check; no qualified two-pose switch |
| 2 | expansion_apical_long_092 | 89 | Repeated bare structural retention; pH release and tag compatibility unestablished |
| 3 | expansion_apical_short_026 | 62 | Repeated bare retention and both tagged core poses retained; tag independence unresolved |
| 4 | expansion_apical_short_015 | 61 | Additional parent lineage; less replication and no tag evaluation |
| 5 | expansion_apical_short_018 | 61 | Additional lineage with mixed geometry: later selected clash failure retained |
| 6 | apical_short_008 | 58 | A41H's matched parent, not another independent lineage; latest bare selected clash failure |

This is an editorial ranking for experiments, not a probability-of-binding score. Rank 1 prioritizes a challenge-relevant pH hypothesis; ranks 2–3 have stronger parent-binding evidence. Rank 6 is useful for interpreting A41H, but its inclusion uses a slot for a near-identical comparator. A narrower four-sequence version would retain ranks 1–4; this six-sequence draft exposes the two weaker choices rather than disguising them as equally supported hits.

Use [submission_draft.csv](submission_draft.csv) for the machine-readable draft and [submission_draft.fasta](submission_draft.fasta) for sequences. They contain **bare designed sequences only**. The published platform-default suffix, known positive binders and presumed-negative controls are excluded. No sequence was edited to improve an audit score.

## What the audit found

The six cores have no cysteines or tryptophans. All have negative mean hydropathy and acidic sequence-estimated pIs, but these descriptors do not demonstrate soluble expression, folding or stability.

- **short_026 has an NNT sequon at residues 26–28** and NT motifs at 27 and 31. The sequon is a context-dependent glycosylation review flag, not evidence that this assay construct is glycosylated. It is retained unchanged.
- **long_092 has an AAAA run at 73–76**, with DG motifs at 9 and 45. These are compositional/chemical review flags, not validated degradation or aggregation predictions.
- **short_026 and short_015 contain methionine**, at 46 and 13 respectively. Oxidation susceptibility requires structural and experimental context.
- A41H, its parent and long_092 have calculated bare-core aromatic/disulfide A280 extinction coefficients of zero. This is a concentration-measurement caveat for bare proteins, not a binding failure; the actual tagged assay construct differs.

The complete positions, composition, pI, charge, hydropathy and provenance are in [audit.json](audit.json). [Biopython ProtParam](https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html) supplies the sequence descriptors; the pI calculation uses a 0–14 search interval to avoid clipping an acidic sequence at the default lower bound. Charges and pIs here are sequence estimates, separate from structure-dependent PROPKA and the pH-release mechanism. There is no calibrated developability score or claim of acceptable expression, immunogenicity, aggregation or shelf stability.

Across different parent lineages, global aligned sequence identity is **22.5–33.9%**. A41H and its parent are **98.3% identical** and count as one lineage. Sequence-aligned backbone comparisons and angle compositions are included as descriptive checks only: the alignment has partial coverage, is not a fold-search algorithm, and does not establish five different fold classes. All candidates come from the same BoltzGen apical campaign. Comparison against the three local controls is not an exhaustive novelty search against known proteins.

## Evidence and exclusions

[audit.json](audit.json) preserves per-round pass counts, sample-selection rules, both selected raw poses and their failures, the latest direction admissions/exclusions, and a reserve list. Pass counts are not pooled across changed selection rules. The [manifest](manifest.json) hashes source reports and all 12 [selected bare models](models/), including the failed short_018 and short_008 representatives. These models come from the latest available bare round for each design; rounds differ between designs. Tagged models and full confidence arrays remain in the [tag-comparison review bundle](../tag_comparison/review_bundle/manifest.json).

A41H's earlier four accepted independently predicted poses showed negative pKa shifts. The newest admitted bare pose adds a −0.608 shift with a 2.12 His-only ratio proxy; the other fresh pose is excluded by preparation movement. Re-preparing old poses does not create independent evidence. Neither tagged variant qualifies for direction interpretation under the frozen raw gate. The tagged controls calibrate core confidence, not flexible-suffix docking; low-confidence tag contacts do not prove interference. See [the latest report](../tag_comparison/RESULTS.md).

E19H stays in reserve: its two selected tagged poses fail interface retention, and its bare direction evidence remains incomplete. R3H, V60H and the H57/H60 pair do not supply consistent substantial release direction. D51H/T70H remains excluded because independent prediction lost the interface; its historical frozen-geometry pKa result does not qualify it. Known binders remain controls, not new designs.

## Recommended next decision

Use this package as the current submission baseline. Before increasing the list, perform a dedicated database-level sequence/structure novelty review on the five lineages, then consider a bounded generation round aimed at genuinely different apical scaffolds if more diversity is needed. Additional near-identical His variants or repeated low-confidence full-tag docking should not be used merely to fill twenty slots. Any new candidates should meet predeclared, control-checked structural rules and carry equally explicit limitations.

Do not claim mouse binding, affinity, pH switching or assay-tag compatibility. Those remain experimental questions. The exact 64-residue public platform suffix was modeled as an assumption about the challenge construct, not independently confirmed by its challenge page. The [methods draft](METHODS.md) states the computational scope and AI assistance.

This draft adds **$0 Modal compute**; latest recorded gross usage remains **$10.75348933 of $30**, subject to billing lag. No submission or upstream issue/PR was made. Reproduce locally with `.venv/bin/python scripts/build_submission_draft.py`; run `.venv/bin/python -m pytest -q` to verify the package and pipeline.
