# TfR1 challenge: computational binder and pH-direction campaign

This project prepares a human TfR1 target, generates and evaluates binder designs on Modal, and screens histidine variants for a predicted acidic-release direction. A computational pKa shift is a hypothesis, not a measurement of binding or proof of a pH switch.

**Campaign reopened at the user's request (25 October):** submission remains on
hold. A [local ProteinMPNN redesign comparison](reports/local_redesign/README.md)
keeps A41H or S90H and likely binding anchors fixed while changing nearby residues.
The eight-design release remains an immutable fallback. This round is limited to
$6 within the original $30 gross cap; no upstream bug report or PR is authorized.
The [completed comparison](reports/local_redesign_panel/RESULTS.md) evaluated six
redesigns in 130 predictions. Three pass both selected raw structural checks.
**A41H with A34K/A38R is a promising reserve:** both prepared poses pass and its
predicted pKa shift is about −1.04 in both. It does not meet the full preset rule
(free pKas 7.12–7.21, plus baseline/control preparation limitations), so no design
was promoted and the conditional tag job was skipped. The observed gross increase
is about **$0.32**, with conservative cumulative accounting **$13.04/$30**. All
cloud apps are stopped; submission remains on hold.

**Previous stopping-point experiment (25 October):** the [12-variant His panel](reports/final_his_panel/RESULTS.md)
completed 160 predictions. Four variants pass both selected raw structural checks,
but **none qualifies as a two-pose pH-direction hypothesis**. All 32 selected
models were prepared uniformly; a positive-control movement failure is disclosed.
The conditional tag job was skipped. Long_037 passes both fresh selected raw poses
again (7/10 samples, following 8/10 previously). The **eight-sequence, seven-lineage
shortlist is unchanged**, with updated evidence and all twelve variants excluded.
This was the previous stopping point, superseded by the local redesign comparison.
The [frozen package](reports/submission_release_2026-10-25/README.md) is retained;
no challenge submission has been made. The final His panel added
about **$0.48** gross; conservative cumulative accounting is about **$12.72/$30**
(the account meter after that round was $12.57). All jobs from that round stopped.

**Earlier diversity round:** the [bounded diversity screen](reports/diversity_generation/RESULTS.md)
generated 100 designs (99 unique), with five full initial passes and two parents
retained after [two-seed confirmation](reports/diversity_confirmation/RESULTS.md).
**diversity_long_037** passes 8/10 samples and **diversity_short_055** passes 9/10;
both selected poses of each also pass preparation guards. All selected controls
pass. The [new database review](reports/diversity_novelty/RESULTS.md) is complete,
with no close-sequence flags and only weak/partial followed-up structural matches.
Short_055 assumes a C31–C42 disulfide; long_037 carries an NTS sequon. Neither has
histidine or established pH release. The round used about **$1.48**; recorded gross
spend is about **$12.24/$30**, with all cloud apps stopped.

**Original-six review:** the [database similarity audit](reports/novelty_review/RESULTS.md)
finds no reported RCSB sequence hits for the six designs under the recorded settings,
with both sequence retrieval controls passing. All 13 Foldseek queries completed;
the known binder is recovered. **short_015 shares a known backbone arrangement with
AF1782/2OO2**: all 61 residues align in both poses at 1.450/1.467 Å TM-align RMSD.
It remains a known-backbone backup, not a novel-fold claim. Other hits include
partial matches to larger proteins. **Local Swiss-Prot BLASTp is complete:** all
seven queries ran against 575,748 canonical sequences (release 2026_03), the human
P02786 control was recovered exactly, and no candidate meets the close-sequence
review threshold. Weak local matches and masking limits are reported. The original
NCBI job was left pending at that review; the completed evidence is the separately archived local
fallback. Database limits prevent a global novelty verdict; sequences and ranks
remain unchanged. No Modal compute was added.

The [ranked submission draft](reports/submission_draft/README.md)
contains **eight bare sequences across seven parent lineages**, with sequence-liability
and diversity audits, methods, source hashes and 16 selected structures. The user
confirmed Track 2 or 3; eight is within the 20-design maximum. A41H ranks first as a
pH-direction hypothesis; its parent is explicitly a related comparator. E19H stays
in reserve. Mixed geometry, preparation exclusions and unresolved tag compatibility
remain visible. This is a review draft, not a submitted application. No new cloud
compute was used. [Draft CSV](reports/submission_draft/submission_draft.csv).

**Earlier construct experiment:** the [matched bare/tagged comparison](reports/tag_comparison/RESULTS.md)
completed **140 predictions** using the documented 64-residue platform-default
suffix. All 12 selected raw control checks pass using binder-core PAE. Both bare
variants pass the two-seed core screen; neither tagged variant does. Tagged A41H
has one 1.909 Å core close pair in a selected pose; tagged E19H misses the
interface-pose limit in both seeds. Candidate suffix placements have low confidence,
so their contacts neither validate tag independence nor prove interference.

**A41H remains the leading hypothesis; no qualified tagged two-pose pH switch is established.**
All 28 selected models were prepared uniformly with the suffix retained; 23 pass
core movement guards and all 56 free/bound PROPKA calculations execute successfully.
One new bare A41H model gives a −0.608 pKa shift and a 2.12 His-only ratio proxy;
the other is excluded by movement. E19H has one usable shift but fails the free-pKa
heuristic. No experimental binding or pH switch is established.

Next: review the frozen eight-sequence package and upload it through the official challenge page. The final His round is complete; no further generation or mutation runs are planned. Tag compatibility and mouse cross-reactivity remain unresolved. The
[earlier preparation review](reports/receptor_preparation_v2/RESULTS.md) and
[24-variant follow-up](reports/expansion_his_direction/RESULTS.md) retain their
original decisions.

The [previous 200-design expansion](reports/next_round/RESULTS.md) supplies four
computational parent scaffolds. Mouse cross-reactivity remains unknown after an
inconclusive matched-conformation comparison. At that earlier stage, recorded gross Modal
usage was approximately **US$10.75 before credits**, within the **US$30 total**;
billing may lag. All cloud jobs are stopped. No upstream issue, PR or challenge
submission has been filed. [Review selected structures and full confidence arrays](reports/tag_comparison/review_bundle/manifest.json).

**Hotspot correction:** feature-level preflight found that our specifications
combined `binding` and `not_binding` in one chain entry, but the pinned parser
reads only the latter in that case. All 75 earlier generated metadata files
confirm zero positive hotspot tokens. Their metrics remain valid; they do not
test the intended hotspot conditioning. Current specifications split the entries,
and the new comparison verifies the actual features before sampling.

**Calibration follow-up:** all 100 earlier control samples were audited, and
80 matched BoltzGen plus 20 Chai-1 predictions were completed. BoltzGen recovered
the held-out 2DS25.5 interface and fold in all 20 samples, strongest with 6WRW.
ipSAE at the predeclared cutoff missed 3DS18; sequence-only Chai reproduced the
2DS25.5 fold but recovered neither published interface. The replacement gate
remains unvalidated. A separately justified exploratory comparison now uses
6WRW, 10 unconstrained plus 10 sheet-constrained 88-residue miniproteins, with
two fresh refolding seeds and matched controls. Full results, PAE export fixes
and reproducibility records are in the results report.

## Current status

| Deliverable | Status |
|---|---|
| Final His experiment | 12 variants, 160 predictions, four raw structural passes, zero qualified two-pose pH hypotheses; conditional tag job skipped |
| Database similarity review | Original-six and added-two reviews complete with retrieval controls; no close-sequence flags; short_015 known-backbone match verified |
| Ranked submission draft | Eight sequences, seven parent lineages; bare CSV/FASTA, sequence audit, methods and 16 selected models; not submitted |
| Matched bare/tagged comparison | 140 predictions; raw controls pass; neither tagged variant clears both selected core screens |
| Tagged-round preparation/direction | 28 models, 56 successful pKa calculations; 23 movement passes; zero qualified two-pose direction passes |
| Uniform receptor preparation v2 | 28 models, 56 warning-free pKa calculations; 25 movement passes; complete prepared control-panel gate fails |
| Assay-tag review | Exact public default suffix found and tested as an explicit challenge assumption; compatibility unresolved |
| New apical expansion | 200 designs; 16 initial full passes; four new parents pass both fresh-seed checks |
| Latest His follow-up | 24 new variants; four raw structural confirmations; zero strict two-pose direction passes; four receptor-protonation scoring exclusions |
| Earlier His panel | 20 variants/comparator; A41H retained as a qualified direction hypothesis |
| Advisory mouse comparison | Inconclusive because AlphaFold-human positive control fails |
| Target structures, human/mouse conservation map, four hotspot sets | Generated |
| Modal/BoltzGen pilot implementation | Both protocols completed successfully |
| 50 real pilot designs | 25 miniproteins + 25 nanobodies downloaded; none pass the parent-confidence gate |
| 6WRW generation comparison | 10 unconstrained + 10 sheet-constrained designs; 230 predictions scored; 0 advance |
| Adaptive follow-up | 25 miniproteins completed; intended apical-tip hotspot labels were not applied; no qualified parents |
| Histidine-pair rebuild, relaxation, and PROPKA pipeline | Ran 8 pairs × 2 models across 4 actual pilot parents |
| pKa-lowering pairs from the pilot | 6 exploratory pairs lower both pKas in both models; 0 qualified candidates because parent confidence is low |

The approved cloud budget is **US$30 total including credits**. The original
pilot used 25 designs per protocol; the earlier 20-design comparison had a separate
**US$5 allowance within that total**, including its failed preflight attempt.
No outcome is guaranteed within that budget.

Local verification: **75 tests passed**. A published 6WRV reference pair was rebuilt,
relaxed and processed by PROPKA in free/bound forms without incomplete-group
warnings. It did not pass the provisional candidate thresholds and is not a
new pilot design. See `reports/day1_3_status.json`.

Actual pilot evidence: `outputs/analysis/pilot_designs.csv` contains all 75 unique
designs. `reports/pilot/all_histidine_pairs.csv` contains all eight scanned pairs.
`reports/pilot/exploratory_pka_lowering.fasta` retains the six diagnostic variants;
`reports/pilot/provisional_switches.fasta` is deliberately empty. The successful
direction models were audited against their actual structures and pKa tables.

## Repository contents

This repository includes the pipeline, pinned local dependencies, public target
inputs, prepared structures, archived run configurations, and compact reports.
The complete 75-design metrics table and sequences are also archived under
`reports/pilot/pilot_designs.csv` and `reports/pilot/pilot_sequences.fasta`.
Historical evidence paths in reports refer to the original local workspace.

Raw generated structures, logs, launch ledgers, installed dependencies, and vendor
checkouts are excluded from Git. Existing raw cloud outputs remain in the
`tfr1-pilot-results` Modal volume; use the download commands below to retrieve them.
A fresh clone does not contain the local launch ledger: reconcile prior usage and
run IDs before launching any additional paid work.

## One-time account step

Completed on this computer. These instructions are retained for a fresh setup.

1. Create a [Modal account](https://modal.com/).
2. In **Usage & Billing**, set the **Workspace budget / usage limit to US$30**. This is the gross usage cap, before credits. A net spend limit is different; see [Modal budgets](https://modal.com/docs/guide/budgets).
3. Authorize this project's client:

```bash
cd /Users/anusha/Documents/code/adaptyv-challenge4
.venv/bin/modal setup
```

Do not paste tokens into chat or save them in this repository.

## Prepared scientific inputs

- `data/raw/manifest.json`: public source URLs and SHA256 checksums.
- `data/prepared/human_apical.pdb` and `.cif`: 11QC chain A, native residues **190–383**, renumbered **1–194** for unambiguous BoltzGen indices.
- `data/prepared/human_11QC_apical_uniprot.pdb`: same crop retaining native numbers.
- `data/prepared/human_6WRV_apical_uniprot.pdb`: comparison crop.
- `data/prepared/human_6WRV_G142_ectodomain.pdb`: comparison receptor with Ser142 changed to Gly, retaining only glycine backbone atoms there.
- `data/prepared/human_11QC_ectodomain_dimer.pdb`: full observed receptor context for collision checking and pKa evaluation.
- `data/prepared/conservation_map.csv`: human/mouse sequence mapping, all divergent positions, SASA in the full dimer, SASA in the crop, glycan distances, crop boundaries and design indices.
- `data/prepared/epitope_divergence.csv`: divergent residues in the stated epitope, including the six highlighted in the brief.
- `data/prepared/hotspot_sets.json`: provisional exposed patches.

The primary 11QC structure already has Gly142. Residue 142 lies outside the chosen apical crop. Human residue 121 is not modeled in 11QC and is not invented.

**Important corrections and limits**

- 6W3H is a circularly permuted apical construct, not a simple contiguous native-sequence fragment. It is downloaded for reference but is not used as the primary target.
- The global BLOSUM62 alignment (gap open −10, extension −0.5) places a gap at human V210; the challenge's stated correspondence places it at L212. Both are recorded separately. Positions 210–212 are excluded from conserved hotspot selection; the sequence alignment does not resolve the structural correspondence.
- Conservation is not proof of cross-reactivity. No mouse structure or mouse complex has been modeled in this stage.
- SASA is computed with the full observed dimer and observed glycans. Crystal/cryo-EM glycans are incomplete; glycan clearance remains unproven.
- A 190–383 crop includes the apical-bearing region. Artificial boundary contacts are discouraged, and generated complexes are checked against the full receptor dimer.
- The nanobody run uses the upstream caplacizumab framework with redesigned CDRs. Framework provenance is retained; this is not a claim of an entirely de novo framework.

## Reproduce local preparation

The `.venv` is already installed here. For another machine:

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.lock.txt
git clone https://github.com/adaptyvbio/protein-design-skills vendor/protein-design-skills
git clone https://github.com/hgbrian/biomodals vendor/biomodals
git clone https://github.com/HannesStark/boltzgen vendor/boltzgen
git -C vendor/boltzgen checkout 247b9bbd8b68a60aba854c2968d6a0cddd21ad6d
.venv/bin/python scripts/fetch_inputs.py
.venv/bin/python scripts/prepare_target.py
.venv/bin/python -m pytest -q
```

Local packages are pinned in `requirements.lock.txt`. Cloud BoltzGen code is pinned to the version in the organizer-recommended wrapper. Transitive cloud dependencies are not fully locked; the successful cloud environment is recorded in each run's `requirements-cloud.txt`. Dependency resolution upgraded the initially requested PyTorch 2.5.1 to 2.14.1. The source pin and the saved environment record are both needed for provenance.

## Launch the pilot

Only after authentication and setting the actual dashboard budget:

```bash
.venv/bin/python scripts/launch_pilot.py --workspace-usage-cap-confirmed 30
```

The launcher uses `modal run --detach`, allowing the job to keep running if the client disconnects. It creates `outputs/cloud_launch.json` **before** image construction. An existing ledger blocks relaunch; do not delete it to retry without reconciling actual Modal usage first.

The wrapper is adapted from `vendor/biomodals/modal_boltzgen.py`, linked from Adaptyv's skills. It uses a CPU-only image build, pinned BoltzGen commit, one L40S GPU, serial protocols, no automatic retries, bounded resources and a four-hour function timeout. A US$20 planning reservation leaves headroom under the approved US$30 cap; only Modal's workspace budget enforces the monetary limit. Image build and startup time are also billable.

Both protocols run on the `conserved_326_347` hotspot set. Miniproteins have 80–110 residues, avoiding the 73–76 window noted in the brief. Nanobody CDR insertions follow the upstream framework recipe. `--num_designs 25` is the generation count per protocol; `--budget 10` asks filtering to retain up to ten, not ten guaranteed successful binders. Each spec passes `boltzgen check` remotely before generation.

Outputs and logs persist in the `tfr1-pilot-results` Modal volume. `pilot-002`
contains the original 50 designs; `pilot-003` contains the 25-design follow-up.
Archived inputs in `reports/run_inputs/` match the remote input hashes.
`pilot-001` is a failed pre-generation attempt, retained for provenance.

```bash
.venv/bin/python scripts/download_pilot.py pilot-002 outputs/pilot002_complete
.venv/bin/python scripts/download_pilot.py pilot-003 outputs/pilot003_complete
.venv/bin/python scripts/read_pilot.py outputs --out outputs/analysis
.venv/bin/python scripts/scan_pilot.py outputs/analysis/scan_manifest.json --out outputs/direction_qualified
```

Downloads are already present here; use a fresh destination if downloading again.
The reader deduplicates repeated copies by structure hash. The qualified scan
correctly selects zero parents from the current results. The performed exploratory
scan used `outputs/analysis_initial/scan_manifest.json`, `--exploratory`,
`--parents-per-protocol 2`, `--pairs 2`, and `--out outputs/direction_initial`.
Its results can be consolidated and plotted with:

```bash
.venv/bin/python scripts/summarize_direction.py --direction outputs/direction_initial --manifest outputs/analysis/scan_manifest.json --out reports/pilot
.venv/bin/python scripts/plot_results.py --designs outputs/analysis/pilot_designs.csv --pairs reports/pilot/all_histidine_pairs.csv --out reports/pilot
```

`read_pilot.py` reads real aggregate CSVs and refolded complexes, preserves raw metrics, extracts sequences, calculates mean Kyte–Doolittle hydropathy, normalizes chain identities from sequence, checks full-dimer clashes, and writes the scan manifest. It fails rather than fabricating missing results. Confidence metrics retain their upstream names/scales. No new Boltz-2 affinity prediction is part of this three-day pilot.

Binder mutation positions use one-based full-chain FASTA indices, not IMGT numbering.
After scanning, `scripts/summarize_direction.py` consolidates the pair-level evidence
and exports only variants meeting the stated provisional thresholds.

## Direction-check interpretation

The scanner:

1. Chooses interface-periphery positions using contacts and buried SASA; nanobody positions are restricted to mapped CDRs.
2. Enumerates pairs, excluding Gly/Pro/Cys and existing His; ranks nearby basic atoms as a heuristic and records oxygen contacts.
3. Rebuilds complete histidine side chains with PDBFixer. Two model seeds provide a limited sensitivity check, not exhaustive rotamer sampling.
4. Uses restrained Amber14/GBn2 OpenMM minimization with a 12 Å cutoff and a pH-7.4 hydrogen assignment. The target is fixed, binder backbone restrained, and terminal oxygen geometry initialized explicitly. This is not constant-pH molecular dynamics. Each rebuild/relaxation has a 180-second process timeout.
5. Fits the complex into the full receptor dimer and screens receptor–binder heavy-atom clashes. This graft is for geometric screening only: its crop junctions need not form valid peptide bonds.
6. Runs PROPKA on the intact trimmed mutant complex and the extracted free binder in the same conformation. It does not use the grafted structure for pKa. The receptor outside residues 190–383, its glycans, and native continuity at the artificial termini are omitted from this calculation.
7. Records `shift = pKa_bound − pKa_free`. Negative shifts favor acidic release for the modeled sites. PROPKA incomplete-group warnings prevent scoring; raw outputs and warnings are retained.

The two-His independent-site estimate is:

```text
Ka(pH)/Ka(reference) = product_i (1 + 10^(pKa_bound_i - pH))
                                 / (1 + 10^(pKa_free_i  - pH))
predicted Kd(6.0)/Kd(7.4) = Ka(7.4)/Ka(6.0)
```

The provisional filter asks for both shifts ≤ −0.5, at least one free pKa between 6 and 7, estimated ratio ≥2, no receptor-interface clash, and agreement across both model seeds. These are transparent pilot heuristics, not validated cutoffs. A half-unit shift is not proof given pKa prediction uncertainty. No passing pair is a valid outcome.

The parent must also pass interface ipTM ≥0.50, backbone refolding RMSD <2.5 Å, and zero receptor clashes. This RMSD gate is exploratory; the brief's standard is stricter at 2.0 Å. `scan_pilot.py --exploratory` permits a diagnostic scan of clash-free, low-confidence parents, but the consolidated exporter still excludes those variants from `provisional_switches.fasta`. The upstream selection of ten designs is not evidence that they pass these gates.

This estimate omits proton linkage at other sites, coupled titration, free-binder conformational changes, complete glycans, and experimental affinity. It cannot establish preserved affinity at pH 7.4. Parent-confidence metrics do not validate a mutated binder. Promising variants still need refolding and later experimental testing.

`tests/fixtures/reference_6WRV_not_pilot.pdb` is a published-complex test fixture. Any calculations under `outputs/reference_smoke` or `outputs/diagnostic` are infrastructure tests, **not pilot designs or competition candidates**.

## Primary sources

- [Competition target and assay](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1)
- [Adaptyv protein-design skills](https://github.com/adaptyvbio/protein-design-skills)
- [Recommended biomodals wrapper](https://github.com/hgbrian/biomodals/blob/main/modal_boltzgen.py)
- [BoltzGen source at the pinned commit](https://github.com/HannesStark/boltzgen/tree/247b9bbd8b68a60aba854c2968d6a0cddd21ad6d)
- [PROPKA](https://github.com/jensengroup/propka)
- [Human UniProt P02786](https://www.uniprot.org/uniprotkb/P02786/entry), [mouse Q62351](https://www.uniprot.org/uniprotkb/Q62351/entry)
- [11QC](https://www.rcsb.org/structure/11QC), [6WRV](https://www.rcsb.org/structure/6WRV), [6W3H](https://www.rcsb.org/structure/6W3H)
- [Modal pricing](https://modal.com/pricing)

## Reproduce the calibration checks

```bash
.venv/bin/python scripts/prepare_controls.py
.venv/bin/python scripts/audit_target_context.py
.venv/bin/python scripts/analyze_controls.py outputs/controls_complete/controls-002
.venv/bin/python scripts/analyze_controls.py outputs/controls_msa_complete/controls-004-msa --out reports/controls_msa
.venv/bin/python scripts/plot_controls.py
```

The analysis commands require the downloaded runs; use `scripts/download_pilot.py`
with the run IDs above. `scripts/launch_controls.py` reconciles live billing,
reserves US$8 of planning headroom under the existing US$30 gross cap, and records
each attempt before cloud image creation. It does not automatically retry.
`--reservation-usd` sets a smaller planning reservation. The matched control
round reserves US$5, limits prediction to 20 minutes, and exports full PAE and
pLDDT. Its optional Chai check reserves US$3 of the same round's remaining
headroom; it is not an additional US$5 allocation.
`scripts/fetch_control_msa.py` submits only the public receptor sequence to the
standard ColabFold MSA server; the saved 512-sequence alignment is in
`data/controls_msa/target.a3m`.

The baseline uses the pinned BoltzGen folding configuration with target-only
templates, empty MSA, three recycling steps, 200 sampling steps, and five diffusion
samples per seed (17 and 29). The MSA diagnostic keeps those settings and explicitly
injects the receptor MSA into the upstream `Input` object, increasing `max_seqs`
from 1 to 512. Preflight checks verify token counts, chain roles, exclusion of
binder templates, and nontrivial MSA depth for the MSA run.

### Reproduce the matched scoring audit

```bash
.venv/bin/python scripts/audit_control_samples.py \
  outputs/controls_complete/controls-002 \
  outputs/controls_msa_complete/controls-004-msa \
  --output reports/controls/all_samples_audit.json
.venv/bin/python scripts/prepare_matched_controls.py
.venv/bin/python scripts/download_pilot.py controls-005-matched outputs/controls_matched_complete
.venv/bin/python scripts/score_matched_controls.py \
  outputs/controls_matched_complete/controls-005-matched --output reports/controls_matched
.venv/bin/python scripts/plot_matched_controls.py
```

The published 6WRW source is archived in `data/raw/6WRW.cif`; its provenance and
the analysis plan fixed before prediction are in `data/controls_matched/manifest.json`.
The official ipSAE code is bundled unmodified with its license and commit/hash
provenance in `third_party/ipsae/`. Primary cutoffs are 10 Å PAE and 10 Å distance;
both directional scores are retained, with their minimum used for the primary
comparison. No affinity or experimental specificity is inferred from these scores.

The local sample reader checks exported CIF coordinates against the selected NPZ
sample before reconstructing any other sample. For matched scoring it restores
predicted sidechains omitted by the upstream resolved-atom export mask. Restoration
uses the pinned canonical atom order in `configs/boltzgen_atom_order.json` and
checks every retained atom. The coordinates come from the prediction, not an
imputation or sidechain rebuild. This is particularly important for the scrambled
control, whose input CIF intentionally retained only backbone atoms.

The bounded Chai cross-check uses the `chai_lab==0.6.1` API from the biomodals
wrapper with Torch 2.5.1, one seed (42), three recycles, 200 diffusion steps,
five samples per control, ESM embeddings, no MSA search, and no templates or
contact restraints. Its worker explicitly saves the returned PAE and pLDDT arrays.
This changes predictor and conditioning simultaneously; it is not experimental
validation or a pure predictor-only ablation.

```bash
.venv/bin/python scripts/download_pilot.py controls-006-chai outputs/controls_chai_complete
.venv/bin/python scripts/analyze_chai_controls.py \
  outputs/controls_chai_complete/controls-006-chai
```

6W3H author chains D/A form the tested receptor/Fc pair; PDB label-chain names
are different. The target retains the shared observed construct span (positions
6–149), with two unresolved internal residues retained in the input sequence.
The upstream writer omits their coordinates; analysis preserves that numbering.
These controls are not new designs. The scrambled sequence and unrelated VHH are
hypothesized negatives, not experimentally established nonbinders.


### Reproduce the 6WRW generation comparison analysis

The prospective specification is in `data/generation_comparison/manifest.json`.
The two arms differ only in a five-residue secondary-structure constraint, with
10 designs of length 88 per arm. Both use target-only templates and empty MSA
for refolding at fresh seeds 101 and 211, five samples per sequence and seed.
The published positive and two presumed negatives are controls, not candidates.

```bash
.venv/bin/python scripts/download_pilot.py generation-008-6wrw outputs/generation008_complete
.venv/bin/python scripts/analyze_generation_comparison.py \
  outputs/generation008_complete/generation-008-6wrw
.venv/bin/python scripts/plot_generation_comparison.py
```

Generation uses the pinned upstream design and inverse-folding stages; refolding
uses the same audited worker as the controls, retaining full PAE and all samples.
The launcher reserves spending before loading Modal and refuses unreviewed repeat
attempts. Both attempts share one $5 experiment ceiling within the $30 total.
Configured RNG seeds are archived; upstream NumPy generators and GPU operations
mean bit-for-bit reproducibility is not asserted.

The original calibration gate remains failed. This experiment has its own
prospective structural triage criteria; passing is not experimental validation.
Inverse-folded reference sidechains are placeholders, so design-contact recall
uses binder backbone contacts against target heavy atoms. Fully predicted atoms
are used for hotspot, PAE and clash measurements. This metric clarification was
recorded before examining new refolding results, with no threshold change.

The hotspot-format defect and its actual-output audit are documented in
`reports/generation_comparison/hotspot_conditioning_audit/audit.json`. Historical
run inputs are preserved, while current `configs/` specifications are corrected.
