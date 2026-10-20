# TfR1 days 1–3: results

The workflow is built and has run on real designs. Target preparation and the 50-design pilot are complete. The direction check found **six exploratory histidine-pair variants with both pKas lowered in both modeled conformations**. **No variant passes the provisional parent-quality gate**, and binding remains unvalidated. The calibration below now shows that confidence failure alone is not sufficient evidence that a sequence cannot bind.

## Matched scoring audit — 21 October 2026, India time

**The new positive control works well, but ipSAE does not provide a universal replacement gate.** We audited all 100 previous diffusion samples and completed 80 new matched predictions: four binders × two bound receptor conformations × two seeds × five samples. No novel designs were generated in this round.

2DS25.5 (the held-out positive from [6WRW](https://www.rcsb.org/structure/6WRW)) recovered its interface in all 20 samples, with whole-binder self-fit Cα RMSD 0.44–0.70 Å. Its ipTM also clears 0.50 in all 20 samples. Thus the earlier suggestion that the confidence threshold is simply wrong for this target is too broad: the predictor can reproduce a complete known binder, and confidence distinguishes that result from the poorly reproduced 3DS18 fold.

| Binder / receptor | Interface recovered | ipTM range | Minimum directional ipSAE | Contact-pair PAE (Å) | Whole-binder self-fit RMSD (Å) |
|---|---:|---:|---:|---:|---:|
| 3DS18 / 6WRV | 10/10 | 0.278–0.297 | 0.000–0.000 | 14.44–14.97 | 14.21–14.61 |
| 2DS25.5 / 6WRV | 10/10 | 0.506–0.529 | 0.059–0.082 | 7.57–8.10 | 0.51–0.69 |
| Scrambled, presumed negative / 6WRV | No reference pose | 0.122–0.141 | 0.000–0.000 | 22.23–25.14 | — |
| Unrelated VHH, presumed negative / 6WRV | No reference pose | 0.097–0.120 | 0.000–0.000 | 24.34–27.49 | — |
| 3DS18 / 6WRW | 9/10 | 0.334–0.352 | 0.000–0.000 | 13.61–14.01 | 14.28–14.70 |
| 2DS25.5 / 6WRW | 10/10 | 0.626–0.647 | 0.202–0.224 | 5.54–5.82 | 0.44–0.70 |
| Scrambled, presumed negative / 6WRW | No reference pose | 0.104–0.143 | 0.000–0.000 | 21.92–25.93 | — |
| Unrelated VHH, presumed negative / 6WRW | No reference pose | 0.101–0.127 | 0.000–0.000 | 23.67–26.90 | — |

The primary ipSAE settings were fixed before inference: PAE cutoff 10 Å, distance cutoff 10 Å, and the smaller of the two directional scores. The official code is bundled unmodified under `third_party/ipsae/`, pinned to commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3`. 3DS18 / 6WRV has minimum interchain PAE 10.51–10.86 Å, so no pairs survive that PAE cutoff and ipSAE is exactly zero. Changing the cutoff after seeing these results would be an exploratory sensitivity analysis, not validation of the predeclared gate.

Contact-pair PAE separates both published positives from these presumed negatives in this panel. This is useful evidence for ranking, but it is not experimentally calibrated specificity, and it cannot establish a correctly folded binder. It differs from upstream `interaction_pae`, which averages all interchain residue pairs rather than just contacting pairs. Neither score predicts affinity or pH switching.

**All-sample audit and RMSD frame check:** baseline 3DS18 / 6WRV recovered the interface in 8/10 earlier samples. The two MSA conditions each contained one recovered sample out of ten, both missed by confidence selection. Both reported pose RMSDs use the same receptor-derived rigid transformation, with no separate interface superposition. The large whole-binder discrepancy is concentrated in residues 1–20; the interface and much of the remaining chain are close to the reference. No biological explanation for that discrepancy has been established.

**Export correction:** previous archives contained PAE summaries, not the full matrices needed for ipSAE. The new worker preserves full PAE and pLDDT without changing sampling. The scrambled-input CIF contained only backbone atoms; the upstream output writer retained that resolved-atom mask despite predicting sidechain coordinates. Local scoring restores the actual predicted sidechains from NPZ using the pinned canonical atom order, checking every retained atom before adding any missing one. No coordinates are imputed. This avoids both an ipSAE auxiliary-metric failure and biased contact counting from missing negative-control sidechains.

**Decision:** the predeclared ipSAE/whole-fold generation gate did not pass. The 20-design comparison remains unlaunched. The completed, sequence-only Chai-1 check below also failed to recover the experimental interfaces. It changes both predictor and template conditioning, so it is not a clean predictor-only comparison. Published binders remain controls, not competition candidates.

**Chai-1 cross-check:** all four controls completed (one seed, five samples each), using `chai_lab==0.6.1`, ESM embeddings, no MSA, no structural templates and no contact restraints. Full PAE and pLDDT arrays were saved explicitly.

| Chai control | Native contact recall | Binder interface RMSD after receptor fit (Å) | Whole-binder self-fit RMSD (Å) | Minimum directional ipSAE |
|---|---:|---:|---:|---:|
| 3DS18 | 0% in all 5 | 33.89–36.04 | 15.49–16.53 | 0 |
| 2DS25.5 | 0% in all 5 | 35.20–36.50 | 0.52–0.65 | 0.055–0.081 |
| Scrambled / unrelated VHH | No experimental reference pose | — | — | 0 for both |

Chai reproduced the 2DS25.5 fold but placed it at a different interface. A rule accepting ipSAE merely because it exceeds these negatives would therefore accept an incorrect experimental pose. This is a concrete counterexample to treating score separation as sufficient calibration. Chai's failure here is specific to this small, sequence-only configuration; it does not establish that Chai generally fails on TfR1. Receptor fits to the experimental crops were approximately 2.8–3.7 Å, another reason to avoid treating grafted clash counts as definitive when the context fit fails. See [Chai metrics](controls_chai/results.json) and [executed settings](controls_chai/run_record/settings.json).

**End-of-round decision:** calibration is complete for this timebox. No further diagnostic jobs or generation jobs are running. Keep the original 75 designs exploratory, do not resume histidine engineering, and do not claim a validated replacement gate. The strongest complete-complex reference in this panel is **2DS25.5 with the 6WRW receptor**, which would be the starting reference for a separately justified, narrow generation experiment. The predeclared broader gate remains failed; no cutoff was changed to make it pass.

![Matched calibration and RMSD audit](controls_matched/matched_calibration.png)

Evidence: [all 100 original samples](controls/all_samples_audit.json), [all 80 matched samples and scores](controls_matched/results.json), [frozen inputs and analysis plan](../data/controls_matched/manifest.json), [executed worker and run records](controls_matched/run_record/). Raw per-sample structures, PAE arrays, and official ipSAE output are retained locally under `outputs/matched_scoring/controls-005-matched/`.

## Earlier calibration checkpoint — 21 October 2026, India time

**Hold further generation and histidine engineering until the scoring setup is better calibrated.** We completed 100 diffusion samples: six baseline conditions × two seeds × five samples, then four target-MSA conditions × two seeds × five samples. Each seed's confidence-selected model was compared with the experimental interface. No new binder designs were generated in this checkpoint.

| Control / receptor | Baseline selected ipTM | Baseline native-contact recall | With target MSA: selected ipTM | With target MSA: native-contact recall |
|---|---:|---:|---:|---:|
| 3DS18 / original 11QC crop | 0.174–0.184 | 0% | 0.290–0.297 | 0% |
| 3DS18 / bound 6WRV crop | 0.297–0.302 | 86–92% | 0.285–0.314 | 0% |
| ATV Fc / 11QC crop | 0.135–0.147 | 0% | Not run | Not run |
| ATV Fc / permuted 6W3H crop | 0.138–0.143 | 0% | Not run | Not run |
| Scrambled 3DS18 / 11QC | 0.145 | No experimental pose | 0.193–0.208 | No experimental pose |
| Unrelated VHH / 11QC | 0.115–0.132 | No experimental pose | 0.121 | No experimental pose |

For 3DS18 with the bound receptor and baseline settings, the binder's interface Cα RMSD after fitting the receptor is **0.77–0.95 Å**. Whole-binder pose RMSD remains high (18.5–18.8 Å), reflecting disagreement outside the recovered interface. Thus an ipTM cutoff of 0.6, or whole-binder RMSD alone, would miss this accurately recovered interface. This does not justify lowering thresholds for novel designs without further validation.

The 512-sequence target MSA was actually present in model features in every MSA condition. Binder templates were explicitly excluded. Adding the MSA failed to improve interface recovery in this small matched comparison; two selected models also distorted the receptor substantially. This does not establish that MSAs are generally harmful, but it refutes the assumption that merely enabling one fixes this setup.

The geometric controls were distinct from blind prediction:

- Both native published complexes are free of <2 Å receptor–binder clashes after mapped placement into the full 11QC context, including observed glycans.
- Replacing the bound receptor coordinates with 11QC introduces 15 close atom pairs for 3DS18, despite a good global alignment. Receptor conformation matters.
- 6W3H uses author receptor chain **D** with Fc chain **A**; label-chain IDs differ. Its circularly permuted target is mapped explicitly. We use construct positions 6–149, preserving two unresolved residues in the sequence, and the observed Fc fragment. This is a secondary construct-sensitive control, not an intact therapeutic antibody.
- The spatial audit adds 10 Å crop-cut and glycan-anchor buffers to the full-dimer SASA map. Original hotspots **327/331 pass** these checks. Positions 300/311 in the alternative patch fail the crop-exposure check. The data do not support claiming the conserved patch is simply inaccessible.

**Decision:** retain the 6WRV bound receptor conformation for further calibration; test an additional independent positive control or predictor and evaluate local interface geometry alongside confidence. The planned beta-sheet-biased versus unrestricted generation comparison remains gated. Do not restart PROPKA screening on new parents yet. Sequence conservation, glycan clearance and mouse binding remain separate requirements.

The controls are published binders ([6WRV](https://www.rcsb.org/structure/6WRV), [6W3H](https://www.rcsb.org/structure/6W3H)), and their structures may occur in predictor training data. The scrambled sequence and unrelated VHH are hypothesized negatives, not experimentally verified nonbinders. Two seeds provide limited sampling. No affinity, pH switch or cross-species performance has been established.

![Calibration confidence and contact recovery](controls/calibration_overview.png)

Evidence: [baseline metrics and poses](controls/results.json), [MSA comparison](controls_msa/results.json), [construct preparation](controls/preparation.json), [spatial audit](controls/spatial_target_audit.json), [spatial residue map](controls/spatial_target_map.csv). Archived execution settings, input hashes, environment versions, preflight checks and confidence arrays are under `controls/run_records/`.

## What was completed

- **Target preparation:** 11QC chain A residues 190–383, a 194-residue crop; native-numbered and sequential-numbered structures; human–mouse alignment and conservation map; four hotspot sets. The older 6WRV Ser142 was corrected to Gly in the comparison ectodomain. Primary 11QC already has Gly142, outside the crop.
- **Initial pilot:** 25 miniproteins and 25 nanobodies at the conserved 326–347 patch, using the Adaptyv-recommended Modal wrapper and pinned BoltzGen source. Both complete pipelines returned successfully.
- **Follow-up:** 25 more miniproteins at the exposed basic apical tip, with additional exclusions for surfaces exposed by cropping. This is a human-first, species-variable patch. This pipeline also completed successfully.
- **Direction check:** eight pairs across four clash-free original-pilot parents, two models per pair. Twelve models scored successfully; four models from two pairs failed PROPKA structure/protonation checks and were left unscored.

## Parent quality

| Batch | Designs | Best interface ipTM | Clash-free | Receptor clashes | Target fit failed |
|---|---:|---:|---:|---:|---:|
| Initial miniproteins | 25 | 0.215 | 2 | 17 | 6 |
| Initial nanobodies | 25 | 0.301 | 2 | 18 | 5 |
| Apical-tip follow-up | 25 | 0.258 | 1 | 20 | 4 |

All **75 designs are readable**, with no import errors. “Clash-free” alone does not establish binding. The provisional parent gate requires interface ipTM ≥0.50, backbone refold RMSD <2.5 Å, and zero full-receptor clashes. The brief's 2.0 Å RMSD standard is stricter; the confidence failure already excludes every design under either RMSD threshold. Upstream rank/retention is not a passing binding prediction.

The full receptor check revealed substantial contacts with omitted receptor regions. The small follow-up did not solve the confidence problem. These results support revisiting target context and design constraints before a large generation campaign; they do not establish that TfR1 is intrinsically undesignable.

## Exploratory pKa evidence

Six pairs lowered both pKas consistently. Three also met the stronger direction criteria: both drops ≥0.5, at least one free pKa between 6 and 7, and an independent-site estimated Kd(6.0)/Kd(7.4) ≥2 in both models.

**Every row below is exploratory and fails the parent-quality gate.** Mutation indices are one-based positions in the full binder sequence, not IMGT numbering. All three variants derive from the original conserved-patch pilot.

| Parent | Mutations | Least negative shifts across seeds | His-only estimated Kd ratio | Parent ipTM |
|---|---|---:|---:|---:|
| Nanobody 03 | N60H,Y62H | -1.17, -1.76 | 4.89–4.98 | 0.120 |
| Nanobody 03 | Y62H,T106H | -1.77, -0.80 | 4.56–4.64 | 0.120 |
| Miniprotein 01 | R44H,Y49H | -0.62, -1.36 | 3.04–3.04 | 0.146 |

The ratios are conditional model estimates, not predicted experimental affinities. Mutant binding strength and folding have not been validated. The six successful pairs' engineered residues are at least 24 Å from the artificial target terminal residues.

PROPKA uses the **intact 194-residue trimmed target** plus the mutant binder. The free binder retains exactly the same coordinates. A graft into the full receptor is used only for clash checks: graft junctions can have invalid peptide geometry, so that graft is not used for pKa. Failed earlier attempts are retained under `outputs/direction_exploratory_initial`; the corrected evidence is under `outputs/direction_initial`. Omitted receptor regions, glycans, artificial termini, coupled titration and conformational changes remain limitations.

## Cost and verification

The latest Modal snapshot reports approximately **US$3.19 gross usage**, within the authorized **US$30 including credits**. Reported billed amount is **US$0.00**, with credits covering compute. Billing can lag and prior amounts can be reconciled; [the latest raw record](controls_chai/billing_after.json) preserves the exact values, and [the earlier pilot snapshot](billing_final.json) is retained. All four successful control runs completed (baseline, MSA, matched BoltzGen, and Chai). The two older failed startup/parser attempts remain in the launch history; neither produced predictions. This round used a US$5 planning ceiling within the original US$30 gross cap. Both new apps are verified stopped with zero tasks.

**25 local tests passed.** Both initial cloud specifications and the follow-up specification validated. Control preflight checks verified residue masks and chain roles; the MSA comparison additionally verified binder-template exclusion and MSA depth 512. Downloaded input hashes match the archived inputs. All 12 scored direction models passed an actual-output audit for complete histidines, matching free/bound binder coordinates, valid PROPKA groups and arithmetic consistency.

## Files to use

- [Conservation map](../data/prepared/conservation_map.csv), [hotspot sets](../data/prepared/hotspot_sets.json), [trimmed target PDB](../data/prepared/human_apical.pdb).
- [All 75 designs and metrics](pilot/pilot_designs.csv), [parent sequences](pilot/pilot_sequences.fasta).
- [All eight histidine pairs and evidence paths](pilot/all_histidine_pairs.csv), [six exploratory variant sequences](pilot/exploratory_pka_lowering.fasta).
- [Direction summary](pilot/direction_summary.json), [output audit](initial_pilot/direction_audit.json), [machine-readable project status](day1_3_status.json).
- [Overview plot](pilot/pilot_overview.png), [PDF plot](pilot/pilot_overview.pdf), [reproduction instructions](../README.md).

`pilot/provisional_switches.fasta` is deliberately empty: zero variants pass all gates. The exploratory FASTA is not a submission recommendation. No account action is needed from you for this completed pilot.
