# Final anchor-preserving His panel — 25 October 2026

**12 variants; 160 mutant/parent/control predictions across two seeds; 4 variants pass both selected raw models; 0 fully qualified two-pose direction hypotheses.** These are computational checks, not biological validation.

## Frozen experiment

The contact audit used the two independently predicted long_037 poses from seeds 2309 and 2411. All charged residues, glycine/proline/cysteine, high-contact anchor proxies and five residues at each binder terminus were preserved. Sites had to be peripheral in both poses (free SASA ≥25 Å², buried SASA 2–65 Å², target distance ≤5 Å). Ten neutral single substitutions and two pairs were fixed before prediction. The pairs include S37, which has a basic target neighbour in both source poses; this geometric rationale does not establish favorable electrostatics or a binding-energy contribution.

All twelve variants, the unchanged parent and three controls received five samples at each fresh seed 2503 and 2609. Entire binder templates were excluded and the target-only template/empty-MSA setup retained. The existing target-integrity-aware confidence selection, structural checks and direction thresholds were unchanged. Favorable unselected samples cannot rescue a failed selected pose.

## Results

| Variant | Raw passing samples / 10 | Both selected raw pass | Both-pose direction heuristic | Fully qualified |
|---|---:|---|---|---|
| diversity_long_037_T6H | 0 | False | False | False |
| diversity_long_037_S37H | 0 | False | False | False |
| diversity_long_037_A22H | 0 | False | False | False |
| diversity_long_037_N40H | 7 | True | False | False |
| diversity_long_037_V47H | 0 | False | False | False |
| diversity_long_037_S90H | 6 | True | False | False |
| diversity_long_037_T74H | 4 | True | False | False |
| diversity_long_037_N27H | 7 | False | False | False |
| diversity_long_037_A43H | 0 | False | False | False |
| diversity_long_037_S72H | 8 | True | False | False |
| diversity_long_037_T6H_S37H | 0 | False | False | False |
| diversity_long_037_A22H_S37H | 0 | False | False | False |

The unchanged parent passes 7/10 raw samples; both selected poses pass: **True**. All six selected raw controls pass: **True**. All six prepared controls pass: **False**.

All 32 confidence-selected models, including failures and controls, received one identical v2 preparation. Execution succeeded for 32; movement checks pass for 25. Warning-free PROPKA completed for 62/64 free/bound forms; warning failures remain excluded. Raw failures cannot be rescued by preparation. Full qualification requires both prepared parent poses and every prepared control to pass, in addition to the mutant-specific criteria.

### Per-pose direction evidence

| Variant | Seed | Status | Bound−free pKa shifts | His-only ratio proxy | Exclusion / failed checks |
|---|---:|---|---|---:|---|
| diversity_long_037_T6H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_T6H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_S37H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_S37H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_A22H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_A22H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_N40H | 2503 | ok | -0.120 | 1.181 | all_shifts_le_minus_0_5, his_only_ratio_ge_2 |
| diversity_long_037_N40H | 2609 | ok | -0.154 | 1.237 | all_shifts_le_minus_0_5, his_only_ratio_ge_2 |
| diversity_long_037_V47H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_V47H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_S90H | 2503 | ok | -0.475 | 1.713 | all_shifts_le_minus_0_5, his_only_ratio_ge_2 |
| diversity_long_037_S90H | 2609 | ok | -0.693 | 1.869 | his_only_ratio_ge_2 |
| diversity_long_037_T74H | 2503 | ok | -0.096 | 1.143 | all_shifts_le_minus_0_5, his_only_ratio_ge_2 |
| diversity_long_037_T74H | 2609 | excluded | — | — | AssertionError: Prepared movement/internal geometry failed |
| diversity_long_037_N27H | 2503 | excluded | — | — | AssertionError: Both raw poses and controls required |
| diversity_long_037_N27H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_A43H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_A43H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_S72H | 2503 | excluded | — | — | AssertionError: Prepared movement/internal geometry failed |
| diversity_long_037_S72H | 2609 | ok | -0.775 | 2.005 | none |
| diversity_long_037_T6H_S37H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_T6H_S37H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_A22H_S37H | 2503 | excluded | — | — | AssertionError: Raw selected screen failed |
| diversity_long_037_A22H_S37H | 2609 | excluded | — | — | AssertionError: Raw selected screen failed |

### Preparation-control limitations

- 2ds25_5_6wrw, seed 2609: AssertionError: Prepared movement/internal geometry failed. Backbone displacement: {"A": {"backbone": {"max_A": 1.0173554420471191, "rms_A": 0.15241256262284525}, "sidechain": {"max_A": 2.9885406494140625, "rms_A": 0.45779587754295453}}, "B": {"backbone": {"max_A": 0.4457177221775055, "rms_A": 0.13035822696983385}, "sidechain": {"max_A": 2.960676670074463, "rms_A": 0.5424626833825977}}, "pass": false}
- scrambled_6wrw, seed 2609: AssertionError: Prepared movement/internal geometry failed. Backbone displacement: {"A": {"backbone": {"max_A": 1.1000992059707642, "rms_A": 0.19700355874196288}, "sidechain": {"max_A": 2.804908037185669, "rms_A": 0.5872894560751656}}, "B": {"backbone": {"max_A": 0.9564506411552429, "rms_A": 0.3368230665104391}, "sidechain": {"max_A": 3.835548162460327, "rms_A": 1.0576857990740953}}, "pass": false}

## Decision

**No variant meets the prospectively frozen admission rule. The conditional tag job is not launched, and none of these twelve variants enters the submission CSV.** The current eight-sequence shortlist remains the baseline; the parent’s new confirmation evidence must be included in its ranking rationale. This is the planned stopping point for computational experimentation, not a reason to relax thresholds or start another generation batch.

## Evidence and limits

The [prospective plan](plan.json), [contact audit](../../data/final_his_panel/contact_audit.json), [all raw sample scores](results.json), [prepared models](preparation/results.json), [direction exclusions](../final_his_direction/results.json) and [review bundle](review_bundle/manifest.json) preserve failures as well as passes. The bundle includes all selected raw/prepared structures and complete five-sample confidence arrays for every case/seed.

PROPKA compares free and bound forms in a fixed bound-like conformation. The introduced-His-only ratio omits other titratable sites, conformational populations and full proton linkage. It is not a predicted experimental KD ratio. Preparation is finite restrained minimization, not proof of equilibrium or global stereochemical validity. No new mouse-binding or assay-tag claim follows from this panel.

The entire final round is limited to $3 gross including credits: $1.50 for this panel and $1.50 reserved for a conditional tag job. The local preparation and analysis add no Modal compute. No challenge submission or upstream issue/PR is performed.

Recorded gross increment: **$0.48**. Latest account meter: $12.57; conservative cumulative accounting, retaining the prior higher snapshot: **$12.72/$30**. The meter may lag. All cloud jobs are stopped; the conditional tag reservation is unused.
