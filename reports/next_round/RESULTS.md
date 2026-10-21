# Apical expansion and histidine follow-up: four new parents and one qualified pH-direction hypothesis

The round generated 200 additional apical miniproteins and completed 599 structural predictions across the initial screen, histidine panel, fresh-seed confirmations and advisory species comparison. Four new sequences pass the unchanged human-target screen in both selected fresh-seed models. The existing short-parent A41H variant retains a consistent computed pKa-lowering direction, but does not pass the strict direction heuristic in both poses. No experimental binding or pH switching is established.

## Retained evidence

| Lead | Confirmation samples passing / 10 | Selected interface ipTM | Current interpretation |
|---|---:|---|---|
| expansion_apical_long_092 | 9 | 0.670 / 0.676 | Highest sample consistency among the new parents; pH direction untested |
| expansion_apical_short_018 | 8 | 0.629 / 0.624 | New computational parent; pH direction untested |
| expansion_apical_short_026 | 8 | 0.606 / 0.614 | New computational parent; pH direction untested |
| expansion_apical_short_015 | 6 | 0.520 / 0.534 | New computational parent with lower sample consistency; pH direction untested |
| apical_short_008_A41H | 9 | 0.628 / 0.627 | Consistent computed release direction; strict free-pKa window failed in one pose |

Both selected models of each listed lead pass the raw structural screen. The four new parents also have zero counted <2 Å clashes after uniform preparation, with movement below the frozen bounds. Sequence identities among the five retained leads are 22.5–33.9% under the recorded global-alignment metric. Descriptive backbone angles suggest different helix/sheet compositions, but do not establish independent fold classes or method-family diversity. All originate from this BoltzGen apical campaign.

The [review bundle](review_bundle/manifest.json) contains all ten selected raw and prepared complexes, raw checks, sequences and hashes. [Lead diversity](lead_diversity.json), [FASTA](retained_human_leads.fasta), and the [working methods narrative](METHODS.md) make the evidence reviewable. These are a provisional research shortlist, not a finalized submission.

## What the experiments resolved

**Expansion:** The 200 new sequences include 13/100 short-arm and 3/100 long-arm initial full passes. All intended hotspot labels were verified in generated features. No exact sequence duplicates occur against the 295 earlier designs. Six diverse full passes entered the fixed confirmation; four passed both selected models, with all selected controls passing. The other two remain failures despite alternative samples passing. This is an 8% initial computational pass rate, not a measured hit rate. See [initial screen](../apical_expansion_fixed/RESULTS.md) and [confirmation](../expansion_validation/RESULTS.md).

**Histidines:** Twenty mutant/comparator sequences, two parents and three controls received five samples each; five selected variants then received two fresh seeds with five samples per seed. D51H is consistently disruptive in the predictor, while T70H retains substantial contact recall but misses confidence. That supports preserving D51 in the main search; it does not prove the salt bridge's causal binding energy. Neither unmodified parent passes both selected models in this confirmation, so a variant-level pass does not establish improved affinity. See [panel](../his_panel/RESULTS.md) and [confirmation](../his_confirmation/RESULTS.md).

**pKa direction:** A41H's free/bound pKas are 6.406/5.500 and 7.302/6.183, respectively. Both shifts are negative, and the independent-His Kd(6.0)/Kd(7.4) proxy is 2.48 and 4.91. One free pKa exceeds the frozen 6–7 window; the strict both-pose pass count remains zero. Its ring contacts mainly target Leu212 and, in one pose, Asn215; direct Arg208 repulsion is not established. L7H/L82H has inconsistent direction and is deprioritized. All saved structures, pKa tables and warning outputs were audited. The proxy is conditional on fixed conformations and is not an experimental affinity-ratio prediction. See [direction evidence and A41H review bundle](../his_direction/RESULTS.md).

**Mouse comparison:** The known human positive works on the crystal target but fails on the AlphaFold human target. Consequently, low scores on AlphaFold mouse cannot be interpreted as species-specific nonbinding. Mouse binding remains unknown for every candidate; there is no validated mouse-positive control in this panel. The comparison records a target-conformation limitation rather than a cross-reactivity result. See [species report](../species_exploration/RESULTS.md).

**Uniform preparation:** The same frozen fixed-receptor, restrained-binder protocol is applied to controls, successes and failures. Execution and displacement failures are retained. Raw decisions are never overwritten by relaxed geometry. The historical parent close pair can be resolved with small backbone movement, but that does not retroactively qualify its raw failed sample. Per-round preparation records retain all results.

## Recommended next work

1. Keep A41H as the leading *pH-direction hypothesis*, with the failed strict heuristic explicit. Prioritize long_092, short_018 and short_026 as new parental scaffolds, while retaining short_015 for diversity with its weaker consistency documented.
2. Apply the existing audited contact/site-selection procedure to the new parents and run a bounded histidine-variant panel before another large generation batch. Preserve candidate-specific anchors; evaluate plausible peripheral and partially buried sites rather than treating periphery as a universal rule. Independently predict mutants before computing free/bound pKa shifts.
3. Freeze the final evidence tiers, sequence/fold diversity and assay-tag/context checks before assembling the submission. Twenty sequences is a maximum, not a quota; use clearly labeled exploratory entries only when justified. None of the present designs is yet an experimentally verified binder or pH switch.
4. Keep the conserved-patch and nanobody arms paused for this campaign. Do not repeat the inconclusive AlphaFold mouse screen without a useful mouse-positive control and target-conformation calibration. No upstream bug report is needed to progress the submission.

No additional cloud generation follows automatically from this report. No upstream issue, pull request or challenge submission was filed.

## Completion and resource record

All 278 planned local preparations completed successfully. 219 meet the frozen displacement bounds; the other 59 prepared models are excluded from diagnostic interpretation, with raw evidence preserved. All ten retained lead models meet the movement bounds. The full [evidence audit](evidence_audit.json) verifies launch/cloud hashes, archived analysis-source hashes, local inputs and the ten-model review bundle.

The gross Modal meter reached $9.96 before credits; the latest snapshot is $9.92. The non-monotonic meter may lag, so the higher observed value is retained for budgeting. The $12 allocation for this round is a conservative reservation, not an additional $12 charge. All cloud apps are verified stopped, no local preparation remains running, and 52 tests pass. All work remains within the $30 total authorization. The first preflight failure and corrected expansion share the original $6 expansion allocation.
