# Apical expansion: 16 complete first-screen passes

All 200 new miniproteins and 203 independent structural predictions completed. All three matched controls pass. The 200 generated feature arrays contain the intended hotspot labels, and all sequences are unique within the new batch. Exact-identity comparison against all 295 earlier generated sequences finds no duplicates; this is not a comprehensive sequence-database novelty search.

| Arm | Designs | Geometry eligible | Complete one-sample raw pass |
|---|---:|---:|---:|
| Apical, 55–72 residues | 100 | 22 | 13 |
| Apical, 80–110 residues | 100 | 15 | 3 |
| Total | 200 | 37 | 16 |

This is an 8% initial computational pass rate, not an experimental hit rate. The short arm is more productive in this batch, but two limited campaigns do not establish a universal scaffold ranking. New seeds and the CLI's ten-sample generation batching differ from the previous 40-per-arm configuration.

The published positive control has interface ipTM 0.632; presumed negatives score 0.138 and 0.135. All candidate raw decisions use the unchanged fold, interface, contact, context and clash criteria plus the prospectively specified target-integrity filter. High-confidence off-pose designs remain failures. Prepared geometry cannot replace raw failures.

The general geometry shortlist contains 30 diverse sequences, including lower-confidence hypotheses. A separate predeclared follow-up selects only complete raw passes, with a maximum of six and round-robin representation of length strata. Existing 80% sequence-identity diversity filtering is applied, including comparison to the prior seven shortlisted leads.

## Fixed fresh-seed follow-up

The six selected sequences are `expansion_apical_long_092`, `expansion_apical_short_026`, `expansion_apical_long_089`, `expansion_apical_short_015`, `expansion_apical_long_020`, and `expansion_apical_short_018`. With three controls, seeds 1409 and 1511 and five samples per seed, this produces 90 predictions. Its additional $1.50 cap comes from the already reserved validation budget within the $30 total. No further generation follows automatically.

Every one of the 203 selected initial-screen models is also processed with the same frozen local preparation diagnostic, including controls and failures. This is recorded separately, and no raw selection is replaced.

## Execution and budget

The first run, `apical-015-6wrw`, stopped during preflight after 22.7 worker seconds and before generating any designs. At 100 designs per arm, the pinned CLI defaults to ten diffusion samples per dataset item. The wrapper incorrectly asserted that dataset length itself equaled final design count. The corrected worker validates dataset length × diffusion samples and retains feature-level checks.

The corrected `apical-017-6wrw` run completed in 1,963 worker seconds. Inputs and frozen analysis-source hashes are verified. A conservative $0.666 allowance for the failed attempt is deducted from the original $6 allocation; the corrected attempt's bound fits within the remaining $5.334. The combined allocation was not increased. The post-screen gross compute meter is approximately $9.79 for the project, before credits; the later validation has a separate reservation.

Both failed and successful run records remain available. No upstream issue, pull request or challenge submission was filed. None of these one-sample passes establishes binding, pH selectivity, mouse cross-reactivity or submission readiness.

## Uniform preparation completion

All 203 initial-screen candidate/control models completed preparation; 203 executed successfully, 155 satisfy the frozen movement bounds, and 48 prepared models are excluded by those bounds. Every raw decision remains unchanged. This diagnostic does not rescue failed interfaces or overwrite the initial ranking. Full per-model displacement and geometry records are in `preparation/results.json`.
