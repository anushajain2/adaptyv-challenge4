# New-parent histidine confirmation

Completed 140 predictions for 14 sequences. Selected structural controls pass. 4/6 new variants pass every selected raw check across 2 seed(s). These are computational structural results, not experimental binding or pH-switching measurements.

| Sequence | Selected ipTM | Passing samples | Selected raw decision |
|---|---|---:|---|
| expansion_apical_short_018 | 0.628 / 0.633 | 6/10 | 1801: clash_clearance |
| expansion_apical_short_026 | 0.602 / 0.613 | 9/10 | Pass |
| expansion_apical_long_092 | 0.663 / 0.668 | 8/10 | Pass |
| expansion_apical_long_092_R3H | 0.657 / 0.650 | 8/10 | Pass |
| expansion_apical_long_092_E53H | 0.656 / 0.659 | 6/10 | 1801: interface_pose |
| expansion_apical_short_018_K3H | 0.618 / 0.622 | 4/10 | 1801: clash_clearance |
| expansion_apical_short_018_E57H_K60H | 0.641 / 0.631 | 9/10 | Pass |
| expansion_apical_short_026_E19H | 0.585 / 0.594 | 8/10 | Pass |
| expansion_apical_short_026_V60H | 0.591 / 0.585 | 9/10 | Pass |
| apical_short_008 | 0.612 / 0.621 | 7/10 | Pass |
| apical_short_008_A41H | 0.620 / 0.624 | 8/10 | Pass |

Target-integrity rejects: 12/140. Passing unselected samples never replace a failed representative. Every sample and failed check is retained in [results.json](results.json).

## Matched parents

- expansion_apical_long_092: passes both selected raw checks. Variant-only success does not demonstrate improved affinity.
- expansion_apical_short_018: does not pass both selected raw checks. Variant-only success does not demonstrate improved affinity.
- expansion_apical_short_026: passes both selected raw checks. Variant-only success does not demonstrate improved affinity.

Only variants passing both selected raw checks, with controls passing, proceed to the separately frozen preparation and PROPKA checks. A41H is an earlier computational comparison, not a measured pH-positive control.

## Uniform preparation
Completed 28/28 planned selected-model preparations; 28 executions succeeded and 26 meet the movement limits. Excluded prepared models remain recorded. Prepared evidence never changes raw calls.

Each phase reserves at most $2 within the $4 round allocation and $30 gross project cap. No new generation, mouse inference, upstream issue/PR or challenge submission is performed.
