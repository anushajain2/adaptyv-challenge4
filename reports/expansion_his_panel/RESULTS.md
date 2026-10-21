# New-parent histidine panel

Completed 160 predictions for 32 sequences. Selected structural controls pass. 15/24 new variants pass every selected raw check across 1 seed(s). These are computational structural results, not experimental binding or pH-switching measurements.

| Sequence | Selected ipTM | Passing samples | Selected raw decision |
|---|---|---:|---|
| expansion_apical_short_018 | 0.624 | 1/5 | Pass |
| expansion_apical_short_026 | 0.616 | 4/5 | Pass |
| expansion_apical_long_092 | 0.668 | 4/5 | 1709: clash_clearance |
| expansion_apical_long_092_E43H | 0.649 | 3/5 | Pass |
| expansion_apical_long_092_R35H | 0.656 | 4/5 | Pass |
| expansion_apical_long_092_R3H | 0.660 | 4/5 | Pass |
| expansion_apical_long_092_N26H | 0.667 | 3/5 | 1709: clash_clearance |
| expansion_apical_long_092_E53H | 0.656 | 5/5 | Pass |
| expansion_apical_long_092_R3H_R35H | 0.643 | 5/5 | Pass |
| expansion_apical_long_092_N26H_R35H | 0.653 | 0/5 | 1709: clash_clearance |
| expansion_apical_long_092_R35H_E53H | 0.653 | 5/5 | Pass |
| expansion_apical_short_018_K3H | 0.623 | 3/5 | Pass |
| expansion_apical_short_018_E4H | 0.613 | 3/5 | 1709: interface_pose |
| expansion_apical_short_018_K60H | 0.614 | 4/5 | Pass |
| expansion_apical_short_018_E57H | 0.649 | 3/5 | 1709: clash_clearance |
| expansion_apical_short_018_I32H | 0.550 | 4/5 | Pass |
| expansion_apical_short_018_K3H_I32H | 0.543 | 3/5 | 1709: clash_clearance |
| expansion_apical_short_018_E4H_I32H | 0.550 | 3/5 | Pass |
| expansion_apical_short_018_E57H_K60H | 0.641 | 5/5 | Pass |
| expansion_apical_short_026_E15H | 0.581 | 5/5 | Pass |
| expansion_apical_short_026_E19H | 0.592 | 3/5 | Pass |
| expansion_apical_short_026_V60H | 0.595 | 4/5 | Pass |
| expansion_apical_short_026_L18H | 0.477 | 0/5 | 1709: interface_confidence |
| expansion_apical_short_026_L25H | 0.404 | 0/5 | 1709: interface_confidence |
| expansion_apical_short_026_E15H_E19H | 0.519 | 3/5 | Pass |
| expansion_apical_short_026_E15H_L18H | 0.431 | 0/5 | 1709: interface_confidence |
| expansion_apical_short_026_E15H_L25H | 0.351 | 0/5 | 1709: interface_confidence |
| apical_short_008 | 0.615 | 2/5 | 1709: interface_pose |
| apical_short_008_A41H | 0.625 | 5/5 | Pass |

Target-integrity rejects: 11/160. Passing unselected samples never replace a failed representative. Every sample and failed check is retained in [results.json](results.json).

Frozen follow-up selects: expansion_apical_long_092_R3H, expansion_apical_long_092_E53H, expansion_apical_short_018_E57H_K60H, expansion_apical_short_018_K3H, expansion_apical_short_026_V60H, expansion_apical_short_026_E19H.
The unchanged parents and A41H comparison sequences are retained in confirmation. At most two new variants per parent advance; ranking uses selected interface ipTM, not prepared geometry or pKa.

## Uniform preparation
Completed 32/32 planned selected-model preparations; 32 executions succeeded and 28 meet the movement limits. Excluded prepared models remain recorded. Prepared evidence never changes raw calls.

Each phase reserves at most $2 within the $4 round allocation and $30 gross project cap. No new generation, mouse inference, upstream issue/PR or challenge submission is performed.
