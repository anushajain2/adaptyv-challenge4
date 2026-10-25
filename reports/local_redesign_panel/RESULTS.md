# Local redesign comparison results

Completed 130 predictions across 13 sequences, with five samples for each of two fresh seeds. The raw control panel passed: **True**. Strict prepared control panel passed: **False**. The separate interface-movement sensitivity panel passed: **False**.

**0 of 6 new sequences qualifies for the preregistered follow-up.** The frozen eight-design release remains unchanged; nothing has been submitted.

The [interpretation note](INTERPRETATION.md) explains why A41H with A34K/A38R remains a promising reserve despite missing the full preset rule. It separates that judgment from the unchanged formal selection outcome.

The strict and sensitivity checks remain separate. Diagnostic pKa values below are shown even when a raw or preparation gate fails; those values cannot qualify a design. Two prediction seeds measure repeatability in the same predictor, not independent model or experimental validation.

| Sequence | Raw samples passing /10 | Both selected raw pass | Strict eligible poses /2 | pKa shifts (seed order) | His-only ratios (seed order) |
|---|---:|---|---:|---|---|
| apical_short_008_A41H | 5 | False | 0 | -0.891, -0.827 | 3.675, 3.180 |
| apical_short_008_A41H_mpnn01 | 6 | False | 1 | -1.570, -0.640 | 7.433, 2.554 |
| apical_short_008_A41H_mpnn02 | 7 | False | 1 | -1.369, -1.178 | 3.030, 5.220 |
| apical_short_008_A41H_mpnn03 | 8 | True | 2 | -1.042, -1.046 | 4.477, 4.451 |
| diversity_long_037_S90H | 6 | True | 2 | -0.344, -0.624 | 1.492, 1.798 |
| diversity_long_037_S90H_mpnn01 | 0 | False | 0 | -0.495, -0.029 | 1.568, 1.040 |
| diversity_long_037_S90H_mpnn02 | 6 | True | 0 | -0.741, -0.457 | 3.003, 1.838 |
| diversity_long_037_S90H_mpnn03 | 7 | True | 2 | -0.591, -1.045 | 1.819, 2.172 |

## Matched original parents

- `apical_short_008`: 10/10 raw samples pass; both selected poses pass: True; 1/2 strict prepared poses eligible.
- `diversity_long_037`: 8/10 raw samples pass; both selected poses pass: True; 1/2 strict prepared poses eligible.

## Exclusions and interpretation

- Control `2ds25_5_6wrw`, seed 3001: raw pass True; strict movement pass False; internal geometry pass True; see the complete preparation/protonation record in `direction.json`.
- Control `scrambled_6wrw`, seed 2909: raw pass True; strict movement pass False; internal geometry pass True; see the complete preparation/protonation record in `direction.json`.
- Control `scrambled_6wrw`, seed 3001: raw pass True; strict movement pass False; internal geometry pass True; see the complete preparation/protonation record in `direction.json`.
- `apical_short_008_A41H_mpnn01`: changes A38K; strict structure-and-direction in both poses False; matched His baseline eligible False; original parent eligible False; baseline-comparison precondition False (this alone does not establish improvement).
- `apical_short_008_A41H_mpnn02`: changes A34E, A38E; strict structure-and-direction in both poses False; matched His baseline eligible False; original parent eligible False; baseline-comparison precondition False (this alone does not establish improvement).
- `apical_short_008_A41H_mpnn03`: changes A34K, A38R; strict structure-and-direction in both poses False; matched His baseline eligible False; original parent eligible False; baseline-comparison precondition False (this alone does not establish improvement).
- `diversity_long_037_S90H_mpnn01`: changes T6S; strict structure-and-direction in both poses False; matched His baseline eligible True; original parent eligible False; baseline-comparison precondition True (this alone does not establish improvement).
- `diversity_long_037_S90H_mpnn02`: changes V57E; strict structure-and-direction in both poses False; matched His baseline eligible True; original parent eligible False; baseline-comparison precondition True (this alone does not establish improvement).
- `diversity_long_037_S90H_mpnn03`: changes V88I; strict structure-and-direction in both poses False; matched His baseline eligible True; original parent eligible False; baseline-comparison precondition True (this alone does not establish improvement).

ProteinMPNN compatibility scores are not protonation-dependent binding energies. The direction calculations use the same bound conformation with and without the receptor; they omit an unbound ensemble and full proton linkage. Histidine neighbour distances are descriptive, not proof of hydrogen bonds or electrostatic repulsion. No experimental release, affinity, mouse cross-reactivity or tag compatibility is established.

The [frozen plan](../local_redesign/plan.json), [sampling records](../local_redesign/mpnn/results.json), [raw results](results.json), [direction and neighbour audit](direction.json), [PROPKA term decomposition](mechanism_audit.json), [evidence audit](evidence_audit.json), and [selected models plus all confidence samples](review_bundle/manifest.json) retain failures and provenance.

Observed gross cost increase: **$0.322**. Conservative cumulative gross accounting: **$13.042/$30**. The round cap is $6; the prediction panel reservation was $2. No credit subtraction is used.

The conditional tag follow-up was skipped because no redesign met the frozen improvement rule.
