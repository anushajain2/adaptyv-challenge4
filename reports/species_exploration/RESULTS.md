# Human/mouse exploration is inconclusive because the matched human conformation control fails

All 81 predictions completed: nine binder sequences on three receptor inputs, with one seed and three samples per input. The human crystal target reproduces the positive-control interface. The same known human binder fails on the AlphaFold human target. Therefore the low scores on the AlphaFold mouse target cannot support a mouse-nonbinding conclusion, and no cross-reactivity claim is made.

| Sequence | Human crystal ipTM | Human AlphaFold ipTM | Mouse AlphaFold ipTM |
|---|---:|---:|---:|
| expansion_apical_short_015 | 0.522 | 0.145 | 0.125 |
| expansion_apical_short_018 | 0.631 | 0.103 | 0.116 |
| expansion_apical_short_026 | 0.601 | 0.102 | No eligible sample |
| expansion_apical_long_092 | 0.668 | 0.117 | 0.124 |
| apical_short_008 | 0.622 | 0.150 | 0.103 |
| apical_short_008_A41H | 0.618 | 0.154 | 0.114 |
| 2DS25.5, positive on human | 0.643 | 0.350 | 0.381 |
| Scrambled, presumed negative on human | 0.127 | 0.134 | 0.148 |
| Unrelated VHH, presumed negative on human | 0.118 | 0.137 | 0.170 |

Scores are from the frozen confidence selection after the 2 Å target-integrity filter. Twenty-two of 81 samples fail that filter; all three mouse short_026 samples fail, so that case has no selected representative. All six human-control checks were required for an interpretable matched comparison: five pass and the human AlphaFold positive fails. On the crystal target, 2DS25.5 recovers 92.9% of native contacts with interface RMSD 0.89 Å and contact PAE 5.50 Å. On AlphaFold human, its selected model recovers 0% of native contacts, with interface RMSD 37.95 Å and contact PAE 12.92 Å. This invalidates a simple species-only interpretation.

## Matched preparation and limitations

The human and mouse AlphaFold monomers use the actual P02786 and Q62351 sequences and audited aligned apical crops, respectively human 190–383 and mouse 192–385. Both crops contain 194 residues. The same binder sequence and reference backbone are carried across all three conditions; binder coordinates are excluded from prediction templates. All inputs, masks, chain sequences, exported PAE dimensions, launch hashes and the prelaunch analysis-source hash were checked. The [mouse mapping audit](../mouse_context_audit/RESULTS.md) records local sequence-alignment uncertainty rather than silently rewriting the earlier conservation map.

The human crystal comparison uses 6WRW. AlphaFold coordinates come from the official [human P02786](https://alphafold.ebi.ac.uk/entry/P02786) and [mouse Q62351](https://alphafold.ebi.ac.uk/entry/Q62351) entries. This is a comparison of target sequence and model conformation together, not an isolated test of species specificity. The core can align well while the binding surface or its response during prediction differs.

There is no validated mouse-positive control in this panel. The mouse reference binder pose is transferred from human, not experimentally observed. Full mouse dimer/glycan screening and mouse pKa calculations were not performed. The human-specific fixed-context relaxation routine was not applied to mouse with an inappropriate residue map. All species scores are raw, advisory evidence; the previous raw human qualification decisions remain unchanged.

## Decision

Retain the four confirmed human parents and the qualified A41H direction hypothesis. Leave mouse cross-reactivity **unknown**. Do not spend the remaining budget repeating this failed AlphaFold comparison. A useful future mouse comparison needs a validated mouse-binding control and a suitable target-conformation model; it is secondary to completing human pH-direction triage for this submission.

Run `species-019-af` completed in 291.4 worker seconds inside its $1.50 reservation. Complete selected metrics are in [comparison.csv](comparison.csv), all 81 samples in [results.json](results.json), and verified input/runtime records in [run_record](run_record). No upstream issue, PR or challenge submission was filed.
