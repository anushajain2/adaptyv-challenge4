# Direction check: one computational hypothesis

**Subsequent result:** independent refolding rejected D51H/T70H in all ten samples. The pKa evidence below is preserved as a historical conditional calculation, not a current advancement recommendation. See the [mutant-refolding report](../mutant_validation/RESULTS.md). Do not treat the review bundle's earlier direction hypothesis as a submission-ready switch.

The confirmed computational parent `apical_long_015` (98 residues) was tested with two His pairs. Each pair was rebuilt on both independently predicted parent poses (seeds 607 and 709), with two local rebuild seeds per pose (11 and 29): **eight local models total**.

**D51H/T70H passes the predeclared direction criteria in all four models.** D51H/K62H fails because the second histidine's pKa shift is too small. No extra pairs were scanned to obtain a positive result.

| Pair | ΔpKa at position 51, bound minus free | ΔpKa at second position | Passes all four models? |
|---|---|---|---|
| D51H/T70H | −0.54 to −0.56 | −1.72 to −2.04 | Yes |
| D51H/K62H | −0.53 to −0.57 | −0.017 to −0.031 | No; second shift fails −0.5 criterion |

For the passing pair, H51 has free pKa 6.41–6.50 and bound pKa 5.85–5.94. H70 has free pKa 7.09–7.12 and bound pKa 5.05–5.41. The independent-His, fixed-conformation calculation gives a Kd(pH 6.0)/Kd(pH 7.4) factor of 14.0–16.6. **That factor is a model-derived triage number, not a measured or calibrated affinity prediction.** It omits other titrating sites, coupled protonation, conformational changes and the energetic cost of the mutations.

All four passing-pair models have zero counted full-context <2 Å clashes. Engineered histidines lie 33.7–35.3 Å from the artificial target terminal residues, beyond the predeclared 10 Å buffer. Both pairs lower both predicted pKas in all models, but sign alone is insufficient: the rejected pair remains rejected despite a calculated factor near 2.

## Evidence checks

All eight models were audited against their saved PDB structures and free/bound pKa tables. Mutant sequences match the requested substitutions; actual histidine heavy atoms are complete. All saved pKa values and shifts match the report. PROPKA emitted no warnings. Recomputed full-context clash counts match the reported zeros. The largest binder backbone pose shift during restrained local relaxation is 0.1805 Å relative to its predicted parent. This small movement confirms the calculation largely holds the parent pose fixed; it does not prove the mutant would retain that pose freely.

A relative output-path error interrupted report bookkeeping after the first two completed models. The path handling was repaired, those exact models were reused, and the remaining six were computed once. `resume_record.json` records the original and resumed code hashes; no scientific threshold or pair selection changed. `evidence_audit.json` records the structure and pKa checks and file hashes.

## What is ready, and what is not

The [review bundle](review_bundle/manifest.json) contains the parent and proposed variant sequences, two independently predicted parent complexes, four locally mutated/relaxed variant complexes, and their free/bound pKa evidence. `direction_hypotheses.fasta` contains only D51H/T70H. The files label the mutant models as locally rebuilt structures, **not independent mutant predictions**.

This is one provisional computational parent and one direction hypothesis. Neither is experimentally validated. There is no evidence yet of mouse binding, preserved mutant affinity at pH 7.4, or measured acidic release. These variants are related sequences, not independent design diversity.

**Best next step:** independently refold the D51H/T70H sequence, alongside its unmodified parent and the same controls, with fresh seeds and the frozen target-filter/structural criteria. Do not supply the mutant binder pose as a template. If it preserves the intended fold, interface and receptor clearance, retain it as a prioritized submission hypothesis with these limitations; then address mouse cross-reactivity separately. Do not scale another broad batch or fill 20 submission slots with near-duplicate variants merely to reach the maximum.

The His scans ran locally and added no Modal GPU job. The latest project gross meter is approximately $7.79 before credits, within the $30 cap. No upstream issue, PR or challenge submission was filed.
