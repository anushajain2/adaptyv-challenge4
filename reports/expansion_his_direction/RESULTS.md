# A41H remains the leading pH-direction hypothesis; no new variant clears the full check

The bounded follow-up is complete: 24 new His variants on three parents, 160 first-screen predictions, 140 fresh-seed confirmation predictions, and 60 uniform selected-model preparations. Fifteen new variants pass the first selected-model screen. Four of six shortlisted variants pass both fresh selected models, with all structural controls passing. A41H also reproduces its structural pass and pKa-lowering direction in two additional poses.

**Zero variants pass the strict direction rule in both new poses.** A41H again fails only the free-pKa 6–7 heuristic in one pose. Four other model calculations are excluded because PROPKA reports incomplete/protonation-failed receptor groups. Those are scoring failures, not evidence of adverse pH behavior or nonbinding.

## Results that should guide the next decision

| Variant | Fresh confirmation raw passes / 10 | Selected raw checks in both seeds | Direction evidence | Interpretation |
|---|---:|---|---|---|
| long_092 R3H | 8 | Pass; matched parent also passes | Opposite signs across poses: +0.259 / −0.643 | Deprioritize as a release hypothesis |
| long_092 E53H | 6 | Fail: interface pose in seed 1801 | Not evaluated | Structural failure remains |
| short_018 K3H | 4 | Fail: clash in seed 1801 | Not evaluated | Structural failure remains |
| short_018 E57H/K60H | 9 | Pass; matched parent fails one selected clash check | Both calculations excluded by receptor protonation warnings | Structurally retained, pH direction unknown |
| short_026 E19H | 8 | Pass; matched parent also passes | −0.854 in one usable pose; second excluded | Secondary hypothesis, not replicated pH evidence |
| short_026 V60H | 9 | Pass; matched parent also passes | −0.084 in one usable pose; second excluded | Little direction signal; terminal-context concern |
| Existing short_008 A41H | 8 | Pass; matched parent also passes raw checks | −1.104 / −1.051 | Leading qualified hypothesis; one strict heuristic failure |

Names in this table abbreviate the full `expansion_apical_*` identifiers. See the [confirmation report](../expansion_his_confirmation/RESULTS.md) for exact names, all samples, selected scores and failures. All ten mutant/comparison models admitted to direction scoring passed the unchanged prepared geometry checks. Preparation movement exceptions are recorded separately and never replace raw decisions.

## Free/bound pKa results

| Variant / seed | Free His pKa | Bound His pKa | Bound − free | His-only Kd(6.0)/Kd(7.4) proxy | Strict model rule |
|---|---:|---:|---:|---:|---|
| long_092 R3H / 1801 | 7.042 | 7.301 | +0.259 | 0.71 | Fail |
| long_092 R3H / 1901 | 6.261 | 5.618 | −0.643 | 1.89 | Fail: proxy below 2 |
| short_026 E19H / 1801 | 7.200 | 6.346 | −0.854 | 3.49 | Fail: free pKa above 7 |
| short_026 V60H / 1801 | 6.673 | 6.589 | −0.084 | 1.14 | Fail |
| Existing A41H / 1801 | 7.235 | 6.131 | −1.104 | 4.84 | Fail: free pKa above 7 |
| Existing A41H / 1901 | 6.373 | 5.322 | −1.051 | 2.56 | Pass |

Invalid calculations are deliberately absent from this numerical table. Their structures, raw tables and warnings remain available. The strict-pass FASTA is empty, and no criterion was changed after inspecting results.

Together with the [earlier direction round](../his_direction/RESULTS.md), A41H now has **four independently predicted, accepted poses with negative pKa shifts**, approximately −0.91 to −1.12. The conditional His-only proxy spans 2.48–4.91. Two of those four poses miss the free-pKa window. This window remains a heuristic; its failures are not proof that the design cannot release. Conversely, repeated agreement in one predictor and PROPKA does not establish actual binding, affinity or a pH switch. A41H passes 17/20 samples across the two fresh-seed confirmation rounds; that descriptive total does not erase any earlier batch's failed representative or parent comparison.

## Why four calculations were not scored

All ten variant/comparison models passed the raw structural screen, the preparation movement limits, and the original prepared pose/contact/context criteria. Six yielded warning-free free/bound calculations. The other four stopped on receptor groups:

- short_018 E57H/K60H, seed 1801: receptor Arg339 (crop residue 150).
- short_018 E57H/K60H, seed 1901: receptor Glu369 and Ser370 (crop residues 180/181).
- short_026 E19H, seed 1901: receptor Asn317 (crop residue 128).
- short_026 V60H, seed 1901: receptor Arg339 (crop residue 150).

A [read-only diagnostic](protonation_warning_diagnosis.json) confirms that expected heavy atoms are present at all warning residues. It also finds abnormal intra-receptor nonbonded distances in some affected models, including a 1.08 Å Glu369 C–OE2 separation in the failed double-mutant model. The frozen preparation fixes receptor coordinates, so such internal receptor geometry is not repaired by binder minimization. Our existing clash screen evaluates binder/receptor/context clashes, not every intra-receptor stereochemical defect. Good receptor Cα fit therefore does not guarantee suitable geometry for protonation. This is a preparation limitation, not proof of the cause of every warning: not every warning group has an associated nonbonded pair below the diagnostic 1.8 Å cutoff.

No receptor atoms were moved selectively, no warnings were ignored, and no alternative rotamers or samples were tried. In particular, E19H's excluded second pose remains unscored even though the first pose is encouraging. The required next repair should be specified uniformly across successes, failures, parents and controls, with original evidence preserved.

## Terminal and assay context

The [context audit](context_audit.json) examines all ten saved complexes. short_026 His19 is about 25–27 Å from binder termini, and its usable pose contacts target Lys374 at 4.36 Å. That is consistent with a plausible local electrostatic hypothesis, not proof of a mechanism. Its one usable free pKa is above 7, and there is no second accepted score.

By comparison, short_026 His60 is only about 4.9–6.2 Å from the binder C-terminal carboxyl atoms. The short_018 double mutant's His60 is about 3.6–6.8 Å from that terminus. long_092 His3 is within roughly 9–10 Å of binder termini. These sites need explicit assay-construct/tag sensitivity analysis before prioritization. A41H remains farther from the C terminus (about 17 Å), but its tag compatibility is also untested.

Free calculations use the bound binder conformation. The receptor is cropped, glycan coverage is incomplete, and the proxy accounts only for introduced histidines rather than full proton linkage and conformational ensembles. Mouse cross-reactivity remains unknown; this round made no additional mouse predictions.

## Recommended next step

Keep **A41H as the main pH-direction hypothesis**, and keep **short_026 E19H as a secondary, incomplete hypothesis**. Preserve the diverse unmodified parent shortlist for the submission. Do not promote the other new variants merely because their structural confidence is high, and do not fill twenty slots with near-duplicate mutants.

Before further generation or another broad His scan, specify and validate one uniform receptor-preparation/protonation procedure that addresses internal receptor geometry while bounding changes to the interface. Apply it equally to candidates, parents and controls; retain these original failures as the primary record and label any revised calculations as a separate analysis. Pair that with assay-tag/context review and then assemble an evidence-tiered submission shortlist. No experimental binding or acidic-release claim is currently justified.

## Verification, artifacts and budget

All 300 prediction samples and their selection decisions are retained. All 60 selected-model preparations executed; 54 meet the frozen movement bounds. Four first-screen variants, one later unmodified parent model and one later negative-control model exceed movement bounds; their prepared evidence is excluded. The scored mutant models have no PROPKA warnings, and their saved atoms, pKa tables, shifts and strict flags have been recomputed and checked by the [evidence audit](evidence_audit.json). The [round audit](round_evidence_audit.json) also checks input/source hashes, the complete sample inventory and the frozen six-variant selection.

The [review bundle](review_bundle/manifest.json) includes the five structurally retained variants/comparison, matched parents, raw checks, structures, valid and invalid pKa evidence, and explicit qualification flags. The [first-screen plan](../expansion_his_panel/PLAN.md) and both run records preserve the prospective rules. The only scoring-implementation clarification, recorded before any first-screen result inspection and before confirmation, validates preparation coverage against actually selected models rather than requiring models for cases with no eligible sample; no threshold changed.

Both cloud runs completed, all apps are verified stopped, and the latest recorded gross Modal usage is **$10.45 before credits** within the **$30 total cap**. Billing may lag. The round's $4 reservation is a bound, not an additional $4 charge. Local verification: **56 tests pass**. No upstream issue, PR or challenge submission was filed.
