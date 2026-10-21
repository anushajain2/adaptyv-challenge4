# Fresh confirmation of the apical diversity round

Five initial full-screen candidates and three controls received five predictions at each of fresh seeds 2309 and 2411 (80 total). All six selected control checks pass. The frozen target-integrity-aware confidence selection was applied before interface evaluation; all samples and failed selected poses remain archived.

| Candidate | Samples passing / 10 | Both selected raw models pass | Selected failure |
|---|---:|---|---|
| diversity_short_004 | 8 | False | seed 2411: clash_clearance |
| diversity_short_023 | 8 | False | seed 2411: clash_clearance |
| diversity_short_041 | 3 | False | seed 2309: clash_clearance |
| diversity_short_055 | 9 | True | None |
| diversity_long_037 | 8 | True | None |

**diversity_long_037 and diversity_short_055 are the two retained computational parents.** Both selected poses of each pass the uniform v2 preparation movement and internal-geometry checks. Preparation completed for all 16 selected models; 14 meet movement guards, with the two scrambled-control poses excluded. Raw clash failures in the other three candidates remain failures even when relaxation improves geometry.

The four retained poses remain below the frozen reciprocal-TM redundancy threshold against every current-draft selected pose and each other; full comparisons are in [confirmed_diversity.json](confirmed_diversity.json). This supports distinctness within this panel, not novel folds or evolutionary unrelatedness. [The review bundle](review_bundle/manifest.json) preserves every selected raw/prepared model and its checks, including the failed candidates and controls.

Short_055 explicitly models a C31–C42 disulfide: raw S–S distances 2.003/1.870 Å and prepared distances 2.039/2.036 Å. OpenMM infers the bond, making oxidation an explicit modeling assumption. Adaptyv [documents support for disulfide-containing proteins](https://docs.adaptyvbio.com/wiki/technologies/cell-free-protein-expression), but formation of this particular bond is unmeasured. Long_037 has an NTS sequon at 24–26 and mildly positive GRAVY (+0.155), which are descriptive review flags rather than established expression failures. See the [disulfide audit](disulfide_audit.json).

Neither retained parent contains histidine. Their pH-dependent release, mouse cross-reactivity and assay-tag compatibility are unestablished. Binding is unmeasured. Both originate from the same BoltzGen method family as the existing draft. Long_037 has the stronger interface-confidence signal; short_055 is a lower-confidence but reproducible additional parent.

The [database review](../diversity_novelty/RESULTS.md) is complete: both sequence services and Foldseek completed with passing retrieval controls, no close-sequence flags and weak/partial followed-up structural matches. Both retained parents have been added to the [eight-sequence draft](../submission_draft/README.md). The timestamped `summary.json` preserves the pre-admission confirmation state; later admission is recorded separately in `submission_update.json`. No submission has been made. No extra cloud generation, mutation panel, upstream issue or PR follows automatically.

The gross Modal meter after the two bounded phases is approximately $12.24 total. All cloud apps are stopped. This round used approximately $1.48 above the prior recorded high-water mark, within its $4 reservation and the $30 gross authorization; meter readings can lag.
