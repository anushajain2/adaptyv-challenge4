# Independent mutant refolding: reject D51H/T70H

**The D51H/T70H mutant fails the independent structural test in both fresh seeds and all ten samples.** The matched controls pass. The earlier favorable PROPKA result therefore remains a conditional calculation on locally rebuilt parent poses; it is not sufficient evidence to advance this variant.

| Sequence | Selected interface ipTM, seeds 809 / 907 | Contact PAE, Å | Selected reference-contact recall | Passing samples |
|---|---|---|---|---:|
| Unmodified apical_long_015 | 0.555 / 0.551 | 7.47 / 7.44 | 77.8% / 77.8% | 6/10 |
| D51H/T70H mutant | 0.108 / 0.103 | 25.40 / 24.42 | 0% / 0% | 0/10 |

Both sequences retain their generated fold reasonably well. All mutant binder self-fit Cα RMSDs are 0.67–1.14 Å, but the selected interface poses are displaced by **44.5 and 49.8 Å** in the receptor-fitted frame. They contact a different part of the crop and clash extensively with omitted receptor regions when restored to the full context. Across all ten mutant samples, interface ipTM ranges only **0.091–0.116**, and generated-reference contact recall ranges 0–5.6%. The failure is not confined to the selected models or a borderline confidence threshold. This is computational rejection, not an experimental demonstration of nonbinding.

## Parent and control interpretation

The parent maintains interface ipTM 0.535–0.555 and reference-contact recall 72.2–83.3% across all ten samples. Its selected seed-907 model passes every criterion. Its selected seed-809 model fails only the unchanged zero-clash rule: target Lys205 NZ and binder Glu66 OE2 are **1.694 Å** apart. This is a real close pair in the saved prediction, not a numbering or grafting error. The parent therefore also fails the strict both-selected-seeds rule in this particular matched experiment. Its earlier passing results remain preserved, but it should be described as a computational lead with variable raw-model geometry, not a reliably validated binder.

All six selected control checks pass. Positive-control interface ipTM is 0.636/0.641; all four selected presumed negatives remain at 0.124–0.139. The target-conformation filter rejects five of 50 samples and changes the selections for the mutant in both seeds and the positive control in seed 907. Every case retains at least one eligible sample. Rejected mutant samples also have low confidence and do not provide a hidden success.

## What the mutation may have disrupted

In both selected parent models, binder Asp51 lies close to target Arg208: its OD2 atom is 2.31 Å from Arg208 NH2 in seed 809 and approximately 2.98–3.08 Å from Arg208 NE/NH2 in seed 907. This geometry is consistent with a possible favorable charged contact. Replacing Asp51 with histidine could remove that interaction; **the contact audit does not establish causality or an energetic contribution**. T70H could also matter, and the two substitutions could interact. The paired mutant alone cannot distinguish these possibilities.

This exposes the limitation of the previous direction calculation: its restrained local rebuilding preserved a parent interface that independent mutant-sequence prediction did not recover. A favorable pKa shift on an assumed bound pose cannot substitute for evidence that the variant can adopt that pose.

## Recommended next step

Keep the unmodified parent as a lead, and **do not advance or submit D51H/T70H as the current switch candidate**. Run a small diagnostic comparison of **D51H alone and T70H alone**, alongside the parent and controls, to distinguish each substitution's structural effect. These singles are diagnostic probes, not a claim that either is an adequate final pH switch. The Asp51 contact makes preserving that residue a sensible redesign priority unless subsequent evidence contradicts it. If one single substitution retains the interface, search for a second peripheral site that avoids disrupting important parent contacts. Do not spend another broad generation batch or loosen this mutant's gates to recover it.

The parent's close pair also needs a consistently applied structural-preparation check before future strict advancement; simply replacing its failed selected sample or selectively minimizing only that model would not validate the original result. No new prediction panel or relaxation was launched after this outcome.

## Execution, files and cost

Five sequences × two fresh seeds × five samples produced **50 completed predictions**. Prediction used target-only templates, no MSA and the same frozen target filter and structural thresholds. Actual template features exclude binder coordinates. Input and analysis-source hashes, exact predicted sequences, PAE dimensions, atom-coordinate correspondence and chain roles all pass checks. Mutant sequence and backbone-reference tests confirm only D51H/T70H was introduced; no locally relaxed mutant binder was used as a template. The local test suite passes **44 tests**.

The first attempt stopped before prediction after 16.2 worker seconds because the input filename contained a double underscore unsupported by the pinned parser. A filename alias corrected this; all structure and mask bytes remained identical. The successful run took 203.8 seconds. Original inputs, failed logs, aliases and both reservations are retained. Conservative allowances for both attempts sum to less than the original **$1.50** experiment cap. Latest reported gross project usage is **$7.96348782**, before credits, within the approved $30 total; all cloud apps are stopped.

`results.json` retains all 50 samples and selections; `candidates.csv` summarizes both sequences; `parent_contact_audit.json` records the atom-level diagnostic; `mutant_validation.png` and `.pdf` show confidence and target integrity for every sample. `provisional_parents.fasta` is empty for this experiment. `review_models/` contains selected models explicitly labeled as failed or diagnostic evidence, not submission-ready candidates.

No upstream issue, PR or challenge submission was filed. No experimental binding, mouse cross-reactivity or pH-dependent affinity has been established.
