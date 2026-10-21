# Confirmation: one parent advances

The prospective 50-prediction confirmation completed with **all six selected control checks passing**. `apical_long_015` passes every candidate criterion in both fresh seeds and is a provisional computational parent. `apical_short_008` remains a backup because one selected model fails the unchanged full-receptor fit cutoff. This supersedes the practical hold on the long lead without rewriting the previous failed experiment.

| Candidate | Length | Selected interface ipTM, seeds 607 / 709 | Contact PAE, Å | Passing samples | Decision |
|---|---:|---|---|---:|---|
| apical_long_015 | 98 | 0.542 / 0.556 | 7.56 / 7.37 | 9/10 | Advance to bounded direction check |
| apical_short_008 | 58 | 0.624 / 0.624 | 6.12 / 6.01 | 9/10 | Backup; selected seed-607 full-context fit fails |

The long lead's selected binder self-fit RMSDs are 0.414 and 0.506 Å; receptor-frame interface RMSDs are 1.232 and 1.156 Å. Both selected models contact the intended hotspots, retain at least 80% of the generated backbone-reference contacts, and have no counted clashes with the modeled target, omitted receptor, partner subunit or observed glycans. These are consistency checks against a generated reference, not experimental pose recovery.

For the shorter lead, the target fits the supplied bound crop at 1.800 Å in seed 607, satisfying the new selection filter. Its fit to the separate full-receptor context is **2.0457 Å**, above the existing 2.0 Å validation limit. These are different comparisons. The cutoff was not rounded or relaxed, and another sample was not substituted after viewing binder validation. This lead remains useful as a backup, not a passing parent from this experiment.

## Selection and controls

The predeclared target-conformation filter rejected **1 of 50 samples**: scrambled control, seed 709, sample 2, with target Cα RMSD 2.526 Å. That was also its original confidence-selected sample. The best remaining eligible sample, index 1, was selected by the unchanged confidence formula and remains negative (ipTM 0.131). No novel or positive-control selection changed. Every case had at least one eligible sample.

The selected positive-control ipTM values are 0.652 and 0.642. All four selected presumed negatives score between 0.126 and 0.131. Both selected positives pass the unchanged fold, receptor-frame interface and contact-recall criteria. The prospective result supports continuing this bounded workflow, but does not prove the filter universally fixes confidence ranking: in these fresh seeds, neither positive control needed a replacement.

All uploaded inputs and archived analysis-source hashes match the prelaunch records. All 50 samples, rejected models and original selections are retained in `results.json`. The previous failure remains in `reports/broad_validation/`; none of its metrics or decisions were changed. `confirmation.png` and `confirmation.pdf` show all samples and selections. The test suite passes **42 tests**.

## Direction check and remaining limitations

Only the long parent advanced. The bounded His-pair experiment tests at most two pairs on both selected prediction poses, using two local rebuild seeds per pose. Results and any direction hypotheses are recorded separately under `reports/confirmation_direction/`.

The parent and its variants still lack experimental binding evidence. The target-conformation filter is a heuristic and can miss local problems; observed glycans are incomplete. Neither conservation nor these human-target models establish mouse binding. Mutated sequences require independent fold/interface evaluation before being proposed as stronger submission candidates. PROPKA shifts alone do not establish pH-dependent binding or retained affinity at pH 7.4.

## Cost and records

The run used 207.4 seconds of measured worker time and stopped normally. Its reservation was $1.50 within the existing $30 gross project cap. The immediate post-run gross meter high-water is approximately **$7.79**, before credits; billing may reconcile. Cloud state and billing snapshots are in `run_record/`. Subsequent direction checks run locally.

No upstream issue, pull request or challenge submission was filed.
