# Fresh-seed evaluation: two reproducible leads, batch held

All 100 planned predictions completed: seven novel leads plus three controls, five samples each at seeds 401 and 503. Two novel designs pass every candidate criterion in both confidence-selected models. **Neither is cleared as a provisional parent because one selected positive-control model fails the frozen batch-control rule.** No histidine scans or additional generation were launched after this result.

| Lead | Length | Selected interface ipTM, seeds 401 / 503 | Contact PAE, Å | Samples passing candidate screen |
|---|---:|---|---|---:|
| apical_short_008 | 58 | 0.619 / 0.623 | 5.95 / 6.04 | 7/10 |
| apical_long_015 | 98 | 0.559 / 0.545 | 7.42 / 7.30 | 8/10 |

Both selected models for each lead retain the generated fold and pose, contact the intended apical hotspots, and have zero counted clashes with modeled target atoms, omitted receptor regions, partner subunit or observed glycans. Their selected binder self-fit RMSDs are 0.41–0.67 Å and receptor-frame interface RMSDs 0.55–1.23 Å. These are self-consistency results against generated references, not experimental binding evidence. The other five leads fail at least the confidence criterion in both selected models.

## Why the batch remains held

All five positive-control samples at seed 401 pass. Four of five pass at seed 503. The confidence-selected sample at seed 503 is the sole failure:

- Interface ipTM: 0.648; native contact recall: 52/56 (92.9%); binder self-fit RMSD: 0.579 Å.
- Target Cα RMSD to the reference: **5.273 Å**; binder interface RMSD after fitting the entire target: **2.553 Å**, above the predeclared 2 Å cutoff.
- Eighteen target residues are displaced by more than 4 Å after the global fit. The largest changes involve native residues 313–326, reaching 26.6 Å.
- A separate diagnostic fit using only the reference target-interface residues gives a binder interface RMSD of **0.833 Å**. This fit uses target atoms, not a superposition of the binder on itself. It explains how high native-contact recall can coexist with failure in the globally fitted frame.

The exported selected CIF agrees with NPZ coordinates, and atom order, chain identities and sequences pass checks. The observed discrepancy is receptor deformation in a predicted sample, not evidence of a numbering or export mismatch. Confidence selection chose that sample even though four other samples for the same control pass. All four selected negative-control models remain below the confidence threshold.

The local-fit diagnostic is **post hoc**. It does not replace the original metric, reverse the batch decision, or establish that this deformation is biologically plausible. The complete ten-sample audit is in `control_deformation_audit.json`.

## Recommended next step

Keep the two leads and perform one small prospective confirmation experiment before histidine engineering. Freeze a target-integrity check before selecting samples by confidence, apply it identically to positives, negatives and novel designs, and declare a case failed if it has no eligible samples. Test the rule using two new seeds with these two leads and the same three controls; retain every sample and the original failed result. Do not select models by their desired binder pose or tune the interface thresholds to make these leads pass.

If the new control test succeeds and both leads remain structurally consistent, perform bounded histidine-pair direction checks. `scripts/check_broad_direction.py` is prepared but deliberately refuses to run against this failed batch. It would compare the same mutations on both selected parent conformations, with two local rebuild seeds per conformation. Its outputs would still be computational direction hypotheses, not measured pH switches.

More broad generation is lower priority now: there are two reproducible apical leads, while the conserved arm and VHH arm produced no geometry-eligible designs. Mouse binding remains untested. Do not fill a 20-design submission with variants merely to meet a count.

## Reproducibility

The raw run is `outputs/broad_validation010_complete/validation-010-6wrw`; the original broad screen and disjoint continuation remain separate. All uploaded input hashes match the launch ledger. `results.json` and `candidates.csv` distinguish candidate-only passes from the failed batch-control decision. `held_leads.fasta` and four confidence-selected structures in `held_lead_models/` are explicitly marked held. `provisional_parents.fasta` is empty.

The Modal app completed and stopped. Snapshots are recorded in `run_record/billing_after.json` and `run_record/billing_final.json`. The conservative observed gross high-water mark is **$7.66787069**. The newest snapshot is slightly lower, so budget accounting retains the high-water mark; billing may reconcile. Credits are not subtracted when checking the approved $30 cap. No upstream bug report, pull request or challenge submission was made.
