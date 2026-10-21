# Bounded apical diversity round

Status: **generation and initial screening complete** (`diversity-023-6wrw`); the preplanned [two-seed confirmation](../diversity_confirmation/RESULTS.md) and [database review](../diversity_novelty/RESULTS.md) are also complete. Two additional computational parents have entered the eight-sequence draft; nothing has been submitted externally.

All 100 designs and 103 one-sample refolds completed successfully. The batch contains **99 unique sequences**, none exactly matching the 495 earlier designs. All intended hotspot labels were verified in every generated feature file, and all cloud inputs and analysis-source hashes match the launch snapshot.

| Arm | Designs | Geometry eligible | Full raw-screen passes | Distinct candidates sent to confirmation |
|---|---:|---:|---:|---:|
| Short, 55–72 residues | 60 | 13 | 4 | 4 |
| Long, 80–110 residues | 40 | 2 | 1 | 1 |

All three controls pass: 2DS25.5 interface ipTM 0.631, scrambled 0.113, unrelated VHH 0.098. The positive also meets the frozen pose/contact/fold checks and all controls meet target integrity. These validate this computational panel, not experimental binding or assay performance.

All five full-screen passes also meet the frozen diversity checks against the current draft and each other. Their first-screen metrics are:

| Candidate | Interface ipTM | Contact PAE (Å) | Largest reciprocal TM to current draft |
|---|---:|---:|---:|
| diversity_long_037 | 0.566 | 7.04 | 0.349 |
| diversity_short_023 | 0.571 | 6.98 | 0.647 |
| diversity_short_004 | 0.670 | 6.20 | 0.472 |
| diversity_short_041 | 0.574 | 7.63 | 0.474 |
| diversity_short_055 | 0.507 | 8.22 | 0.359 |

The [diversity audit](diversity_audit.json) records every comparison and the exact chosen order. These were one-sample candidates at this screening stage; subsequent confirmation retained two, as reported above. A 5/100 computational pass rate is not a measured binding hit rate. The gross meter after generation is approximately **$12.02 total**, and no Modal app remained active before confirmation was launched.


The prospective [plan](plan.json) specifies 100 additional apical miniproteins: 60 of length 55–72 and 40 of length 80–110. Both arms retain the audited 6WRW receptor crop, native hotspots 208/210/211, corrected positive/negative hotspot encoding and pinned BoltzGen commit. Generation and inverse-folding seeds are fresh. Refolding uses seed 2203, one sample per design, target-only templates and empty MSA. The published 2DS25.5 positive and two presumed negatives are included.

The existing target-integrity, interface confidence, pose, contact and full-context clash gates are unchanged. The initial screen is not a qualified-parent test. Advancement requires all three controls to pass, a full raw-screen pass, no exact match among 495 previous designs, and less than 80% global sequence identity to the current draft and already selected new candidates. Binder-only TM-align additionally excludes redundant models with both length-normalized TM scores at least 0.80 against any current draft selected pose or an already selected new candidate. This is a descriptive diversity heuristic, not proof of different fold classes or novel folds.

At most six candidates, chosen round-robin across length arms in the frozen ranking order, may receive one confirmation at seeds 2309 and 2411, five samples per seed, with the same controls. Confidence selection first applies the existing target-integrity filter; no alternative seeds or substitutions after outcomes. Both selected models must pass. Every selected candidate and control, including failures, then receives the same v2 structural preparation. Raw failures cannot be rescued by relaxation. Database similarity and sequence-liability review are required before any submission change.

The generation/screen reservation is $2.50 and conditional confirmation $1.50, totaling at most $4 within the $30 gross authorization. Maximum calculated resource-plus-miscellaneous allowances are $2.4774 and $1.4292. One L40S per phase, no retries, and hard function timeouts of 3150/1650 seconds bound execution. All previous allocations are retained conservatively; the full conservative accounting ceiling is $29.9635, not actual spend. The highest earlier recorded gross meter is about $10.75; the latest lower reading is not used to free budget. [Launch evidence](run_record/launch_ledger.json).

This round seeks additional distinct computational parents. It does not establish affinity, pH-dependent release, mouse cross-reactivity or assay-tag compatibility. The conserved-patch and nanobody arms remain paused. No upstream issue/PR or challenge submission is part of this round.
