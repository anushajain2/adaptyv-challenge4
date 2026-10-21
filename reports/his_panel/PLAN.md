# Bounded His panel and confirmation

The first screen contains 25 sequences: two unmodified parents, three matched controls, and 20 mutant sequences. Five samples are generated at seed 1009 for each sequence. This is a structural-retention screen, not an affinity or pH-binding prediction.

Additional sites must satisfy the existing peripheral-contact/SASA heuristic in both selected parent poses. Positions making a possible oppositely charged contact within 3.5 Å in either pose are excluded from the main search. The five highest-ranked eligible sites per parent generate singles and up to four pairs with Cα separation 5–16 Å in both poses. These rules are heuristics, not an energetic ranking. D51H, T70H, and D51H/T70H on the long parent are explicitly included as diagnostics/comparator; Asp51 is otherwise preserved. The full audit is in `data/his_panel/contact_audit.json`.

The unchanged target-integrity filter selects one representative before inspecting binder geometry. Apply the original raw-model gates and preserve every sample. Apply one frozen structural preparation to every selected candidate and control; do not use prepared geometry to replace original gate decisions.

If controls pass, select at most three passing variants from each parent, ranked by interface ipTM, excluding the already-failed double-mutant comparator. Confirm that shortlist, both parents, and all controls at seeds 1103 and 1201, five samples per seed. Parents remain included irrespective of their initial screen outcome. No alternate-seed rescue is allowed.

Only variants passing both confirmation seeds receive direction triage on their independently predicted interfaces. The pKa criteria were frozen separately in `direction_plan.json` before panel results were available. The full first-screen and confirmation allocation is at most $3, including conservative startup allowances. A separate 200-design apical expansion has a $6 allocation within the $30 total.
