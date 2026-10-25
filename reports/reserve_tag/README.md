# Fresh confirmation and tag sensitivity of the MPNN reserve

The user authorized this follow-up on 26 October 2026. It is a new experiment on
an unqualified reserve, not the automatic follow-up skipped by the previous round.
Nothing is being submitted and the frozen eight-design fallback remains unchanged.

The [frozen plan](plan.json) compares four related sequences and three controls:
the original short_008 parent, A41H, the A41H/A34K/A38R reserve, its matched H41A
ablation, known binder 2DS25.5, scrambled control, and unrelated VHH. Each is modeled
bare and with the same 64-residue platform-default suffix, in fresh seeds
3509/3607 with five samples per seed: **140 predictions**.

The suffix is the previously archived [Adaptyv platform default](https://docs.adaptyvbio.com/wiki/technologies/construct-tag-design),
including linkers, GFP11 and Twin-Strep. The exact challenge-specific construct is
not independently confirmed. The documentation could not be freshly retrieved
with the web tool on this turn; this experiment explicitly uses the archived
default, whose sequence and provenance are retained in the plan. Tagged sequences
are modeling inputs only, never submission strings.

Sample selection and geometry thresholds are unchanged from the earlier tag
comparison. Confidence is restricted to the binder core so the suffix cannot
dilute the metric. The full binder and suffix are excluded from templates.
Raw control separation is checked independently in each condition. Every selected
full construct receives the same numerical preparation, including failures.

The primary readout is descriptive replication of the reserve's pKa shift across
fresh poses and tag conditions, with structural eligibility and controls reported.
Free pKas, shifts and the His-only release proxy remain continuous readouts. The
original 6–7 free-pKa band and all other legacy qualification checks are also
reported explicitly; they are not replaced or tuned to admit the reserve.
The earlier interface-movement sensitivity rule is retained as a separate
diagnostic and cannot automatically promote a candidate.

PROPKA uses the complete suffix chemistry for tagged free/bound forms. Only the
designed His41 contributes to the reported His-only proxy. The H41A comparator
tests modeled structure/pose preservation; its lack of a designed His score is
not evidence of absent experimental pH sensitivity. The tag and other titratable
residues remain present. Neither the predictor nor this fixed-conformation pKa
calculation measures binding, pH release, or capture-surface compatibility.

The run has a **$2 gross cap** within the existing **$30 total**, with one L40S
container, no automatic retries, and bounded runtime. The launch ledger keeps a
$1 meter-lag allowance above the higher of the prior conservative accounting and
the current gross meter. Credits are not deducted from spending.

Reproduction entry points are `prepare_reserve_tag.py`, `launch_reserve_tag.py`,
`analyze_reserve_tag.py`, `prepare_reserve_tag_models.py`, and
`analyze_reserve_tag_direction.py`. Completed directories are never overwritten.
The run record archives frozen inputs and code. The original parent/reference
coordinates are copied into the input `references` directory and hashed.
