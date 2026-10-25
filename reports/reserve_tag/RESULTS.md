# Fresh reserve confirmation and tag sensitivity — 26 October 2026

**The A41H/A34K/A38R reserve replicated its bare-model pKa-lowering signal.**
It passes the core structural screen in 9/10 bare and 10/10 tagged predictions,
and all four confidence-selected reserve models pass the strict preparation checks.
This strengthens the computational case for retaining it as a challenger to A41H.
It does **not** establish an experimentally functional pH switch or clear the full
legacy qualification rule. No candidate was promoted and nothing was submitted.

The frozen [plan](README.md) used seven sequences, two construct conditions, two
fresh seeds (3509/3607), and five samples per seed: **140 predictions**. All completed.
Selection used binder-core cross-PAE with the original target-integrity guard;
pKa results did not influence selection. All 28 selected models received the same
preparation once. All 56 free/bound PROPKA calculations completed without recorded
protonation failure; 23/28 models pass the strict core-backbone movement guard.

## Structural comparison

Counts are computational screen outcomes, not experimental hit rates. The selected
prepared column includes raw eligibility, prepared geometry, movement, and PROPKA
execution, but does not assert full control-panel or pH-direction qualification.

| Sequence | Bare core passes /10 | Tagged core passes /10 | Selected prepared bare /2 | Selected prepared tagged /2 |
| --- | ---: | ---: | ---: | ---: |
| Original short_008 | 10 | 5 | 2 | 2 |
| A41H | 9 | 6 | 2 | 1 |
| **A41H + A34K/A38R reserve** | **9** | **10** | **2** | **2** |
| Reserve H41A ablation, retaining K34/R38 | 8 | 5 | 1 | 1 |

The reserve's selected contact PAE is 6.07/6.25 Å bare and 5.44/5.21 Å tagged.
Its tagged-versus-bare core interface CA RMSD, after receptor fitting, is
1.38/1.29 Å. The bound target conformation and predictor are shared across these
tests; fresh seeds provide within-model replication, not independent validation.

All **12 selected raw control checks pass**. Known-positive core cross-PAE is
13.19–14.27 Å versus 25.48–27.36 Å for the presumed negatives, separated by the
unchanged 19.827851 Å cutoff. All four selected positive controls also pass strict
preparation. The full prepared-control requirement fails in both conditions:
scrambled bare/3509, scrambled tagged/3607, and VHH tagged/3509 exceed movement
limits. The first and third also fail the interface-motion sensitivity rule.
These failures remain exclusions; they do not erase the observed raw separation.

A41H tagged/3607 exceeds the CA movement limit (Asp14, 1.016 Å), so its pKa result
below is diagnostic only. H41A bare/3607 fails strict movement because a remote
receptor Phe316 oxygen moves 1.155 Å; it passes the separately reported interface
sensitivity rule. H41A tagged/3509 fails the raw clash check (Arg208 NH2–Glu49 OE1,
1.933 Å), even though preparation removes that close pair. No raw failure is rescued.

## Histidine readouts

Only core His41 contributes to the His-only ratio. Full suffix chemistry is retained
for tagged free and bound calculations. Ratios are model-derived proxies, **not
measured affinities**. The free form is the extracted binder in the same prepared
conformation, not a sampled unbound ensemble.

| Sequence / condition | Seed | Free pKa | Bound pKa | Bound − free | His-only Kd6/Kd7.4 proxy | Strict structure eligible |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| A41H bare | 3509 | 7.305 | 6.325 | −0.980 | 4.09 | Yes |
| A41H bare | 3607 | 7.240 | 6.518 | −0.722 | 2.86 | Yes |
| **Reserve bare** | **3509** | **7.320** | **6.200** | **−1.121** | **4.92** | **Yes** |
| **Reserve bare** | **3607** | **7.357** | **6.396** | **−0.961** | **3.93** | **Yes** |
| A41H tagged | 3509 | 6.380 | 6.235 | −0.145 | 1.22 | Yes |
| A41H tagged | 3607 | 7.319 | 6.016 | −1.304 | 6.11 | **No** |
| **Reserve tagged** | **3509** | **6.398** | **4.573** | **−1.825** | **3.08** | **Yes** |
| **Reserve tagged** | **3607** | **7.242** | **6.101** | **−1.141** | **5.06** | **Yes** |

The two fresh bare reserve shifts agree in direction and approximate size with
the previous −1.042/−1.046 results. Relative to matched fresh bare A41H poses,
the reserve has a more negative shift in both seeds and a larger ratio proxy.
Two seeds are insufficient to establish a statistically reliable sequence advantage.
All bare free pKas remain above the preset 6–7 heuristic. That heuristic failure
is recorded unchanged; it is not proof the sequence cannot release experimentally.

The tagged reserve retains a negative shift in both poses, but the magnitudes
and free pKas vary. Seed 3509 lacks the positive backbone hydrogen-bond pKa term
seen in the other reserve free forms. Its larger shift also gives a *smaller*
ratio proxy than the matched bare pose. Thus ranking by shift alone is misleading,
and the tag cannot be credited with a reliable improvement.

PROPKA decomposition again attributes most of the reserve shift to its
volume/desolvation term: −1.182/−1.154 pKa units bare, and −1.530/−1.013 tagged.
Bare Coulomb changes are +0.062/+0.194, opposing rather than causing the lowering.
This supports a burial-based explanation within this model; it does not demonstrate
direct repulsion from the added Lys/Arg. These terms are pKa offsets, not binding
energies. The decomposition is part of the same calculation, not independent evidence.

The H41A construct is a structural comparator. Its lack of the engineered His41
score is tautological and does not show absence of experimental pH sensitivity.

## Tag uncertainty and qualification

The documented 64-residue platform-default suffix includes GFP11, linkers and
Twin-Strep. Its exact use in this challenge remains unconfirmed. Both selected
tagged reserve poses place suffix residues near the receptor (6/68 residue pairs),
but their contact PAE is 26.91/27.12 Å and mean suffix pLDDT is 36.16/36.21 out of
100. **No residue contact pair repeats across the two seeds.** The positive
control also has inconsistent suffix contacts. This is uncertain tag placement,
not reliable evidence of interference or compatibility with the assay surface.

The full preset rule admits **zero** new candidates: the reserve misses the
two-seed free-pKa criterion, no tagged candidate is contact-free in both selected
poses, and the full prepared-control panel does not pass. Strict and sensitivity
decisions remain separately archived. No thresholds were changed after results.

## Decision and next research priority

Retain the reserve as a better-supported computational challenger; keep the
frozen eight-sequence submission package unchanged. Do not claim validated binding,
pH release, tag compatibility, or mouse cross-reactivity. Do not launch further
generation merely to improve a favorable confidence score.

If continuing computational work, the next useful question is whether the reserve's
advantage over A41H survives plausible **free-state conformations and His41 side-chain
arrangements**. A bounded comparison with criteria fixed before looking at outcomes
would address a larger remaining uncertainty than more seeds of this same bound
setup. It would still require experimental testing to establish pH-dependent binding.
That additional study has not been launched or charged by this experiment.

## Evidence, checks and spending

The [raw results](results.json), [prepared/direction results](direction.json),
[tag audit](tag_contact_audit.json), and [review bundle](review_bundle/manifest.json)
retain all selected models and all 28 five-sample confidence arrays. The
[run audit](run_record/verification.json) verifies inputs, pinned code, successful
stages, target-only templates and intended MSA settings. The
[evidence audit](evidence_audit.json) verifies complete coverage, suffix retention,
eight independently recomputed reported His comparisons, and unchanged frozen
release/CSV hashes. The full repository test suite passes: **96 tests**.

Observed gross meter increase: **$0.43**, below the $2 run cap. Conservative
cumulative accounting: **$13.47/$30**, leaving **$16.53**. Credits are included.
The later meter snapshot is lower; the unexplained decrease is not treated as
a refund. Both snapshots and the conservative calculation are retained in the
[budget record](run_record/budget_reconciliation.json). All cloud apps are stopped.
No upstream issue/PR or challenge submission was made.
