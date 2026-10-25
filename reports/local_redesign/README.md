# Local redesign around fixed histidines

Submission remains on hold. The user approved a change from substitution-only
panels to local sequence redesign, with at most $6 new gross Modal spending inside
the existing $30 cap. The frozen eight-sequence release is preserved as a fallback.

This experiment asks whether changing the surroundings of A41H in
`apical_short_008` or S90H in `diversity_long_037` can preserve the modeled interface
and improve the release-direction hypothesis. It does not establish binding,
affinity, pH release, tag compatibility, or mouse cross-reactivity.

The [prospective plan](plan.json) precedes sequence sampling and new predictions.
The [movement audit](movement_audit.json) covers all 32 selected models from the
previous panel. The positive control's 1.017 Å maximum displacement is a receptor
carbonyl oxygen 12.016 Å from the binder. S72H's 1.208 Å movement is a binder oxygen
at the interface, 4.797 Å from the receptor. The scrambled control also has a
movement failure at its predicted interface. These are distinct findings; the
original failures remain failures.

The numerical preparation protocol and strict movement limits are unchanged.
A separate sensitivity analysis reports global backbone RMS, all C-alpha movement,
and backbone movement at residues within 6 Å of the other chain. It cannot rescue
historical outcomes or automatically promote new sequences. This distinction was
declared before the new predictions.

## Local sampling and selection

[ProteinMPNN](https://github.com/dauparas/ProteinMPNN) was pinned to commit
`8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`, vanilla `v_48_020`. The complete receptor
crop and fixed binder residues provide sequence/backbone context. ProteinMPNN
does not model histidine protonation or optimize pH-dependent binding energy.

Both independent source poses must place a candidate position within 8 Å of the
histidine ring. Charged residues, existing histidines, Gly/Pro/Cys, the terminal
five residues, and high-contact or salt-contact anchor proxies remain fixed.
This leaves positions 34/36/38/40 around A41H and 6/57/88/89 around S90H.
Anchor status is a geometric proxy, not measured energetic importance.

CPU sampling generated 256 sequences across two poses per family, temperatures
0.2/0.4, and fixed random seeds. New Cys and His are excluded. After deduplication,
21 A41H-family and 11 S90H-family sequences met the preset identity, mutation-count
and sequon restrictions. They were ranked by mean mutable-position negative log
likelihood on both backbones using four common decoding orders. The matched
histidine baseline was scored on the same positions. Likelihood is a compatibility
score, not a binding-energy estimate.

The frozen greedy Hamming-distance rule selected only three per family; the panel
was not padded with near duplicates. All 256 samples, exclusions, scores and
selected sequences are retained in [the sampling record](mpnn/results.json).

| Fixed histidine family | Additional mutations | Ensemble NLL |
|---|---|---:|
| A41H baseline | None | 1.007 |
| A41H mpnn01 | A38K | 0.828 |
| A41H mpnn02 | A34E, A38E | 0.957 |
| A41H mpnn03 | A34K, A38R | 1.076 |
| S90H baseline | None | 0.801 |
| S90H mpnn01 | T6S | 0.757 |
| S90H mpnn02 | V57E | 0.846 |
| S90H mpnn03 | V88I | 1.159 |

The comparison contains these six redesigns, both His baselines, both unmodified
parents, and three controls. Seeds 2909/3001 with five samples each give 130 fresh
predictions, retaining the target-only template, empty MSA, target-integrity
selection, and original raw screen. Every selected model receives the same
preparation. Direction is assessed separately with strict and sensitivity labels,
paired baseline comparisons, and a geometric inventory of histidine neighbours.
At most two qualifying improvements can receive a subsequent bare/tag comparison;
no qualifying improvement means no tag job.

The rationale is informed by [Computational design of pH-sensitive binders](https://www.biorxiv.org/content/10.1101/2025.09.29.678932v1),
which explores sequence redesign around introduced histidines among several
strategies. This small fixed-backbone experiment is not a reproduction of that
campaign. In particular, it lacks its large-scale sampling and experimental tests.

## Reproduction

Install the main project dependencies and `requirements-local-redesign.txt`.
Clone the official ProteinMPNN repository into `vendor/ProteinMPNN` and check out
the pinned commit. The plan records the checkpoint SHA-256. Scripts refuse to
overwrite completed experimental folders; reproduce in a fresh checkout/output
location rather than removing historical results.

1. `scripts/audit_preparation_motion.py`
2. `scripts/prepare_local_redesign.py`
3. `scripts/run_local_mpnn.py`
4. `scripts/build_local_redesign_panel.py`
5. `scripts/launch_local_redesign.py` (paid, budget guarded)
6. `scripts/analyze_his_panel.py` with the downloaded run and new input/report paths
7. `scripts/prepare_local_redesign_models.py`
8. `scripts/analyze_local_redesign_direction.py`

The prediction run's `run_record` archives the exact code and input hashes.
The original CSV and frozen release are not modified by this experiment.
