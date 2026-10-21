# Broad screen: 21 October 2026

All 200 unique designs were generated and scored. Two pass the complete one-sample screen, and seven pass every geometry criterion. All seven eligible designs are apical miniproteins. Fresh-seed evaluation is complete: two leads reproduce their candidate-level passes, but the batch is held after one selected positive-control model deforms its receptor. See [the fresh-seed report](../broad_validation/RESULTS.md). These are computational leads, not confirmed binders or submission-ready pH switches.

| Arm | Designs | Geometry eligible | Complete one-sample pass | Advanced |
|---|---:|---:|---:|---:|
| Apical, 55–72 residues | 40 | 5 | 1 | 5 |
| Conserved, 55–72 residues | 40 | 0 | 0 | 0 |
| Apical, 80–110 residues | 40 | 2 | 1 | 2 |
| Conserved, 80–110 residues | 40 | 0 | 0 | 0 |
| Apical nanobody | 40 | 0 | 0 | 0 |

The complete initial passes are `apical_short_008` (interface ipTM 0.620, contact PAE 5.79 Å) and `apical_long_015` (0.554, 7.19 Å). The other five eligible leads have lower confidence and enter fresh-seed evaluation under the predeclared geometry-based shortlist rule, not a relaxed final gate. The shortlist contains seven sequences after diversity filtering; it was not padded to 30.

The published 2DS25.5 positive control passes: interface ipTM 0.639, native contact recall 96.4%, binder self-fit Cα RMSD 0.615 Å, and receptor-fitted interface Cα RMSD 0.888 Å. Both presumed negative controls remain below the confidence threshold. These limited controls support this screening configuration; they do not establish a general classifier or a calibrated nanobody threshold.

## Conserved-patch result

Actual generated feature files contain both intended hotspot labels for all 80 conserved-arm designs. Nevertheless, none contacts both hotspots after refolding, and none satisfies the intended-region contact fraction. Seven contact one hotspot; 73 contact neither. The highest-scoring example has ipTM 0.670 but binds a different predicted region. Its high score is not evidence of a conserved-patch hit.

A separate audit of the generated binder backbones finds no contacts to either conserved hotspot at the existing contact cutoff. That backbone-only observation excludes placeholder inverse-folded sidechains and does not assert that full sidechains could never reach a hotspot. It does show that hotspot labels in the input were insufficient to enforce the intended placement in this batch. See `conserved_placement_audit.json`. The failed arm remains reported as failed; its off-target predictions are not relabeled as prospective successes.

## Execution and cost

The original Modal run completed all generation and inverse folding, then stopped with a generic cancellation after 72 of 203 predictions. The logs do not establish the cause. A separate run completed only the missing 131 predictions. Both raw runs are preserved; a new local directory joins their disjoint outputs without replacing predictions. The continuation uses the same nominal seed 307 in a fresh process, so its random trajectory differs from an uninterrupted run. Per-design source provenance is included in `results.json` and `run_record/execution_provenance.json`.

All 203 predictions have complete finite PAE and pLDDT arrays, matching sequences and chain assignments. The analyzer reconstructs actual predicted sidechains and checks coordinate ordering. All 200 designs have the intended generated hotspot features. The local test suite passes 36 tests.

The post-screen Modal meter reports **$7.58348721 gross project usage**, before credits. Billing may settle later. The interrupted run and continuation fit within the original $6 batch reservation using conservative resource bounds. Fresh-seed evaluation reserves at most $4 more within the approved $30 gross project cap. No upstream issue or PR has been filed.

## Next decision

All seven leads and three controls were evaluated with seeds 401 and 503, five samples per seed. Each seed's representative was selected by confidence before inspecting its geometry. A miniprotein must pass the complete structural screen in both selected representatives, with controls passing, to become a provisional computational parent. The batch-control condition failed; no histidine engineering or further generation followed. The fresh-seed report proposes a bounded confirmation test for target-integrity-aware sample selection. Conservation alone does not establish mouse binding; PROPKA alone does not establish pH-dependent affinity.

Files: `results.json` contains all metrics and provenance; `ranked_designs.csv` and `all_designs.fasta` cover all 200 designs; `shortlist_for_evaluation.fasta` contains the seven leads; `screening.png` and `screening.pdf` visualize the complete screen.
