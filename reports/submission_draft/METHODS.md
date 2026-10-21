# Methods accompanying the draft

## Scope and provenance

These are computational hypotheses for human TfR1 binding and, for A41H, an acidic-release direction. No binding, pH-selectivity, mouse-cross-reactivity or expression measurements are available. The six ranked sequences represent five generated parent lineages; A41H is a substitution in our own generated short_008 parent. Published binders were calibration controls, not submitted sequences or starting binder templates.

The retained designs originate from corrected BoltzGen apical campaigns against the bound 6WRW receptor crop, native residues 190–383. Native/crop mapping is explicit. Full-dimer and observed-glycan context is restored for steric auditing while retaining the modeled bound crop. Observed glycans are not a complete glycan ensemble. The crop omits native position 142; its genotype correction is handled in full-target preparation, not an invented residue within the crop.

BoltzGen revision `247b9bbd8b68a60aba854c2968d6a0cddd21ad6d` and Adaptyv-derived Modal wrappers underpin generation/refolding. Binding and exclusion labels are separate parser entries, and actual feature tokens are checked. The earliest 75-design pilot had zero positive hotspot labels and is not evidence for successful hotspot conditioning. Corrected miniprotein generation uses native hotspots 208, 210 and 211. Conserved-patch and nanobody branches did not yield retained parents. See [generation methods](../next_round/METHODS.md) and [campaign results](../RESULTS.md).

## Structural evidence

Refolding templates contain the receptor only; binder coordinates are excluded. This campaign uses a single-sequence/no-MSA configuration supported by its own limited control comparison, not a general claim that MSAs harm prediction. The published 2DS25.5 complex is the structural positive. Scrambled sequence and unrelated VHH controls are presumed negatives, not experimentally established nonbinders. Chai-1 did not provide independent interface confirmation.

Each round defines model selection before geometry inspection. Historical confirmation uses the frozen confidence-based selection after a target-integrity check. The latest bare/tagged comparison selects the minimum mean target/core cross-PAE among target-integrity-eligible samples; core PAE replaces whole-chain ipTM to avoid suffix dilution. These rules are reported separately, without pooling their pass rates into one success estimate.

Structural screens cover binder fold, receptor-frame interface pose, contacts to the generated design reference, intended hotspots, context fit and counted clashes. Reference-contact recall for a novel design is self-consistency, not experimentally observed interface recovery. Selected failed poses remain failures; favorable alternative samples do not replace them.

## Direction and construct context

Histidine variants are independently predicted before pKa evaluation. Free and bound PROPKA calculations use the same binder conformation, and the introduced-His-only proton-linkage ratio is a conditional mechanism proxy rather than an affinity prediction. The frozen heuristics require shifts at most −0.5, a free His pKa within 6–7, ratio proxy at least 2, and separation from artificial receptor termini. Preparation warnings and movement failures exclude interpretation; minimization never rescues a raw structural failure.

The original preparation fixed receptor coordinates. A separately frozen v2 protocol restrained both backbones while allowing receptor sidechains and hydrogen coordinates to move, applying the same procedure to parents, variants and controls. Revised calculations on old predictions are sensitivity analyses. Some control/movement checks still fail, and zero qualified two-pose pH switches were obtained. All of these exclusions are preserved.

A matched bare/tagged panel tests the published 64-residue platform-default linker/split-GFP11/Twin-Strep suffix as an explicit challenge-construct assumption. Core scoring excludes suffix contacts; full-construct pKa retains the suffix where modeled. Flexible suffix placements have low confidence, so tag-independence failures cannot be read as demonstrated experimental interference. Bare sequences alone appear in the draft upload.

## Ranking, diversity and sequence review

Ranking prioritizes the repeated A41H direction hypothesis, then repeatedly retained parent interfaces, additional lineages and a matched parent comparator. It is a transparent editorial choice rather than a trained scoring function. All candidates have unresolved mouse cross-reactivity; the matched human-AlphaFold positive failed in the mouse exploration, preventing an interpretable mouse classification.

Sequence descriptors use installed Biopython, whose version is recorded in the audit. Motif scans and composition flags are descriptive only and did not change any sequence. Global alignment uses match 2, mismatch −1, gap-open −4 and gap-extension −0.5, with identity divided by the longer sequence. Backbone RMSD is measured after binder-only alignment on sequence-matched Cα pairs, with both coverage fractions recorded. It is distinct from receptor-fitted interface RMSD used in the structural screen. The first selected model in each latest bare round is used even when it fails a screen. No exhaustive known-protein similarity or structural-novelty search has yet been performed.

## Reproducibility and assistance

The repository contains source snapshots, pinned model revision, run configurations, source hashes, selected structures, metrics, raw rejection reasons and preparation evidence. The draft's [manifest](manifest.json) links its exact inputs and exported sequences. The [audit](audit.json) records per-round evidence and package software versions. Different prediction seeds are computational replication within one modeling system, not independent biological validation.

Anusha directed the project using OpenAI Codex for code development, literature/documentation retrieval, execution and analysis. BoltzGen, PROPKA, OpenMM and the documented wrappers supplied the relevant scientific computations. This statement does not imply use of Claude for this project or experimental confirmation. The code is available at [anushajain2/adaptyv-challenge4](https://github.com/anushajain2/adaptyv-challenge4). No external submission has been made by this workflow.
