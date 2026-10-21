# Working method narrative for the TfR1 challenge

This is a draft methods artifact, not a submitted application or a claim of experimental binding.

## Target and generation

The campaign uses the bound human TfR1 apical conformation from 6WRW, native residues 190–383, with explicit mapping between native and crop numbering. Generated designs are checked against the intact receptor dimer and observed glycans while preserving the modeled bound crop during context restoration. Observed glycans do not constitute an exhaustive glycan ensemble.

The new bounded expansion contains 100 miniproteins of 55–72 residues and 100 of 80–110 residues, conditioned on native residues 208, 210 and 211. It uses the pinned BoltzGen revision `247b9bbd8b68a60aba854c2968d6a0cddd21ad6d`. Binding and exclusion labels occupy separate parser entries. The configured and actual generated feature arrays are audited, because correct labels in a YAML file alone do not prove they reached the model. Generation seeds, sequence-design seeds, input hashes and source snapshots are archived.

The initial expansion preflight stopped before generation: the wrapper mistakenly equated dataset entries with design count after the CLI switched to ten diffusion samples per entry. The corrected check multiplies batch count by diffusion samples. The failed attempt and correction share the original $6 allocation; no upstream bug report or pull request was filed.

## Structural screening

Sequence-conditioned refolding supplies only receptor templates, excludes binder templates, and uses the previously calibrated no-MSA configuration. 2DS25.5 is a positive structural control; a sequence shuffle and unrelated VHH are presumed negative controls, not experimentally demonstrated nonbinders. Confidence scores are not affinity measurements.

Target integrity is checked using all 194 target Cα atoms, with a maximum 2 Å RMSD to the supplied receptor conformation. In multisample runs, select the highest confidence representative among eligible samples before inspecting binder geometry. Preserve every sample and every rejection.

The unchanged candidate screen combines interface confidence, binder fold, receptor-frame interface pose, generated-reference backbone contact recall, intended hotspot and region contacts, full-context receptor fit, and zero counted atom-pair clashes below 2 Å. Generated-reference agreement demonstrates self-consistency rather than recovery of an experimentally observed pose. New generation is screened once; it is not automatically promoted to a confirmed parent. Sequence diversity filtering excludes exact previous sequences and near-duplicates of the earlier shortlisted leads.

## Histidine engineering

Both existing apical parents enter the diagnostic panel. The main substitution search requires peripheral-contact/SASA eligibility in two selected parent poses and excludes possible oppositely charged anchors within 3.5 Å. These are geometric heuristics, not energetic alanine-scanning results. D51H and T70H are explicit diagnostic exceptions, and the failed D51H/T70H pair is retained as a comparator.

The first panel evaluates 20 variant/comparator sequences, two parents and three controls, five samples each. A predeclared shortlist contains at most three passing variants per parent. Five variants, both parents and controls then receive five samples at each of two fresh seeds. This procedure does not select alternate seeds or replace failed representatives after viewing geometry.

## Preparation and pKa direction

A single preparation protocol fixes the receptor and restrains the binder backbone in Amber14/GBn2 minimization. Execution, atom displacement and geometry are reported separately from raw decisions. Prepared geometry never retroactively turns a raw-model failure into a pass.

Only confirmed variant interfaces undergo PROPKA direction triage. The actual independently predicted mutant structures are used; histidines are not inserted afterward into an assumed parental interface. Free and bound calculations retain the same binder conformation. Introduced-His pKa shifts and an independent-site proton-linkage proxy compare pH 6.0 with 7.4. Receptor cropping, artificial termini, rotamer/conformation uncertainty and omitted titratable-site coupling limit interpretation. The proxy is not an experimental affinity-ratio prediction.

The strict direction rule is preserved even when a candidate offers useful evidence outside a heuristic boundary. Short-parent A41H retains a consistent release-direction hypothesis despite failing the free-pKa 6–7 window in one pose. It is labeled accordingly rather than silently promoted to a strict pass.

## Submission interpretation

Prioritize diverse, traceable computational hypotheses; twenty is a maximum, not a quota. Report human binding, acidic release and mouse cross-reactivity as unvalidated unless measured experimentally. Human-to-mouse sequence conservation, proximity to a basic residue, and a favorable pKa calculation do not establish these properties. The local map flags species-variable residues near A41H and an alignment caveat. An audited matched human-crystal / human-AlphaFold / mouse-AlphaFold comparison uses actual species sequences and identical binders. Its AlphaFold-human positive control fails, so mouse binding remains unknown. No mouse classification is inferred from the low confidence scores.
