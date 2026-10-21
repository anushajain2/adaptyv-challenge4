# Broad campaign: 21 October 2026

The user authorized the 200-design campaign after reviewing the staged recommendation. This is a new search after the hotspot-conditioning correction, not a rerun of the original 75 designs. No upstream issue or PR is authorized.

- 80 apical miniproteins: 40 of 55–72 residues and 40 of 80–110 residues.
- 80 conserved-patch miniproteins, with the same length allocations.
- 40 apical nanobodies with redesigned CDRs on the existing 7EOW framework, configured for 120–130 residues.

Target: bound 6WRW apical crop, native residues 190–383. Hotspots are 208/210/211 for apical arms and 327/331 for conserved arms. The archived full-dimer accessibility audit on the bound crop passes all five positions. The conserved arm uses native 326–347 as its intended contact region, rather than the published apical binder's footprint. Exclusion masks and positive hotspots use separate parser entries. Actual features are checked before generation and again in generated files.

One fresh seed (307), one sample per design, no MSA, target-only templates, full PAE export. Two hundred novel designs plus three matched controls make 203 screening predictions. Sequence diversity and equal-weight within-class percentile ranks of interface ipTM, contact-pair PAE, hotspot fraction and fold consistency order the leads. Geometry-eligible candidates enter an at-most-30 sequence-diverse shortlist with representation across arms. The miniprotein confidence threshold is descriptive for nanobodies; there is no positive VHH control on this target. A one-sample pass is not a qualified parent.

The first batch reserves $6 gross, within a $16 generation/screening allocation; $7 is reserved for deeper evaluation and $3 for contingency. The project cap remains $30 gross including credits. One L40S, no automatic retries, 6,600-second work bound, 6,750-second function timeout and 600-second startup bound. Launch checks live billing, prior high-water billing, other apps and input hashes. Expansion requires a cost/quality checkpoint after this batch.

If there are geometry-eligible leads and controls pass, evaluate a diverse shortlist with fresh seeds. Do not claim mouse binding from conservation, or pH switching from PROPKA. Keep published controls out of any novel-design submission. Preserve original failed results without changing thresholds retroactively.

Diversity clarification, recorded during generation before refolding: use the actual generated design mask to extract redesigned VHH CDR residues for the 80% identity comparison; use full miniprotein sequences. Always remove exact full-sequence duplicates. This prevents a shared framework from collapsing distinct CDR designs.
