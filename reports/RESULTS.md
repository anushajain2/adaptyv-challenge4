# TfR1 days 1–3: results

The workflow is built and has run on real designs. Target preparation and the 50-design pilot are complete. The direction check found **six exploratory histidine-pair variants with both pKas lowered in both modeled conformations**. **No variant qualifies as a reliable binder candidate**, because every parent design failed the interface-confidence screen. The stronger candidate milestone remains unmet.

## What was completed

- **Target preparation:** 11QC chain A residues 190–383, a 194-residue crop; native-numbered and sequential-numbered structures; human–mouse alignment and conservation map; four hotspot sets. The older 6WRV Ser142 was corrected to Gly in the comparison ectodomain. Primary 11QC already has Gly142, outside the crop.
- **Initial pilot:** 25 miniproteins and 25 nanobodies at the conserved 326–347 patch, using the Adaptyv-recommended Modal wrapper and pinned BoltzGen source. Both complete pipelines returned successfully.
- **Follow-up:** 25 more miniproteins at the exposed basic apical tip, with additional exclusions for surfaces exposed by cropping. This is a human-first, species-variable patch. This pipeline also completed successfully.
- **Direction check:** eight pairs across four clash-free original-pilot parents, two models per pair. Twelve models scored successfully; four models from two pairs failed PROPKA structure/protonation checks and were left unscored.

## Parent quality

| Batch | Designs | Best interface ipTM | Clash-free | Receptor clashes | Target fit failed |
|---|---:|---:|---:|---:|---:|
| Initial miniproteins | 25 | 0.215 | 2 | 17 | 6 |
| Initial nanobodies | 25 | 0.301 | 2 | 18 | 5 |
| Apical-tip follow-up | 25 | 0.258 | 1 | 20 | 4 |

All **75 designs are readable**, with no import errors. “Clash-free” alone does not establish binding. The provisional parent gate requires interface ipTM ≥0.50, backbone refold RMSD <2.5 Å, and zero full-receptor clashes. The brief's 2.0 Å RMSD standard is stricter; the confidence failure already excludes every design under either RMSD threshold. Upstream rank/retention is not a passing binding prediction.

The full receptor check revealed substantial contacts with omitted receptor regions. The small follow-up did not solve the confidence problem. These results support revisiting target context and design constraints before a large generation campaign; they do not establish that TfR1 is intrinsically undesignable.

## Exploratory pKa evidence

Six pairs lowered both pKas consistently. Three also met the stronger direction criteria: both drops ≥0.5, at least one free pKa between 6 and 7, and an independent-site estimated Kd(6.0)/Kd(7.4) ≥2 in both models.

**Every row below is exploratory and fails the parent-quality gate.** Mutation indices are one-based positions in the full binder sequence, not IMGT numbering. All three variants derive from the original conserved-patch pilot.

| Parent | Mutations | Least negative shifts across seeds | His-only estimated Kd ratio | Parent ipTM |
|---|---|---:|---:|---:|
| Nanobody 03 | N60H,Y62H | -1.17, -1.76 | 4.89–4.98 | 0.120 |
| Nanobody 03 | Y62H,T106H | -1.77, -0.80 | 4.56–4.64 | 0.120 |
| Miniprotein 01 | R44H,Y49H | -0.62, -1.36 | 3.04–3.04 | 0.146 |

The ratios are conditional model estimates, not predicted experimental affinities. Mutant binding strength and folding have not been validated. The six successful pairs' engineered residues are at least 24 Å from the artificial target terminal residues.

PROPKA uses the **intact 194-residue trimmed target** plus the mutant binder. The free binder retains exactly the same coordinates. A graft into the full receptor is used only for clash checks: graft junctions can have invalid peptide geometry, so that graft is not used for pKa. Failed earlier attempts are retained under `outputs/direction_exploratory_initial`; the corrected evidence is under `outputs/direction_initial`. Omitted receptor regions, glycans, artificial termini, coupled titration and conformational changes remain limitations.

## Cost and verification

Modal currently reports approximately **US$2.40 gross usage**, within the authorized **US$30 including credits**. Reported billed amount is **US$0.00**, with credits covering compute. Billing can lag; [the raw record](billing_final.json) preserves the exact values. Both cloud apps have completed; no generation job remains running.

**14 local tests passed.** Both initial cloud specifications and the follow-up specification validated. Downloaded input hashes match the archived inputs. All 12 scored direction models passed an actual-output audit for complete histidines, matching free/bound binder coordinates, valid PROPKA groups and arithmetic consistency.

## Files to use

- [Conservation map](../data/prepared/conservation_map.csv), [hotspot sets](../data/prepared/hotspot_sets.json), [trimmed target PDB](../data/prepared/human_apical.pdb).
- [All 75 designs and metrics](../outputs/analysis/pilot_designs.csv), [parent sequences](../outputs/analysis/pilot_designs.fasta).
- [All eight histidine pairs and evidence paths](pilot/all_histidine_pairs.csv), [six exploratory variant sequences](pilot/exploratory_pka_lowering.fasta).
- [Direction summary](pilot/direction_summary.json), [output audit](initial_pilot/direction_audit.json), [machine-readable project status](day1_3_status.json).
- [Overview plot](pilot/pilot_overview.png), [PDF plot](pilot/pilot_overview.pdf), [reproduction instructions](../README.md).

`pilot/provisional_switches.fasta` is deliberately empty: zero variants pass all gates. The exploratory FASTA is not a submission recommendation. No account action is needed from you for this completed pilot.
