# Similarity review of two new computational parents

**Review complete within the stated database limits.** Both retained parents have zero RCSB sequence hits and zero Swiss-Prot close-sequence flags under the frozen settings. Foldseek returned structural matches; the followed-up matches are weak or partial, and do not establish whole-protein duplication. This is not a novel-fold verdict.

All five binder-only structure queries completed against both databases, and every retrieval control passed. The published 2DS25.5 control recovers 6WRW; the local Swiss-Prot search recovers human P02786 exactly. No new NCBI job or Modal compute was used.

## Sequence results

| Design | RCSB hits | Swiss-Prot HSPs | Smallest E-value | Close flags | Masked intervals |
|---|---:|---:|---:|---:|---|
| diversity_short_055 | 0 | 7 | 2.93 | 0 | `[]` |
| diversity_long_037 | 0 | 4 | 3.21 | 0 | `[{"from": 42, "to": 54}, {"from": 55, "to": 73}]` |

The close-sequence review heuristic is identity ≥50%, query coverage ≥80% and E ≤0.001 per HSP. It is not a competition cutoff. BLAST+ 2.17.0 uses BLOSUM62, word size 3, SEG soft masking, composition adjustment 2 and up to 100 targets against the verified UniProtKB/Swiss-Prot canonical 2026_03 collection (575,748 sequences). Isoforms from the separate supplement are excluded. Masking and short/compositionally biased queries limit sensitivity.

## Structural results

The table reports the largest query-normalized TM score among the first ten server-ranked hits, separately for each pose/database. Target-normalized TM and coverage prevent a partial match to a large protein from appearing to be a whole-protein duplicate.

| Query | Database | Target | Query TM | Target TM | Query coverage | Target coverage |
|---|---|---|---:|---:|---:|---:|
| diversity_short_055_seed2309 | afdb-swissprot | AF-P54498-F1-model_v6 | 0.554 | 0.140 | 89.1% | 16.1% |
| diversity_short_055_seed2309 | pdb100 | 4hxg-assembly1.cif.gz_F | 0.601 | 0.085 | 89.1% | 9.3% |
| diversity_short_055_seed2411 | afdb-swissprot | AF-Q9M6E8-F1-model_v6 | 0.561 | 0.086 | 92.2% | 9.6% |
| diversity_short_055_seed2411 | pdb100 | 4hxg-assembly1.cif.gz_F | 0.607 | 0.085 | 89.1% | 9.3% |
| diversity_long_037_seed2309 | afdb-swissprot | AF-A6W1C5-F1-model_v6 | 0.513 | 0.295 | 76.8% | 38.0% |
| diversity_long_037_seed2309 | pdb100 | 1ua8-assembly1.cif.gz_A | 0.496 | 0.309 | 73.7% | 40.2% |
| diversity_long_037_seed2411 | afdb-swissprot | AF-O29261-F1-model_v6 | 0.527 | 0.305 | 81.1% | 39.7% |
| diversity_long_037_seed2411 | pdb100 | 8otz-assembly1.cif.gz_Y | 0.547 | 0.145 | 78.9% | 16.5% |

For short_055, the largest PDB query score is about 0.607 against part of 4HXG chain F, but target-normalized TM is only about 0.085 and target coverage about 9.3%. This is a small-region match. Long_037 has only modest followed-up scores; its largest PDB query score is about 0.547, with target TM about 0.145. No followed-up candidate match has both length-normalized TM scores above 0.5. These observations are limited to the retrieved top-ten hits, not exhaustive alignments over all known proteins.

The server advertises PDB100 **20240101** and AlphaFold/Swiss-Prot **v6**. The old structural snapshot, finite search sensitivity, limited follow-up and predicted AlphaFold structures prevent broad novelty claims. A separate comparison of the confirmed poses against the current draft is archived in [confirmed_diversity.json](../diversity_confirmation/confirmed_diversity.json); it tests redundancy within this panel, not global fold novelty.

## Provenance and decision

The [frozen plan](plan.json), [full analysis](analysis.json), raw/gzipped server responses, local BLAST JSON and hashes preserve the review. [Anonymous public-byte verification](public_provenance_check.json) preceded API upload. Queries contain binder atoms only. Both selected poses are retained. No sequences were changed to improve these results.

The review identifies no close-sequence flag requiring exclusion of either new parent. Retain long_037 and short_055 as additional computational binding hypotheses, with their preparation checks and weaker/stronger confidence differences explicit. Neither contains histidine; pH release, mouse binding, tag compatibility and experimental affinity remain unestablished.

Method references: [Foldseek](https://github.com/steineggerlab/foldseek), [RCSB search API](https://search.rcsb.org/), [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), [UniProt release files](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).
