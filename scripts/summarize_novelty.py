"""Concise review and per-candidate metadata, preserving database/search limits."""
import json
from pathlib import Path
from novelty_search import ROOT, OUT, save, sha, now
from analyze_novelty import human_tfr1_retrieved


def main():
    data=json.loads((OUT/'analysis.json').read_text());plan=json.loads((OUT/'plan.json').read_text())
    verification=json.loads((OUT/'independent_pdb_verification.json').read_text())
    names=list(dict.fromkeys(c['candidate'] for c in plan['structures'] if c['kind']=='candidate'))
    blast={r['query']:r for r in data['swissprot_sequence'] if r.get('status')=='complete'}
    rcsb={r['query']:r for r in data['rcsb_sequence'] if r.get('status')=='complete'}
    summaries=[]
    for name in names:
        structural=[]
        for r in data['foldseek']:
            if r.get('candidate')!=name:continue
            followed=[h for h in r['hits'] if 'local_TM_align' in h]
            best=max(followed,key=lambda h:h['local_TM_align']['query_normalized_TM']) if followed else None
            structural.append({'query':r['query'],'database':r['database'],'hit_count':r['hit_count'],
                               'near_server_result_limit':r['near_server_result_limit'],
                               'best_query_TM_among_top10':best})
        b=blast.get(name)
        summaries.append({'name':name,'rcsb_sequence_hits':rcsb.get(name,{}).get('total_count'),
                          'swissprot_status':'complete' if b else 'pending',
                          'swissprot_hsps':len(b['hits']) if b else None,
                          'swissprot_close_sequence_flags':sum(h['close_sequence_review_flag'] for h in b['hits']) if b else None,
                          'swissprot_min_evalue':min((h['evalue'] for h in b['hits']),default=None) if b else None,
                          'swissprot_query_masking':b.get('masking') if b else None,
                          'structural_search':structural,
                          'original_pdb_verification':[r for r in verification['rows'] if r['candidate']==name]})
    expected=set(plan['blast']['query_ids'])
    complete=expected==set(blast)
    bcontrol=blast.get('human_tfr1_crop_control')
    bpass=bool(bcontrol and human_tfr1_retrieved(bcontrol['hits']))
    status='COMPLETE_WITH_DATABASE_LIMITS' if complete and bpass else 'STRUCTURE_AND_PDB_SEQUENCE_COMPLETE_SWISSPROT_PENDING_OR_FAILED'
    result={'status':status,'as_of':'2026-10-21','candidates':summaries,
            'controls':{'foldseek_6WRW_recovered':data['foldseek_positive_retrieval_pass'],
                        'rcsb_6WRW_recovered':any(h['target']=='6WRW_2' and h['sequence_identity']==1 for h in rcsb['published_2ds25_5_control']['hits']),
                        'rcsb_human_crop_recovered':any(h['sequence_identity']==1 and h['query_coverage']==1 for h in rcsb['human_tfr1_crop_control']['hits']),
                        'swissprot_human_crop_recovered':bpass if complete else None},
            'all_13_structural_queries_complete':len(data['foldseek'])==26 and all(r['status']=='complete' for r in data['foldseek']),
            'swissprot_all_seven_queries_complete':complete,'new_modal_compute_usd':0,
            'swissprot_sequence_source':data.get('swissprot_sequence_source'),
            'swissprot_result_file':data.get('swissprot_result_file'),
            'submission_sequence_or_rank_changed':False,
            'decision':'Keep A41H as the pH-direction hypothesis and long_092/short_026 as binding candidates. Label short_015 as a known-backbone backup. Do not claim novel folds or evolutionary unrelatedness. Do not remove a candidate solely because a common fold is shared.'}
    save(OUT/'summary.json',result)
    if complete:
        flagged=sum(r['swissprot_close_sequence_flags'] for r in summaries)
        btext=f"Swiss-Prot BLASTp completed all seven queries, including retrieval of human TfR1 accession P02786 ({'pass' if bpass else 'FAIL'}). Across the six candidates, **{flagged} alignments meet the predeclared close-sequence review rule** (identity ≥50%, query coverage ≥80%, E-value ≤0.001). Flags apply to individual HSPs, not unions of disjoint alignments. This is a review heuristic, not a competition cutoff. The composition-corrected, low-complexity-masked search is not proof that no related sequence exists."
        if data.get('swissprot_sequence_source')=='local_UniProtKB_SwissProt_canonical':
            lp=json.loads((OUT/'local_swissprot/plan.json').read_text())
            lr=json.loads((OUT/'local_swissprot/run.json').read_text())
            btext+=f"\n\nThe completed evidence comes from **local BLAST+ 2.17.0 against UniProtKB/Swiss-Prot canonical release {lp['release']}** ({lr['database_sequences']:,} sequences; {lr['database_residues']:,} residues). The original NCBI RID remained WAITING after more than an hour; it is still pending and is not described as completed. A [prospective fallback plan](local_swissprot/plan.json) preserved all seven queries and the screening thresholds before the local search. Official archive checksums were verified; commands, SHA-256 hashes, release metadata, logs and raw JSON results are archived in [local_swissprot/](local_swissprot/). The separate isoform supplement is excluded. This local database need not match the remote service's Swiss-Prot snapshot."
        btext+='\n\n| Query | Reported HSPs | Smallest E-value | Close-sequence flags | Reported masked query intervals |\n|---|---:|---:|---:|---|\n'
        for s in summaries:
            e=s['swissprot_min_evalue'];es='—' if e is None else f'{e:.3g}'
            btext+=f"| {s['name']} | {s['swissprot_hsps']} | {es} | {s['swissprot_close_sequence_flags']} | `{json.dumps(s['swissprot_query_masking'],separators=(',',':'))}` |\n"
    else:
        btext='Swiss-Prot BLASTp is still pending or incomplete. Its queued job is not counted as a no-hit search. The completed sequence conclusion below is limited to RCSB experimental PDB sequences.'
    table=[]
    for s in summaries:
        maxima={db:max((r['best_query_TM_among_top10']['local_TM_align']['query_normalized_TM'] for r in s['structural_search'] if r['database']==db and r['best_query_TM_among_top10']),default=None) for db in ['pdb100','afdb-swissprot']}
        fmt=lambda x:'No hit' if x is None else f'{x:.3f}'
        table.append(f"| {s['name']} | {s['rcsb_sequence_hits']} | {fmt(maxima['pdb100'])} | {fmt(maxima['afdb-swissprot'])} |")
    text=f'''# Database similarity review

**The six sequences have no reported RCSB sequence hits under the recorded settings, but sequence distinctness does not imply a novel fold.** The clearest structural caveat is short_015: both poses closely resemble the known protein AF1782 (2OO2). Retain it as a known-backbone binding backup, without a fold-novelty claim. No sequence or submission rank was changed by this review.

Status: **{status}**. Searches cover six distinct sequences from five parent lineages, both selected bare poses per design, and public retrieval controls. This review added **$0 Modal compute**.

## Sequence results

The RCSB sequence service returned HTTP 204/no hits for all six candidates at E-value ≤10 and identity cutoff 0. Both controls worked: the human TfR1 crop was recovered at 100% query identity/coverage, and the designed positive matched 6WRW entity 2 at 100% identity over all 88 residues. Controls validate retrieval in this panel; they do not quantify sensitivity for these short, compositionally biased designs.

{btext}

The exact requests, retrieval timestamps and server responses are archived under [rcsb/](rcsb/) and [blast/](blast/); any completed local fallback is identified above. Current RCSB searches exclude computed models. Swiss-Prot is a curated sequence collection, not all known sequences or unpublished designs. Alignment identity, coverage and significance must be read together; short local similarities alone do not establish copied sequence. Low-complexity masking reduces sensitivity for these compositionally biased designs; a successful receptor control does not establish sensitivity for every candidate.

## Structural results

Foldseek completed all 13 binder-only queries against PDB100 and AlphaFold/Swiss-Prot. The known 2DS25.5 control retrieves 6WRW, with local TM-align scores approximately 0.993 under both length normalizations and full coverage. The input contains no receptor atoms or assay suffix, so receptor matches cannot dominate the search. Failed selected candidate poses remain included.

The table shows the **largest query-normalized TM score among the first ten server-ranked hits**, across the two poses. It is not the best score over the entire database, a novelty probability, or a binding metric.

| Design | RCSB sequence hits | PDB best query TM | AlphaFold/Swiss-Prot best query TM |
|---|---:|---:|---:|
{chr(10).join(table)}

### short_015: a real known-backbone similarity

Independent comparison to the original [2OO2 coordinates](https://www.rcsb.org/structure/2OO2), author chain A, reproduces the match in both poses. All **61/61 design residues** align to **61/76 modeled AF1782 residues**. Query-normalized TM scores are **0.858/0.858**; target-normalized scores are **0.709/0.710**; TM-align-reported backbone RMSDs are **1.450/1.467 Å**. Global sequence identity is **27.6%**, using the longer modeled sequence as denominator.

![Short_015 and the known AF1782 backbone](short015_known_fold.png)

This supports a shared backbone arrangement, not sequence copying, evolutionary homology, the same interface, or TfR1 binding by AF1782. Foldseek's corresponding E-values are weak (1.437/2.454), illustrating why a simple short helical fold can have close geometry without strong homology evidence. The design was generated de novo in this campaign; de novo generation does not require an unprecedented fold. Competition eligibility remains the organizers' decision.

### Other lineages and partial matches

- **long_092** has moderate similarity to known alpha/beta structures. Against [CcmK2, 3SSQ chain A](https://www.rcsb.org/structure/3SSQ), original-coordinate TM scores are 0.597/0.587 normalized by the design and 0.546/0.538 normalized by the modeled known chain. About 86.5% of the design and 76.2% of the known chain align. This is not evidence of an unprecedented topology.
- **short_026** matches part of [metacaspase 9, 8A53 chain C](https://www.rcsb.org/structure/8A53): query TM about 0.694, but target TM only 0.175–0.177. It covers 96.8–98.4% of the short design and only 19.7–20.0% of the 305-residue modeled chain. A region match is not a whole-protein duplicate.
- **A41H** has no Foldseek hits in one pose and only weak partial hits in the other; the strongest followed-up PDB query TM is 0.428. **short_018** also has relatively weak followed-up matches. Absence or weakness of retrieved hits does not establish a novel fold, especially with short models and finite search sensitivity.
- **A41H and its parent** remain one lineage. Across the four cross-pose comparisons, their reciprocal TM scores are approximately 0.947–0.971. The five parent lineages are distinguishable in the internal structural comparisons, but this is not a validated count of fold classes.

## Scope, reproducibility and limits

The [frozen plan](plan.json) precedes search results. [summary.json](summary.json) gives compact per-design results, [analysis.json](analysis.json) retains all returned-hit metrics and local top-ten alignments, and [original-PDB verification](independent_pdb_verification.json) documents the named-chain checks. Full Foldseek responses are losslessly compressed alongside tickets and statuses; follow-up coordinate/sequence arrays are archived. Both length normalizations and alignment coverage are retained. All structures are compared in a binder-only fitted frame, separate from the receptor-fitted interface RMSD used earlier. The plotted transformation is checked against TM-align's distance markers.

The live server advertises **PDB100 20240101** and **AlphaFold/Swiss-Prot v6**; the full database parameters are in [database_inventory.json](database_inventory.json). The old PDB structural snapshot limits recent-structure coverage, even though the separate RCSB sequence search is current. Long_092 approaches the server's hit-count limit in AlphaFold/Swiss-Prot. Only the first ten hits per query/database receive local TM-align, and no exhaustive all-PDB or all-UniProt structural alignment was attempted. AlphaFold matches are predicted structures. No universal novelty verdict follows from these results.

Queries use the exact sequences and binder coordinates already available in the public project repository; [the provenance check](public_provenance_check.json) records source URLs and hashes. No new sequence designs, emails, upstream issues/PRs or challenge submission were created by this review.

## Recommendation

Keep **A41H as the leading pH-direction hypothesis**, and keep **long_092 and short_026** as binding candidates. Keep short_015 as a transparently labeled known-backbone backup; do not spend compute merely to make a fold look novel. The larger unresolved problems remain actual human binding, assay-tag compatibility, pH release and mouse cross-reactivity.

If expanding the shortlist, use one bounded apical generation batch to seek additional structurally distinct, well-supported parents, and apply this same sequence/structure audit before adding them. Avoid filling twenty slots with near-identical His variants. The current draft remains a computational submission, not a validated binder panel.

Reproduce interpretation locally with `scripts/analyze_novelty.py`, `scripts/verify_novelty_matches.py`, `scripts/plot_novelty_review.py`, then `scripts/summarize_novelty.py`. Network acquisition commands are in `scripts/novelty_search.py`; submission commands reuse archived tickets and do not silently resubmit jobs. The local Swiss-Prot fallback uses `scripts/local_swissprot.py` and the official URLs/checksums in its frozen plan; large downloaded archives and executables stay in ignored `outputs/`.

Method references: [RCSB search API](https://search.rcsb.org/), [Foldseek](https://github.com/steineggerlab/foldseek), [NCBI BLAST URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html), [TM-align bindings](https://github.com/jvkersch/tmtools). These sources describe the methods; the candidate conclusions above come from this repository's archived calculations.
'''
    (OUT/'RESULTS.md').write_text(text)
    print(status)


if __name__=='__main__':main()
