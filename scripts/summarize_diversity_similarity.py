"""Close the two-parent similarity review without claiming novel folds."""
import json
from novelty_search import ROOT, sha, save, now
from analyze_novelty import human_tfr1_retrieved

OUT=ROOT/'reports/diversity_novelty'


def main():
    plan=json.loads((OUT/'plan.json').read_text());data=json.loads((OUT/'analysis.json').read_text())
    expected={(q['id'],d['path']) for q in plan['structures'] for d in plan['foldseek']['databases']}
    assert {(r['query'],r['database']) for r in data['foldseek'] if r['status']=='complete'}==expected
    assert data['foldseek_positive_retrieval_pass']
    rc={r['query']:r for r in data['rcsb_sequence']};bs={r['query']:r for r in data['swissprot_sequence']}
    assert set(rc)=={q['id'] for q in plan['sequences']} and all(r['status']=='complete' for r in rc.values())
    assert set(bs)==set(plan['blast']['query_ids']) and human_tfr1_retrieved(bs['human_tfr1_crop_control']['hits'])
    assert any(h['target']=='6WRW_2' and h['sequence_identity']==1 for h in rc['published_2ds25_5_control']['hits'])
    assert any(h['sequence_identity']==1 and h['query_coverage']==1 for h in rc['human_tfr1_crop_control']['hits'])
    for q in plan['structures']:
        assert sha(ROOT/q['source'])==q['source_sha256'] and sha(ROOT/q['query_file'])==q['query_sha256']
    rows=[]
    for q in [q for q in plan['sequences'] if q['kind']=='candidate']:
        structural=[]
        for r in data['foldseek']:
            if r.get('candidate')!=q['id']:continue
            hits=[h for h in r['hits'] if 'local_TM_align' in h]
            best=max(hits,key=lambda h:h['local_TM_align']['query_normalized_TM']) if hits else None
            reciprocal=max((min(h['local_TM_align']['query_normalized_TM'],h['local_TM_align']['target_normalized_TM']) for h in hits),default=0)
            structural.append({'query':r['query'],'database':r['database'],'reported_hits':r['hit_count'],
                               'best_query_TM_hit_among_top10':best,'max_reciprocal_TM_among_top10':reciprocal})
        rows.append({'id':q['id'],'sequence':q['sequence'],'rcsb_sequence_hits':rc[q['id']]['total_count'],
            'swissprot_HSPs':len(bs[q['id']]['hits']),
            'swissprot_min_evalue':min((h['evalue'] for h in bs[q['id']]['hits']),default=None),
            'swissprot_close_flags':sum(h['close_sequence_review_flag'] for h in bs[q['id']]['hits']),
            'swissprot_masking':bs[q['id']]['masking'],'structural':structural})
    summary={'status':'COMPLETE_WITH_DATABASE_LIMITS','completed_utc':now(),'candidates':rows,
             'reviewed_candidate_ids':[r['id'] for r in rows],
             'controls':{'Foldseek_6WRW_retrieval':True,'RCSB_published_binder_retrieval':True,
                         'RCSB_human_crop_retrieval':True,'SwissProt_P02786_retrieval':True},
             'new_modal_compute_usd':0,'new_NCBI_RID_submitted':False,
             'interpretation':'Weak/partial retrieved structural matches and no close-sequence flags at the recorded settings; finite top-ten follow-up and incomplete databases do not prove novel folds or lack of homology.'}
    save(OUT/'summary.json',summary)
    lines=['# Similarity review of two new computational parents','',
        '**Review complete within the stated database limits.** Both retained parents have zero RCSB sequence hits and zero Swiss-Prot close-sequence flags under the frozen settings. Foldseek returned structural matches; the followed-up matches are weak or partial, and do not establish whole-protein duplication. This is not a novel-fold verdict.','',
        'All five binder-only structure queries completed against both databases, and every retrieval control passed. The published 2DS25.5 control recovers 6WRW; the local Swiss-Prot search recovers human P02786 exactly. No new NCBI job or Modal compute was used.','',
        '## Sequence results','',
        '| Design | RCSB hits | Swiss-Prot HSPs | Smallest E-value | Close flags | Masked intervals |',
        '|---|---:|---:|---:|---:|---|']
    for r in rows:
        lines.append(f"| {r['id']} | {r['rcsb_sequence_hits']} | {r['swissprot_HSPs']} | {r['swissprot_min_evalue']:.3g} | {r['swissprot_close_flags']} | `{json.dumps(r['swissprot_masking'])}` |")
    lines+=['','The close-sequence review heuristic is identity ≥50%, query coverage ≥80% and E ≤0.001 per HSP. It is not a competition cutoff. BLAST+ 2.17.0 uses BLOSUM62, word size 3, SEG soft masking, composition adjustment 2 and up to 100 targets against the verified UniProtKB/Swiss-Prot canonical 2026_03 collection (575,748 sequences). Isoforms from the separate supplement are excluded. Masking and short/compositionally biased queries limit sensitivity.','',
        '## Structural results','',
        'The table reports the largest query-normalized TM score among the first ten server-ranked hits, separately for each pose/database. Target-normalized TM and coverage prevent a partial match to a large protein from appearing to be a whole-protein duplicate.','',
        '| Query | Database | Target | Query TM | Target TM | Query coverage | Target coverage |',
        '|---|---|---|---:|---:|---:|---:|']
    for r in rows:
        for s in r['structural']:
            h=s['best_query_TM_hit_among_top10']
            if not h:continue
            a=h['local_TM_align'];target=h['target'].split(' ')[0]
            lines.append(f"| {s['query']} | {s['database']} | {target} | {a['query_normalized_TM']:.3f} | {a['target_normalized_TM']:.3f} | {a['query_coverage']:.1%} | {a['target_coverage']:.1%} |")
    lines+=['',
        'For short_055, the largest PDB query score is about 0.607 against part of 4HXG chain F, but target-normalized TM is only about 0.085 and target coverage about 9.3%. This is a small-region match. Long_037 has only modest followed-up scores; its largest PDB query score is about 0.547, with target TM about 0.145. No followed-up candidate match has both length-normalized TM scores above 0.5. These observations are limited to the retrieved top-ten hits, not exhaustive alignments over all known proteins.','',
        'The server advertises PDB100 **20240101** and AlphaFold/Swiss-Prot **v6**. The old structural snapshot, finite search sensitivity, limited follow-up and predicted AlphaFold structures prevent broad novelty claims. A separate comparison of the confirmed poses against the current draft is archived in [confirmed_diversity.json](../diversity_confirmation/confirmed_diversity.json); it tests redundancy within this panel, not global fold novelty.','',
        '## Provenance and decision','',
        'The [frozen plan](plan.json), [full analysis](analysis.json), raw/gzipped server responses, local BLAST JSON and hashes preserve the review. [Anonymous public-byte verification](public_provenance_check.json) preceded API upload. Queries contain binder atoms only. Both selected poses are retained. No sequences were changed to improve these results.','',
        'The review identifies no close-sequence flag requiring exclusion of either new parent. Retain long_037 and short_055 as additional computational binding hypotheses, with their preparation checks and weaker/stronger confidence differences explicit. Neither contains histidine; pH release, mouse binding, tag compatibility and experimental affinity remain unestablished.','',
        'Method references: [Foldseek](https://github.com/steineggerlab/foldseek), [RCSB search API](https://search.rcsb.org/), [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), [UniProt release files](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).']
    (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    save(OUT/'manifest.json',{'status':summary['status'],'files':{str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'},
        'scripts':{f'scripts/{s}':sha(ROOT/'scripts'/s) for s in ['diversity_similarity.py','summarize_diversity_similarity.py','analyze_novelty.py','novelty_search.py']}})
    print(summary['status'],summary['reviewed_candidate_ids'])


if __name__=='__main__':main()
