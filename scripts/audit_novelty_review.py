"""Verify completed search evidence and write a content-hash manifest."""
import json
import subprocess
from pathlib import Path

from novelty_search import ROOT, OUT, sha, save
from analyze_novelty import load_fold, flatten_hits
from prepare_controls import read_structure, sequence


def main():
    plan=json.loads((OUT/'plan.json').read_text());summary=json.loads((OUT/'summary.json').read_text())
    assert sha(OUT/'draft_audit_before_review.json')==plan['draft_audit_sha256']
    assert len(plan['structures'])==13 and len(plan['sequences'])==8
    for c in plan['structures']:
        assert sha(ROOT/c['query_file'])==c['query_sha256']
        assert sha(ROOT/c['source'])==c['source_sha256']
        s=read_structure(ROOT/c['query_file']);assert [ch.id for ch in s[0]]==['B']
        assert sequence(s[0]['B'])==c['sequence']
        result=load_fold(OUT/'foldseek'/c['id']/'result.json')
        assert result['mode']==plan['foldseek']['mode']
        assert result['queries'][0]['sequence']==c['sequence']
        assert {r['db']for r in result['results']}=={'pdb100','afdb-swissprot'}
    assert summary['all_13_structural_queries_complete']
    assert all(summary['controls'][k] for k in ['foldseek_6WRW_recovered','rcsb_6WRW_recovered','rcsb_human_crop_recovered'])
    if summary['swissprot_all_seven_queries_complete']:
        assert summary['controls']['swissprot_human_crop_recovered']
        if summary.get('swissprot_sequence_source')=='local_UniProtKB_SwissProt_canonical':
            lp=json.loads((OUT/'local_swissprot/plan.json').read_text())
            lr=json.loads((OUT/'local_swissprot/run.json').read_text())
            assert lr['state']=='COMPLETE' and lr['plan_sha256']==sha(OUT/'local_swissprot/plan.json')
            assert lp['original_plan_sha256']==sha(OUT/'plan.json')
            assert lr['result_sha256']==sha(OUT/'local_swissprot/results.json')
            assert lp['query_sha256']==sha(OUT/'local_swissprot/queries.fasta')
            assert (OUT/'local_swissprot/queries.fasta').read_text()==json.loads((OUT/'blast/request.json').read_text())['QUERY']
    # Preserve this historical review's payload when later rounds extend the draft.
    commit=json.loads((OUT/'public_provenance_check.json').read_text())['commit']
    for name in ['submission_draft.csv','submission_draft.fasta']:
        rel='reports/submission_draft/'+name
        prior=subprocess.check_output(['git','show',f'{commit}:{rel}'],cwd=ROOT)
        archived=OUT/('reviewed_'+name)
        assert prior==(archived.read_bytes() if archived.exists() else (ROOT/rel).read_bytes())
    analysis=json.loads((OUT/'analysis.json').read_text())
    for r in analysis['foldseek']:
        for h in r['hits']:
            if 'followup_target_file' in h:assert sha(ROOT/h['followup_target_file'])==h['followup_target_sha256']
    scripts=['novelty_search.py','local_swissprot.py','analyze_novelty.py','verify_novelty_matches.py','plot_novelty_review.py','summarize_novelty.py','audit_novelty_review.py']
    manifest={'status':summary['status'],'checks':{'query_sequences_and_binder_only_scope':True,
              'all_structure_queries_and_databases_present':True,'PDB_retrieval_controls':True,
              'SwissProt_complete_and_control_pass':summary['controls']['swissprot_human_crop_recovered'],
              'reviewed_submission_CSV_and_FASTA_unchanged':True,'archived_pre_review_audit_hash':True,
              'followup_coordinate_hashes':True},
              'files':{str(p.relative_to(ROOT)):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='manifest.json'},
              'scripts':{f'scripts/{name}':sha(ROOT/'scripts'/name) for name in scripts}}
    save(OUT/'manifest.json',manifest)
    print('Evidence verified:',manifest['status'],';',len(manifest['files']),'files hashed; submission sequences unchanged.')


if __name__=='__main__':main()
