"""Freeze one bounded apical batch and its prospective advancement rule."""
import copy
import json
from Bio import SeqIO
from novelty_search import ROOT, sha, save, now
from receptor_preparation_v2 import PROTOCOL


def main():
    out=ROOT/'data/diversity_generation';out.mkdir(exist_ok=False)
    src=ROOT/'data/apical_expansion_fixed'
    m=json.loads((src/'manifest.json').read_text())
    arms=[]
    for i,a in enumerate(m['arms']):
        b=copy.deepcopy(a);b.update(id='diversity_'+('short' if i==0 else 'long'),num_designs=60 if i==0 else 40,
                                  generation_seed=20261401+i,inverse_fold_seed=20261501+i)
        (out/f"{b['id']}.yaml").write_bytes((src/f"{a['id']}.yaml").read_bytes());arms.append(b)
    for name in m['input_sha256']:
        if not name.endswith('.yaml'):(out/name).write_bytes((src/name).read_bytes())
    prior=[];sources={}
    for path in ['reports/pilot/pilot_sequences.fasta','reports/generation_comparison/all_designs.fasta',
                 'reports/broad_campaign/all_designs.fasta','reports/apical_expansion_fixed/all_designs.fasta']:
        sources[path]=sha(ROOT/path)
        prior.extend({'id':r.id,'sequence':str(r.seq),'source':path} for r in SeqIO.parse(ROOT/path,'fasta'))
    assert len(prior)==495,len(prior)
    draft=json.loads((ROOT/'reports/submission_draft/audit.json').read_text())
    catalog={'prior_designs':prior,'prior_sources_sha256':sources,'current_draft':draft['candidates'],
             'draft_audit_sha256':sha(ROOT/'reports/submission_draft/audit.json')}
    save(out/'prior_catalog.json',catalog)
    m.pop('repair',None)
    m.update(prepared_utc=now(),arms=arms,scope='100 additional apical miniproteins: 60 short and 40 long; seek distinct parents',
             refolding={**m['refolding'],'seeds':[2203]},max_prediction_seconds=3000,experiment_cap_usd=2.5,
             budget_allocation={'generation_and_screen_usd':2.5,'conditional_confirmation_usd':1.5,'round_total_usd':4},
             structural_preparation={**PROTOCOL,'cohort':'Every confidence-selected confirmation candidate and control, including failures, in this diversity round.'},
             diversity_policy={'exact_duplicates':'Exclude exact matches to all 495 prior generated sequences.',
                'sequence':'Exclude >=80% global identity to any current draft sequence or already chosen new candidate.',
                'structure':'Among full raw-screen passes, exclude min(query_TM,target_TM)>=0.80 against any current draft selected bare pose or already chosen new candidate. Binder-only TM-align; descriptive redundancy heuristic, not a validated fold-class test.',
                'selection':'Round-robin short/long arms in frozen within-class rank order; at most six. All comparisons archived; no substitutions after confirmation.',
                'confirmation':'Only if all three screening controls pass. Seeds 2309 and 2411, five samples each, all three controls, unchanged target-integrity-aware selection and raw structural gates. Both selected seeds must pass.',
                'admission':'A two-seed pass is a provisional parent. Apply uniform preparation, sequence-liability and database-similarity review before changing the submission; no pH or mouse-binding claim. No automatic extra generation.'},
             prior_campaign='reports/broad_campaign/results.json')
    m['screen']['decision']='One 100-design batch; at most six distinct full-screen passes enter frozen confirmation. No automatic expansion.'
    m['input_sha256']={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    save(out/'manifest.json',m)
    report=ROOT/'reports/diversity_generation';report.mkdir(exist_ok=False)
    save(report/'plan.json',m)
    print('Frozen 100-design round, 495 prior sequences, 6 current draft sequences, $4 total reservation.')


if __name__=='__main__':main()
