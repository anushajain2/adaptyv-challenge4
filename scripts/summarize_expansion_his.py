"""Create factual stage reports without changing any screen or selection decision."""
import argparse,json
from pathlib import Path
from prepare_controls import ROOT
from prepare_expansion_his import PARENTS,choose_confirmation


def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['panel','confirmation']);a=p.parse_args();folder='expansion_his_'+a.phase;out=ROOT/'reports'/folder;m=json.loads((ROOT/'data'/folder/'manifest.json').read_text());d=json.loads((out/'results.json').read_text());cases={c['id']:c for c in m['cases']}
    new=[c for c in d['candidate_summaries'] if cases[c['id']].get('parent_id') in PARENTS];passing=[c for c in new if c['all_selected_seeds_full_screen']]
    lines=[f'# New-parent histidine {a.phase}',f"Completed {len(d['rows'])} predictions for {len(m['cases'])} sequences. Selected structural controls {'pass' if d['controls_pass'] else 'FAIL'}. {len(passing)}/{len(new)} new variants pass every selected raw check across {len(d['prediction_seeds'])} seed(s). These are computational structural results, not experimental binding or pH-switching measurements.", '| Sequence | Selected ipTM | Passing samples | Selected raw decision |', '|---|---|---:|---|']
    for c in d['candidate_summaries']:
        rs=c['selected_models'];score=' / '.join(f"{r['metrics']['design_to_target_iptm']:.3f}" for r in rs) or 'No eligible model';fail='; '.join(f"{r['seed']}: {', '.join(r['screen']['failed'])}" for r in rs if not r['screen']['pass'])
        decision='Pass' if c['all_selected_seeds_full_screen'] else fail or 'Missing eligible model or batch control failure'
        lines.append(f"| {c['id']} | {score} | {c['passing_samples']}/{len(d['prediction_seeds'])*d['samples_per_seed']} | {decision} |")
    lines+=['',f"Target-integrity rejects: {sum(not r['target_integrity_pass'] for r in d['rows'])}/{len(d['rows'])}. Passing unselected samples never replace a failed representative. Every sample and failed check is retained in [results.json](results.json)."]
    if a.phase=='panel':
        chosen=choose_confirmation(d,m)
        lines+=['',('Frozen follow-up selects: '+', '.join(chosen)+'.') if chosen else 'No automatic confirmation: no eligible variants or the selected controls fail.','The unchanged parents and A41H comparison sequences are retained in confirmation. At most two new variants per parent advance; ranking uses selected interface ipTM, not prepared geometry or pKa.']
    else:
        summaries={c['id']:c for c in d['candidate_summaries']}
        lines+=['','## Matched parents','']
        for parent in PARENTS:lines.append(f"- {parent}: {'passes both selected raw checks' if summaries[parent]['all_selected_seeds_full_screen'] else 'does not pass both selected raw checks'}. Variant-only success does not demonstrate improved affinity.")
        lines+=['','Only variants passing both selected raw checks, with controls passing, proceed to the separately frozen preparation and PROPKA checks. A41H is an earlier computational comparison, not a measured pH-positive control.']
    prep=out/'preparation/results.json'
    if prep.exists():
        pr=json.loads(prep.read_text());rs=pr['rows'];lines+=['','## Uniform preparation',f"Completed {len(rs)}/{len(pr['plan'])} planned selected-model preparations; {sum(r['status']=='ok' for r in rs)} executions succeeded and {sum(r.get('displacement',{}).get('pass',False) for r in rs)} meet the movement limits. Excluded prepared models remain recorded. Prepared evidence never changes raw calls."]
    lines+=['','Each phase reserves at most $2 within the $4 round allocation and $30 gross project cap. No new generation, mouse inference, upstream issue/PR or challenge submission is performed.']
    (out/'RESULTS.md').write_text('\n\n'.join(lines[:2])+'\n\n'+'\n'.join(lines[2:])+'\n')
    print('Wrote',out/'RESULTS.md')

if __name__=='__main__':main()
