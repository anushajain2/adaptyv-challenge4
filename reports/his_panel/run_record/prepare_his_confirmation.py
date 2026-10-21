"""Freeze at most six screen-passing variants plus both parents and controls."""
import hashlib
import json
from prepare_controls import ROOT
from prepare_next_round import finish


def choose(result, manifest):
    if not result['controls_pass']:return None
    lookup={c['id']:c for c in manifest['cases']}
    chosen=[]
    for parent in ('apical_long_015','apical_short_008'):
        eligible=[r for r in result['candidate_summaries'] if lookup[r['id']].get('parent_id')==parent
                  and r['all_selected_seeds_candidate_screen'] and r['id']!='apical_long_015_D51H_T70H']
        eligible.sort(key=lambda r:(-min(x['metrics']['design_to_target_iptm'] for x in r['selected_models']),r['id']))
        chosen.extend(r['id'] for r in eligible[:3])
    return chosen


def main():
    src=ROOT/'data/his_panel';out=ROOT/'data/his_confirmation'
    source=ROOT/'reports/his_panel/results.json';r=json.loads(source.read_text());m=json.loads((src/'manifest.json').read_text())
    chosen=choose(r,m)
    if chosen is None:raise SystemExit('Controls fail: no automatic confirmation launch.')
    out.mkdir(exist_ok=False)
    ids=chosen+['apical_long_015','apical_short_008']
    m['cases']=[c for c in m['cases'] if c['id'] in ids or c['arm']=='control']
    for c in m['cases']:
        for ext in ('.cif','.npz','.pdb'):
            name=c['id']+ext;assert hashlib.sha256((src/name).read_bytes()).hexdigest()==m['input_sha256'][name]
            (out/name).write_bytes((src/name).read_bytes())
    m.update(selected_ids=ids,source_run=r['run'],source_results_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
             scope='Fixed fresh-seed confirmation of shortlisted His variants and both unmodified parents',
             refolding={**m['refolding'],'seeds':[1103,1201]},
             selection_details={'shortlisted_variant_ids':chosen,'max_per_parent':3,'ranking':'Selected raw-screen-pass interface ipTM descending; ties lexical id; excluded failed double-mutant comparator'})
    finish(out,m)
    print(json.dumps({'variants':chosen,'cases':len(m['cases']),'predictions':len(m['cases'])*10},indent=2))
if __name__=='__main__':main()
