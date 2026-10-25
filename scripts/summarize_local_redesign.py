"""Write a concise evidence table without overwriting the frozen submission."""
import json
from novelty_search import ROOT,save


def main():
    folder=ROOT/'reports/local_redesign_panel';raw=json.loads((folder/'results.json').read_text());direction=json.loads((folder/'direction.json').read_text())
    m=json.loads((ROOT/'data/local_redesign_panel/manifest.json').read_text());cases={c['id']:c for c in m['cases']}
    summaries={r['id']:r for r in raw['candidate_summaries']};variants=direction['variants']
    text=['# Local redesign comparison results','',
        f"Completed {len(raw['rows'])} predictions across {len(cases)} sequences, with five samples for each of two fresh seeds. "
        f"The raw control panel passed: **{direction['raw_controls_pass']}**. "
        f"Strict prepared control panel passed: **{direction['strict_prepared_controls_pass']}**. "
        f"The separate interface-movement sensitivity panel passed: **{direction['sensitivity_prepared_controls_pass']}**.",'',
        f"**{sum(v['fully_qualified_for_followup'] for v in variants)} of {len(variants)} new sequences qualifies for the preregistered follow-up.** "
        'The frozen eight-design release remains unchanged; nothing has been submitted.','',
        'The [interpretation note](INTERPRETATION.md) explains why A41H with A34K/A38R remains a promising reserve despite missing the full preset rule. It separates that judgment from the unchanged formal selection outcome.','',
        'The strict and sensitivity checks remain separate. Diagnostic pKa values below are shown even when a raw or preparation gate fails; those values cannot qualify a design. '
        'Two prediction seeds measure repeatability in the same predictor, not independent model or experimental validation.','',
        '| Sequence | Raw samples passing /10 | Both selected raw pass | Strict eligible poses /2 | pKa shifts (seed order) | His-only ratios (seed order) |',
        '|---|---:|---|---:|---|---|']
    for family in m['prospective_plan']['families']:
        ids=[family['id']]+[v['id'] for v in variants if v['baseline_id']==family['id']]
        for name in ids:
            rows=sorted([r for r in direction['rows'] if r['id']==name],key=lambda r:r['seed'])
            def values(key):
                result=[]
                for r in rows:
                    d=r['direction']
                    if d['status']!='diagnostic':result.append('excluded');continue
                    value=d[key];value=value[0] if isinstance(value,list) else value
                    result.append(f'{value:.3f}')
                return ', '.join(result)
            s=summaries[name]
            text.append(f"| {name} | {s['passing_samples']} | {s['all_selected_seeds_full_screen']} | {sum(r['strict_eligible'] for r in rows)} | {values('shifts')} | {values('his_only_Kd6_over_Kd7_4')} |")
    text += ['', '## Matched original parents','']
    for family in m['prospective_plan']['families']:
        name=family['unmodified_parent'];s=summaries[name];rows=[r for r in direction['rows'] if r['id']==name]
        text.append(f"- `{name}`: {s['passing_samples']}/10 raw samples pass; both selected poses pass: {s['all_selected_seeds_full_screen']}; {sum(r['strict_eligible'] for r in rows)}/2 strict prepared poses eligible.")
    text += ['', '## Exclusions and interpretation','']
    for r in direction['rows']:
        if r['kind']=='novel_design' or r['strict_eligible']:continue
        text.append(f"- Control `{r['id']}`, seed {r['seed']}: raw pass {r['raw_pass']}; strict movement pass {r.get('displacement',{}).get('pass')}; internal geometry pass {r.get('internal_geometry',{}).get('pass')}; see the complete preparation/protonation record in `direction.json`.")
    for v in variants:
        text.append(f"- `{v['id']}`: changes {', '.join(v['mutations'])}; strict structure-and-direction in both poses {v['strict_structure_and_direction_both']}; "
                    f"matched His baseline eligible {v['baseline_strict_eligible_both']}; original parent eligible {v['unmodified_parent_strict_eligible_both']}; baseline-comparison precondition {v['prospective_improvement_heuristic']} (this alone does not establish improvement).")
    text += ['', 'ProteinMPNN compatibility scores are not protonation-dependent binding energies. The direction calculations use the same bound conformation with and without the receptor; they omit an unbound ensemble and full proton linkage. '
             'Histidine neighbour distances are descriptive, not proof of hydrogen bonds or electrostatic repulsion. No experimental release, affinity, mouse cross-reactivity or tag compatibility is established.','',
             'The [frozen plan](../local_redesign/plan.json), [sampling records](../local_redesign/mpnn/results.json), [raw results](results.json), '
             '[direction and neighbour audit](direction.json), [PROPKA term decomposition](mechanism_audit.json), [evidence audit](evidence_audit.json), and [selected models plus all confidence samples](review_bundle/manifest.json) retain failures and provenance.','']
    budget=folder/'run_record/budget_reconciliation.json'
    if budget.exists():
        b=json.loads(budget.read_text());text += [f"Observed gross cost increase: **${b['observed_panel_increment_usd']:.3f}**. Conservative cumulative gross accounting: **${b['conservative_current_total_usd']:.3f}/$30**. The round cap is $6; the prediction panel reservation was $2. No credit subtraction is used.",'']
    if not direction['followup_ids']:text += ['The conditional tag follow-up was skipped because no redesign met the frozen improvement rule.','']
    (folder/'RESULTS.md').write_text('\n'.join(text))
    save(folder/'completion.json',{'predictions':len(raw['rows']),'new_sequences':len(variants),'raw_controls_pass':raw['controls_pass'],
        'new_both_raw_pass':sum(summaries[v['id']]['all_selected_seeds_full_screen'] for v in variants),
        'strict_structure_and_direction_both':sum(v['strict_structure_and_direction_both'] for v in variants),
        'qualified_improvements':sum(v['fully_qualified_for_followup'] for v in variants),'followup_ids':direction['followup_ids'],
        'challenge_submitted':False,'frozen_submission_changed':False})


if __name__=='__main__':main()
