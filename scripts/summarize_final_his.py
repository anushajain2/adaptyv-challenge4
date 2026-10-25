"""Write a result-led final-panel report without promoting failed hypotheses."""
import json
from novelty_search import ROOT,sha,save

def main():
    out=ROOT/'reports/final_his_panel';raw=json.loads((out/'results.json').read_text());direction=json.loads((ROOT/'reports/final_his_direction/results.json').read_text());prep=json.loads((out/'preparation/results.json').read_text())
    summaries={c['id']:c for c in raw['candidate_summaries']}
    v=direction['variants'];qualified=[c for c in v if c['fully_qualified_direction_pass_both_seeds']]
    lines=['# Final anchor-preserving His panel — 25 October 2026','',
        f"**{len(v)} variants; {len(raw['rows'])} mutant/parent/control predictions across two seeds; {sum(c['raw_pass_both_seeds'] for c in v)} variants pass both selected raw models; {len(qualified)} fully qualified two-pose direction hypotheses.** These are computational checks, not biological validation.",'',
        '## Frozen experiment','',
        'The contact audit used the two independently predicted long_037 poses from seeds 2309 and 2411. All charged residues, glycine/proline/cysteine, high-contact anchor proxies and five residues at each binder terminus were preserved. Sites had to be peripheral in both poses (free SASA ≥25 Å², buried SASA 2–65 Å², target distance ≤5 Å). Ten neutral single substitutions and two pairs were fixed before prediction. The pairs include S37, which has a basic target neighbour in both source poses; this geometric rationale does not establish favorable electrostatics or a binding-energy contribution.','',
        'All twelve variants, the unchanged parent and three controls received five samples at each fresh seed 2503 and 2609. Entire binder templates were excluded and the target-only template/empty-MSA setup retained. The existing target-integrity-aware confidence selection, structural checks and direction thresholds were unchanged. Favorable unselected samples cannot rescue a failed selected pose.','',
        '## Results','',
        '| Variant | Raw passing samples / 10 | Both selected raw pass | Both-pose direction heuristic | Fully qualified |',
        '|---|---:|---|---|---|']
    for c in v:
        lines.append(f"| {c['id']} | {summaries[c['id']]['passing_samples']} | {c['raw_pass_both_seeds']} | {c['direction_pass_both_seeds']} | {c['fully_qualified_direction_pass_both_seeds']} |")
    parent=summaries['diversity_long_037']
    lines+=['',f"The unchanged parent passes {parent['passing_samples']}/10 raw samples; both selected poses pass: **{parent['all_selected_seeds_full_screen']}**. All six selected raw controls pass: **{raw['controls_pass']}**. All six prepared controls pass: **{direction['prepared_controls_pass']}**.",'',
        f"All {len(prep['rows'])} confidence-selected models, including failures and controls, received one identical v2 preparation. Execution succeeded for {sum(p['status']=='ok' for p in prep['rows'])}; movement checks pass for {sum(p.get('displacement',{}).get('pass',False) for p in prep['rows'])}. Warning-free PROPKA completed for {sum(v['status']=='ok' for p in prep['rows'] for v in p.get('propka',{}).values())}/64 free/bound forms; warning failures remain excluded. Raw failures cannot be rescued by preparation. Full qualification requires both prepared parent poses and every prepared control to pass, in addition to the mutant-specific criteria.",'',
        '### Per-pose direction evidence','',
        '| Variant | Seed | Status | Bound−free pKa shifts | His-only ratio proxy | Exclusion / failed checks |',
        '|---|---:|---|---|---:|---|']
    for c in v:
        for r in c['models']:
            d=r['direction'];shifts=', '.join(f'{x:.3f}' for x in d.get('shifts',[])) or '—';ratio=f"{d['his_only_Kd6_over_Kd7_4']:.3f}" if 'his_only_Kd6_over_Kd7_4' in d else '—'
            reason=r.get('exclusion',', '.join(k for k,x in d.get('checks',{}).items() if not x)) or 'none'
            lines.append(f"| {c['id']} | {r['seed']} | {d['status']} | {shifts} | {ratio} | {reason} |")
    failed_controls=[r for r in direction['rows'] if r['kind']!='novel_design' and not r['prepared_eligible']]
    lines+=['','### Preparation-control limitations','']
    for r in failed_controls:
        lines.append(f"- {r['id']}, seed {r['seed']}: {r.get('exclusion','ineligible')}. Backbone displacement: {json.dumps(r.get('displacement',{}),sort_keys=True)}")
    if not failed_controls:lines.append('All selected prepared controls pass the frozen checks.')
    lines+=['','## Decision','']
    if not direction['tag_selected_ids']:
        lines+=['**No variant meets the prospectively frozen admission rule. The conditional tag job is not launched, and none of these twelve variants enters the submission CSV.** The current eight-sequence shortlist remains the baseline; the parent’s new confirmation evidence must be included in its ranking rationale. This is the planned stopping point for computational experimentation, not a reason to relax thresholds or start another generation batch.']
    else:
        lines+=['The following variants qualify for the separately bounded matched bare/tagged follow-up: '+', '.join(direction['tag_selected_ids'])+'. Tag results are required before any submission admission.']
    lines+=['','## Evidence and limits','',
        'The [prospective plan](plan.json), [contact audit](../../data/final_his_panel/contact_audit.json), [all raw sample scores](results.json), [prepared models](preparation/results.json), [direction exclusions](../final_his_direction/results.json) and [review bundle](review_bundle/manifest.json) preserve failures as well as passes. The bundle includes all selected raw/prepared structures and complete five-sample confidence arrays for every case/seed.','',
        'PROPKA compares free and bound forms in a fixed bound-like conformation. The introduced-His-only ratio omits other titratable sites, conformational populations and full proton linkage. It is not a predicted experimental KD ratio. Preparation is finite restrained minimization, not proof of equilibrium or global stereochemical validity. No new mouse-binding or assay-tag claim follows from this panel.','',
        'The entire final round is limited to $3 gross including credits: $1.50 for this panel and $1.50 reserved for a conditional tag job. The local preparation and analysis add no Modal compute. No challenge submission or upstream issue/PR is performed.']
    cost=json.loads((out/'run_record/budget_reconciliation.json').read_text())
    lines+=['',f"Recorded gross increment: **${cost['observed_panel_increment_usd']:.2f}**. Latest account meter: ${cost['after_gross_usd']:.2f}; conservative cumulative accounting, retaining the prior higher snapshot: **${cost['conservative_current_total_usd']:.2f}/$30**. The meter may lag. All cloud jobs are stopped; the conditional tag reservation is unused."]
    (out/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print('Wrote final panel report; qualifying tag candidates:',direction['tag_selected_ids'])
if __name__=='__main__':main()
