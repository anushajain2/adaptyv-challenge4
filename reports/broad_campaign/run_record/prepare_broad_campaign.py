"""Prepare the authorized 200-design screen with arm-specific epitopes."""
import copy
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]

def prepare():
    source = ROOT/'data/generation_comparison'
    out = ROOT/'data/broad_campaign'
    out.mkdir(exist_ok=False)
    old = json.loads((source/'manifest.json').read_text())
    rows = list(csv.DictReader((source/'spatial_map.csv').open()))
    lookup = {int(r['native_position']):r for r in rows}
    for n in (208,210,211,327,331):
        assert lookup[n]['bound_spatial_pass'] == 'True'
    for n in (327,331): assert lookup[n]['conserved'] == 'True'
    for name in ('target_6wrw.cif','spatial_map.csv'):
        (out/name).write_bytes((source/name).read_bytes())
    for c in old['controls']:
        for ext in ('.pdb','.cif','.npz'):
            name=c['id']+ext
            (out/name).write_bytes((source/name).read_bytes())
    (out/'7eow.cif').write_bytes((ROOT/'configs/7eow.cif').read_bytes())
    base = yaml.safe_load((source/'unconstrained.yaml').read_text())
    nano = yaml.safe_load((ROOT/'configs/apical_tip_basic_direction_probe_nanobody-anything.yaml').read_text())['entities'][1]
    for item, value in zip(nano['file']['design_insertions'], ('1..3','1..3','3..9')):
        item['insertion']['num_residues'] = value
    arms=[]
    for i,(name,epitope,length,protocol) in enumerate([
        ('apical_short','apical','55..72','protein-anything'),
        ('conserved_short','conserved','55..72','protein-anything'),
        ('apical_long','apical','80..110','protein-anything'),
        ('conserved_long','conserved','80..110','protein-anything'),
        ('apical_nanobody','apical','120..130','nanobody-anything'),
    ]):
        spec=copy.deepcopy(base)
        hot=[208,210,211] if epitope=='apical' else [327,331]
        spec['entities'][0]['file']['binding_types'][1]['chain']['binding']=','.join(str(n-189) for n in hot)
        if protocol=='nanobody-anything': spec['entities'][1]=copy.deepcopy(nano)
        else: spec['entities'][1]={'protein':{'id':'B','sequence':length}}
        (out/f'{name}.yaml').write_text(yaml.safe_dump(spec,sort_keys=False))
        arms.append({'id':name,'epitope':epitope,'protocol':protocol,'num_designs':40,
                     'binder_length_range':[int(n) for n in length.split('..')],
                     'generation_seed':20261023+i,'inverse_fold_seed':20261123+i,
                     'hotspots_native':hot,'hotspots_crop':[n-189 for n in hot],
                     'reference_contact_region_crop':old['reference_contact_region_crop'] if epitope=='apical' else list(range(137,159)),
                     'region_definition':'published apical footprint' if epitope=='apical' else 'native 326-347 conserved segment'})
    manifest={k:old[k] for k in ('boltzgen_commit','target_sequence','target_native_map','receptor','controls','spatial_context')}
    manifest.update(prepared_utc=datetime.now(timezone.utc).isoformat(),scope='Authorized 200-design exploratory screen; no experimental binding claims',arms=arms,
        refolding={**old['refolding'],'seeds':[307],'diffusion_samples_per_seed':1},
        max_prediction_seconds=6600,experiment_cap_usd=6,total_cap_usd=30,
        budget_allocation={'generation_and_triage_usd':16,'deep_evaluation_usd':7,'contingency_usd':3},
        screen={**old['screen'],'replicates':'ONE seed/sample; screening only, not a qualified parent',
            'decision':'Cost/quality checkpoint after 200. No automatic expansion. Shortlist at most 30 diverse leads for deep evaluation.',
            'ranking':'Within binder class: geometry eligibility first, then equal-weight percentile ranks of interface ipTM descending, contact-pair PAE ascending, hotspot fraction descending, binder fold RMSD ascending. Sequence diversity selection at 80% identity; preserve epitope representation.',
            'nanobody_confidence':'Record ipTM>=0.5 as descriptive only; no validated positive VHH control on this target.',
            'min_hotspots_contacted':2,
            'control_policy':'One-sample matched 2DS25.5 plus presumed negatives; fail blocks advancement, does not trigger automatic reruns.'})
    manifest['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({'designs':sum(a['num_designs'] for a in arms),'arms':arms},indent=2))
if __name__=='__main__': prepare()
