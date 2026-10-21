"""Feature-checked generation and one-sample screening inputs for varied binders."""
import argparse
import json
import subprocess
from pathlib import Path
from generation_comparison_worker import config_path, run_stage


def design_count(batches, samples):
    assert batches>0 and samples>0
    return batches*samples


def preflight(out):
    import hydra
    import numpy as np
    from omegaconf import OmegaConf
    from boltzgen.data import const
    source=out/'inputs'
    manifest=json.loads((source/'manifest.json').read_text())
    checks=[]
    for arm in manifest['arms']:
        name=arm['id']; count=str(arm['num_designs'])
        subprocess.run(['boltzgen','configure',str(source/f'{name}.yaml'),'--output',str(out/name),
            '--protocol',arm['protocol'],'--num_designs',count,'--budget',count,'--devices','1',
            '--num_workers','0','--steps','design','inverse_folding'],check=True)
        cfg=OmegaConf.load(config_path(out/name,'design'))
        data=hydra.utils.instantiate(cfg.data)
        feat=data.predict_set.get_sample(source/f'{name}.yaml')
        assert not feat.get('exception',False), name
        mask=feat['design_mask'].numpy().astype(bool)
        assert not mask[:194].any() and mask[194:].any()
        low,high=arm['binder_length_range']
        assert low <= len(mask)-194 <= high, (name,len(mask)-194)
        if arm['protocol']=='protein-anything': assert mask[194:].all()
        else: assert not mask[194:].all(), 'Expected fixed VHH framework'
        bind=feat['binding_type'].numpy()
        hot=(np.flatnonzero(bind[:194]==const.binding_type_ids['BINDING'])+1).tolist()
        assert hot==arm['hotspots_crop'], (name,hot)
        assert not (bind[194:]==const.binding_type_ids['BINDING']).any()
        assert design_count(len(data.predict_set), int(cfg.diffusion_samples))==arm['num_designs']
        checks.append({'arm':name,'target_tokens':194,'binder_tokens':len(mask)-194,
                      'designed_tokens':int(mask.sum()),'hotspots_crop':hot,'diffusion_batches':len(data.predict_set),'diffusion_samples':int(cfg.diffusion_samples),'design_count':design_count(len(data.predict_set),int(cfg.diffusion_samples))})
    (out/'generation_preflight.json').write_text(json.dumps(checks,indent=2))


def assemble(out):
    import numpy as np
    from prepare_controls import make_complex,read_structure,residues,sequence,save_input_cif,save_pdb,fit_mapped
    source=out/'inputs'; manifest=json.loads((source/'manifest.json').read_text())
    target=read_structure(source/'target_6wrw.cif')
    refs=out/'references'; refs.mkdir()
    cases=[]; audit=[]
    for arm in manifest['arms']:
        folder=out/arm['id']/'intermediate_designs_inverse_folded'
        files=sorted(p for p in folder.glob('*.cif') if not p.stem.endswith('_native'))
        assert len(files)==arm['num_designs'], (arm['id'],len(files))
        generated=sorted((out/arm['id']/'intermediate_designs').glob('*.npz'))
        assert len(generated)==arm['num_designs']
        for path in generated:
            with np.load(path,allow_pickle=False) as d:
                hot=(np.flatnonzero(d['binding_type'][:194]==1)+1).tolist()
                assert hot==arm['hotspots_crop'], (path,hot)
                assert not d['design_mask'][:194].any() and d['design_mask'][194:].any()
                assert arm['binder_length_range'][0]<=len(d['design_mask'])-194<=arm['binder_length_range'][1]
            audit.append({'source':str(path.relative_to(out)),'hotspots_crop':hot})
        for index,path in enumerate(files):
            raw=read_structure(path); seqs={c.id:sequence(c) for c in raw[0]}
            targets=[c for c,s in seqs.items() if s==manifest['target_sequence']]
            assert len(targets)==1 and len(seqs)==2
            b=next(c for c in seqs if c!=targets[0])
            normalized=make_complex(residues(raw[0][targets[0]]),residues(raw[0][b]))
            fitted,rms=fit_mapped(normalized,{i:i for i in range(1,195)},target)
            ref=make_complex(residues(target[0]['A']),residues(fitted[0]['B']))
            name=f"{arm['id']}_{index:03d}"; seq=sequence(ref[0]['B'])
            assert arm['binder_length_range'][0]<=len(seq)<=arm['binder_length_range'][1]
            assert set(seq)<=set('ACDEFGHIKLMNPQRSTVWY')
            save_pdb(ref,refs/f'{name}.pdb')
            save_input_cif(ref,refs/f'{name}.cif',manifest['target_sequence'],seq)
            np.savez_compressed(refs/f'{name}.npz',design_mask=np.array([False]*194+[True]*len(seq)))
            cases.append({'id':name,'arm':arm['id'],'epitope':arm['epitope'],'protocol':arm['protocol'],
                'kind':'novel_design','target_sequence':manifest['target_sequence'],'binder_sequence':seq,
                'target_native_map':manifest['target_native_map'],'generated_source':str(path.relative_to(out)),
                'generation_target_fit_rmsd_A':rms})
    for control in manifest['controls']:
        case=dict(control,arm='control');cases.append(case)
        for ext in ('.cif','.npz','.pdb'):(refs/(case['id']+ext)).write_bytes((source/(case['id']+ext)).read_bytes())
    (out/'cases.json').write_text(json.dumps(cases,indent=2))
    (out/'actual_feature_audit.json').write_text(json.dumps(audit,indent=2))
    for seed in manifest['refolding']['seeds']:
        dest=out/f'seed-{seed}';dest.mkdir()
        for case in cases:
            for ext in ('.cif','.npz'):(dest/(case['id']+ext)).write_bytes((refs/(case['id']+ext)).read_bytes())
        (dest/'expected.json').write_text(json.dumps({c['id']:c for c in cases},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('operation',choices=['preflight','stage','assemble']);p.add_argument('output',type=Path)
    p.add_argument('--arm');p.add_argument('--stage',choices=['design','inverse_folding']);a=p.parse_args()
    if a.operation=='preflight':preflight(a.output)
    elif a.operation=='assemble':assemble(a.output)
    else:run_stage(a.output,a.arm,a.stage)
