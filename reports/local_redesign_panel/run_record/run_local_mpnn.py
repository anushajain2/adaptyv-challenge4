"""CPU ProteinMPNN sampling with verified fixed residues and two-backbone ranking."""
import copy
import json
import random
import subprocess
import sys
import time
import numpy as np
from novelty_search import ROOT, sha, save

ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'


def sequons(seq):
    return {i for i in range(len(seq)-2) if seq[i] == 'N' and seq[i+1] != 'P' and seq[i+2] in 'ST'}


def eligible_sequence(seq, family):
    original = family['sequence']
    assert len(seq) == len(original)
    changes = [i+1 for i,(a,b) in enumerate(zip(original,seq)) if a != b]
    assert set(changes) <= set(family['mutable_positions']), 'Fixed residue mutated'
    assert seq[family['histidine_position']-1] == 'H'
    assert all(seq[p-1] not in 'CHX' for p in changes)
    return 1 <= len(changes) <= 4 and not (sequons(seq)-sequons(original))


def select_diverse(rows, maximum=6):
    selected = []
    for row in sorted(rows, key=lambda r:(r['ensemble_nll'], r['sequence'])):
        if all(sum(a!=b for a,b in zip(row['sequence'],other['sequence'])) >= 2 for other in selected):
            selected.append(row)
        if len(selected) == maximum: break
    return selected


def main():
    import torch
    torch.set_num_threads(2)
    vendor = ROOT/'vendor/ProteinMPNN'; sys.path.insert(0,str(vendor))
    from protein_mpnn_utils import parse_PDB, tied_featurize, ProteinMPNN, _scores
    plan_path = ROOT/'data/local_redesign/plan.json'; plan = json.loads(plan_path.read_text()); cfg = plan['proteinmpnn']
    assert subprocess.check_output(['git','-C',str(vendor),'rev-parse','HEAD'],text=True).strip() == cfg['commit']
    weights = vendor/cfg['model']; assert sha(weights) == cfg['weights_sha256']
    checkpoint = torch.load(weights, map_location='cpu', weights_only=True)
    model = ProteinMPNN(num_letters=21,node_features=128,edge_features=128,hidden_dim=128,
                        num_encoder_layers=3,num_decoder_layers=3,augment_eps=0,k_neighbors=checkpoint['num_edges'])
    model.load_state_dict(checkpoint['model_state_dict']); model.eval()
    out = ROOT/'reports/local_redesign/mpnn'; out.mkdir(exist_ok=False)
    started = time.monotonic(); result = {'plan_sha256':sha(plan_path),'script_sha256':sha(__file__),
        'torch_version':torch.__version__,'device':'cpu','families':[]}

    def features(source, family, batch):
        path = ROOT/source['path']; assert sha(path) == source['sha256']
        parsed = parse_PDB(str(path), input_chain_list=['A','B'])[0]
        assert parsed['seq_chain_B'] == family['sequence'] and len(parsed['seq_chain_A']) == 194
        name = parsed['name']; fixed = [i for i in range(1,len(family['sequence'])+1) if i not in family['mutable_positions']]
        f = tied_featurize([copy.deepcopy(parsed) for _ in range(batch)],torch.device('cpu'),{name:(['B'],['A'])},{name:{'B':fixed,'A':list(range(1,195))}})
        assert all(order == ['B','A'] for order in f[6])
        expected = np.zeros(f[1].shape[1]); expected[np.array(family['mutable_positions'])-1] = 1
        assert np.all(f[2].numpy()==1) and np.all((f[4]*f[10]).numpy()==expected)
        return f

    with torch.no_grad():
        for family in plan['families']:
            samples = []; unique = set()
            for source, seed in zip(family['models'], cfg['sampling_seeds']):
                torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
                f = features(source,family,4)
                X,S,mask,_,chain_M,chain_enc = f[:6]
                for temperature in cfg['temperatures']:
                    for batch in range(cfg['samples_per_temperature_per_backbone']//4):
                        sampled = model.sample(X,torch.randn(chain_M.shape),S,chain_M,chain_enc,f[12],mask=mask,
                            temperature=temperature,omit_AAs_np=np.array([a in 'CHX' for a in ALPHABET],dtype=np.float32),
                            bias_AAs_np=np.zeros(21,dtype=np.float32),chain_M_pos=f[10],omit_AA_mask=f[11],
                            pssm_coef=f[15],pssm_bias=f[16],pssm_multi=0.,pssm_log_odds_flag=False,pssm_bias_flag=False,bias_by_res=f[18])['S']
                        assert torch.equal(sampled[:,len(family['sequence']):],S[:,len(family['sequence']):]), 'Target sequence changed'
                        for encoded in sampled:
                            seq = ''.join(ALPHABET[i] for i in encoded[:len(family['sequence'])].tolist())
                            eligible = eligible_sequence(seq,family)
                            samples.append({'sequence':seq,'eligible':eligible,'source_seed':source['seed'],'sampling_seed':seed,'temperature':temperature,'batch':batch})
                            if eligible: unique.add(seq)
                print(f"Sampled {family['id']} source {source['seed']}: {len(unique)} unique eligible",flush=True)
            scored = []
            # Common decoding orders across sequences; report baseline on exactly the same masks.
            for seq in sorted(unique | {family['sequence']}):
                values = []
                for source in family['models']:
                    f = features(source,family,1); X,S,mask,_,chain_M,chain_enc = f[:6]
                    S = S.clone(); S[0,:len(seq)] = torch.tensor([ALPHABET.index(a) for a in seq])
                    scores = []
                    for seed in cfg['score_decoding_seeds']:
                        torch.manual_seed(seed)
                        log_probs = model(X,S,mask,chain_M*f[10],f[12],chain_enc,torch.randn(chain_M.shape))
                        scores.append(float(_scores(S,log_probs,mask*chain_M*f[10]).item()))
                    values.append({'source_seed':source['seed'],'nll_by_decoding_seed':scores})
                scored.append({'sequence':seq,'ensemble_nll':float(np.mean([r['nll_by_decoding_seed'] for r in values])),
                               'scores':values,'mutations':[f'{a}{i+1}{b}' for i,(a,b) in enumerate(zip(family['sequence'],seq)) if a!=b]})
            chosen = select_diverse([r for r in scored if r['sequence'] != family['sequence']])
            for i,row in enumerate(chosen,1): row['id'] = f"{family['id']}_mpnn{i:02d}"
            record = {'id':family['id'],'samples':samples,'scored_sequences':scored,'selected':chosen,
                      'baseline':next(r for r in scored if r['sequence']==family['sequence'])}
            save(out/f"{family['id']}.json",record); result['families'].append(record)
            save(out/'results.json',result)
            print(f"Selected {len(chosen)} for {family['id']}",flush=True)
    result['elapsed_seconds'] = time.monotonic()-started;save(out/'results.json',result)
    print(json.dumps({'selected':sum(len(f['selected']) for f in result['families']),'elapsed_seconds':result['elapsed_seconds']}))


if __name__ == '__main__': main()
