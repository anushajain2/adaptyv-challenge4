"""Sequence identity and descriptive backbone-angle composition of retained leads."""
import hashlib,itertools,json,math
from Bio.PDB.Polypeptide import PPBuilder
from prepare_controls import ROOT,read_structure
from analyze_generation_comparison import aligned_identity


def main():
    out=ROOT/'reports/next_round';result=json.loads((ROOT/'reports/expansion_validation/results.json').read_text())
    leads=[c for c in result['candidate_summaries'] if c['all_selected_seeds_full_screen']]
    hist=json.loads((ROOT/'reports/his_confirmation/results.json').read_text())
    leads.append(next(c for c in hist['candidate_summaries'] if c['id']=='apical_short_008_A41H'))
    rows=[]
    for c in leads:
        models=[]
        for r in c['selected_models']:
            path=ROOT/r['raw_scoring_output']/'model.pdb';s=read_structure(path)
            angles=[(math.degrees(phi),math.degrees(psi)) for pp in PPBuilder().build_peptides(s[0]['B']) for phi,psi in pp.get_phi_psi_list() if phi is not None and psi is not None]
            models.append({'seed':r['seed'],'backbone_angles_evaluated':len(angles),'alpha_region_fraction':sum(-100<=p<=-30 and -80<=q<=-5 for p,q in angles)/len(angles),'beta_region_fraction':sum(-180<=p<=-80 and (q>=80 or q<=-120) for p,q in angles)/len(angles),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        rows.append({'id':c['id'],'sequence':c['sequence'],'length':len(c['sequence']),'models':models,'pH_evidence':'A41H direction hypothesis; strict gate failed in one pose' if c['id'].endswith('A41H') else 'No pH-direction evaluation for this new parent'})
    comparisons=[{'a':a['id'],'b':b['id'],'global_sequence_identity':aligned_identity(a['sequence'],b['sequence'])} for a,b in itertools.combinations(rows,2)]
    report={'leads':rows,'pairwise_sequence_identity':comparisons,'limitations':'Backbone-angle bins are descriptive, not DSSP assignments or independent fold-class clustering. Sequence diversity does not establish structural-class or method-family diversity. All candidates originate from BoltzGen and this same apical campaign.'}
    (out/'lead_diversity.json').write_text(json.dumps(report,indent=2))
    (out/'retained_human_leads.fasta').write_text(''.join(f">{r['id']} COMPUTATIONAL_LEAD_NOT_VALIDATED_BINDING\n{r['sequence']}\n" for r in rows))
    print('Maximum pairwise sequence identity',max(c['global_sequence_identity'] for c in comparisons))
    for r in rows:print(r['id'],[round(m['alpha_region_fraction'],2) for m in r['models']])
if __name__=='__main__':main()
