"""Read real BoltzGen outputs, preserve raw metrics and normalize refolded complexes."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from prepare_target import ROOT, align_map, save_pdb, write_csv
from direction_check import full_context, TargetFitError

KD = dict(zip("ACDEFGHIKLMNPQRSTVWY",
              [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))

def parent_quality(metrics, clashes):
    try:
        return (clashes==0 and float(metrics["design_to_target_iptm"])>=0.5
                and float(metrics["bb_rmsd"])<2.5)
    except (KeyError,TypeError,ValueError):
        return False


def normalize(path,out,allow_context_failure=False):
    structure=MMCIFParser(QUIET=True).get_structure("design",path)
    prepared=MMCIFParser(QUIET=True).get_structure("target",ROOT/"configs/human_apical.cif")
    target_seq=seq1("".join(r.resname for r in prepared[0]["A"] if r.id[0]==" "))
    chains={c.id:seq1("".join(r.resname for r in c if r.id[0]==" ")) for c in structure[0]}
    target=[c for c,s in chains.items() if s==target_seq]
    if len(target)!=1 or len(chains)!=2:
        raise ValueError(f"Expected one exact target and one binder; got chain lengths { {c:len(s) for c,s in chains.items()} }")
    binder=next(c for c in chains if c!=target[0])
    normalized=copy.deepcopy(structure)
    for c in list(normalized[0]):
        normalized[0].detach_child(c.id)
    for old,new in [(target[0],"A"),(binder,"B")]:
        c=copy.deepcopy(structure[0][old])
        c.id=new
        residues=[r for r in c if r.id[0]==" "]
        for r in list(c):
            c.detach_child(r.id)
        for i,r in enumerate(residues,1):
            r.id=(" ",i," ")
            c.add(r)
        normalized[0].add(c)
    save_pdb(normalized,out)
    try:
        rms,clashes=full_context(normalized,out.with_name(out.stem+"_full_context.pdb"))
    except TargetFitError as e:
        if not allow_context_failure:
            raise
        # A failed scientific screen must not discard an otherwise readable design.
        rms,clashes=e.rms,None
    return chains[binder],rms,clashes


def cdr_positions(sequence):
    """Align to the framework; retain CDR spans including inserted residues."""
    scaffold=MMCIFParser(QUIET=True,auth_chains=False).get_structure("scaffold",ROOT/"configs/7eow.cif")
    reference=seq1("".join(r.resname for r in scaffold[0]["B"] if r.id[0]==" "))
    _,mapping=align_map(reference,sequence)
    # Anchoring on framework flanks handles variable-length CDRs.
    allowed=[]
    for left,right in [(25,35),(51,60),(97,119)]:
        if left not in mapping or right not in mapping:
            raise ValueError("Cannot map nanobody CDR framework boundaries")
        allowed.extend(range(mapping[left]+1,mapping[right]))
    return allowed


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("pilot",type=Path)
    p.add_argument("--out",type=Path,default=ROOT/"outputs/analysis")
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    metric_files=sorted(args.pilot.rglob("aggregate_metrics*.csv"))
    if not metric_files:
        raise SystemExit("No real aggregate_metrics*.csv found; no results fabricated.")
    rows=[]
    failures=[]
    manifest=[]
    seen=set()
    for metrics in metric_files:
        protocol=next((p for p in ("protein-anything","nanobody-anything") if p in str(metrics)), "unknown")
        with metrics.open() as f:
            for row in csv.DictReader(f):
                filename=row.get("file_name") or (row.get("id","")+".cif")
                path=metrics.parent/"refold_cif"/filename
                if path.suffix!=".cif":
                    path=path.with_suffix(".cif")
                if not path.exists():
                    failures.append({"metrics":str(metrics),"file":filename,"error":"Refolded complex missing"})
                    continue
                digest=hashlib.sha256(path.read_bytes()).hexdigest()
                if (protocol,digest) in seen:
                    continue
                seen.add((protocol,digest))
                name=f"{protocol}_{path.stem}"
                try:
                    seq,rms,clashes=normalize(path,args.out/f"{name}.pdb",allow_context_failure=True)
                    positions=cdr_positions(seq) if protocol=="nanobody-anything" else list(range(1,len(seq)+1))
                    rows.append({"design":name,"protocol":protocol,"sequence":seq,
                                 "length":len(seq),"mean_kyte_doolittle":sum(KD[a] for a in seq)/len(seq),
                                 "target_fit_rmsd_A":rms,"full_dimer_clashes":clashes,
                                 "geometry_status":"target_fit_failed" if clashes is None else
                                  "receptor_clashes" if clashes else "clash_free",
                                 "design_to_target_iptm":row.get("design_to_target_iptm"),
                                 "backbone_refolding_rmsd_A":row.get("bb_rmsd"),
                                 "min_design_to_target_pae_A":row.get("min_design_to_target_pae"),
                                 "parent_quality_pass":parent_quality(row,clashes),
                                 "raw_metrics_json":json.dumps(row,sort_keys=True)})
                    manifest.append({"design":name,"complex":str(args.out/f"{name}.pdb"),
                                     "sequence":seq,
                                     "protocol":protocol,"allowed_histidine_positions":positions,
                                     "full_dimer_clashes":clashes,
                                     "parent_quality_pass":parent_quality(row,clashes),
                                     "metrics":row})
                except Exception as e:
                    failures.append({"file":str(path),"error":str(e)})
    if rows:
        write_csv(args.out/"pilot_designs.csv",rows)
        (args.out/"pilot_designs.fasta").write_text("".join(f">{r['design']}\n{r['sequence']}\n" for r in rows))
    (args.out/"scan_manifest.json").write_text(json.dumps(manifest,indent=2))
    (args.out/"read_errors.json").write_text(json.dumps(failures,indent=2))
    print(json.dumps({"readable_designs":len(rows),"errors":len(failures)}))
    if not rows:
        raise SystemExit("No readable designs; inspect read_errors.json")


if __name__=="__main__":
    main()
