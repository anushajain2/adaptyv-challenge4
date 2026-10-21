"""Consolidate real pilot pair models, separating sign evidence from strict triage."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from prepare_target import ROOT, write_csv
from read_pilot import parent_quality


def terminal_distance(path,positions):
    """Distance from engineered His heavy atoms to artificial target termini."""
    import numpy as np
    from direction_check import read_structure
    structure=read_structure(path)
    ends=np.array([a.coord for i in [1,194] for a in structure[0]["A"][i]
                   if a.element!="H"])
    histidines=np.array([a.coord for i in positions for a in structure[0]["B"][i]
                        if a.element!="H"])
    return float(np.linalg.norm(histidines[:,None]-ends[None],axis=2).min())


def number(row,key):
    value=row.get(key)
    return float(value) if value not in (None,"") else None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--direction",type=Path,default=ROOT/"outputs/direction")
    p.add_argument("--manifest",type=Path,default=ROOT/"outputs/analysis/scan_manifest.json")
    p.add_argument("--out",type=Path,default=ROOT/"reports/pilot")
    args=p.parse_args()
    parents={r["design"]:r for r in json.loads(args.manifest.read_text())}
    pairs=[]
    counts=Counter()
    evaluated={p.parent.name for p in args.direction.glob("*/enumerated_pairs.json")}
    evaluated.update(p.parent.name for p in args.direction.glob("*/direction_results.csv"))
    if not evaluated <= parents.keys():
        raise ValueError("Direction results include parents absent from the real pilot manifest")
    parent_count=len(evaluated)
    for path in sorted(args.direction.glob("*/direction_results.csv")):
        parent=path.parent.name
        if parent not in parents:
            raise ValueError(f"Results are not traceable to a real pilot parent: {parent}")
        with path.open() as f:
            rows=list(csv.DictReader(f))
        groups=defaultdict(list)
        for row in rows:
            groups[row["positions"]].append(row)
            counts[f"model_{row['status']}"]+=1
            counts[f"direction_{row['direction']}"]+=1
        for positions,models in groups.items():
            seeds={m["seed"] for m in models}
            if len(seeds)!=len(models):
                raise ValueError(f"Duplicate model seeds: {parent} / {positions}")
            ok=[m for m in models if m["status"]=="ok"]
            complete=len(ok)==len(models) and len(seeds)>=2
            lower=complete and all(number(m,"shift_1")<0 and number(m,"shift_2")<0 for m in ok)
            direction_pass=complete and all(m["provisional_candidate"]=="True" for m in ok)
            quality=parent_quality(parents[parent]["metrics"],parents[parent].get("full_dimer_clashes"))
            strict=direction_pass and quality
            pos=list(map(int,positions.split(",")))
            distances=[terminal_distance(model,pos) for m in ok
                       if (model:=path.parent/m.get("variant","")/"complex.pdb").exists()]
            sequence=list(parents[parent]["sequence"]) if "sequence" in parents[parent] else None
            # Older manifests keep the full chain in the raw upstream metrics.
            if sequence is None:
                from direction_check import read_structure
                from Bio.SeqUtils import seq1
                s=read_structure(parents[parent]["complex"])
                sequence=list(seq1("".join(r.resname for r in s[0]["B"] if r.id[0]==" ")))
            changes=[]
            for site in pos:
                changes.append(f"{sequence[site-1]}{site}H")
                sequence[site-1]="H"
            ratio=[number(m,"his_only_Kd6_over_Kd7_4") for m in ok]
            metrics=parents[parent]["metrics"]
            pairs.append({
                "variant":parent+"__"+"_".join(changes),"parent":parent,
                "protocol":parents[parent]["protocol"],"mutations":",".join(changes),
                "sequence":"".join(sequence),"models_successful":len(ok),
                "models_attempted":len(models),"seeds":",".join(sorted(seeds)),
                "both_pKas_lowered_all_models":lower,
                "direction_criteria_pass_all_models":direction_pass,
                "parent_quality_pass":quality,
                "provisional_candidate_all_models":strict,
                "min_predicted_Kd6_over_Kd7_4":min(ratio) if ratio else None,
                "max_predicted_Kd6_over_Kd7_4":max(ratio) if ratio else None,
                "least_negative_shift_1":max(number(m,"shift_1") for m in ok) if ok else None,
                "least_negative_shift_2":max(number(m,"shift_2") for m in ok) if ok else None,
                "parent_iptm":metrics.get("design_to_target_iptm",metrics.get("design_iptm")),
                "parent_backbone_rmsd":metrics.get("bb_rmsd"),
                "min_distance_to_target_terminal_residues_A":min(distances) if distances else None,
                "pka_target_context":"intact_194_residue_apical_crop",
                "evidence_csv":str(path),
                "interpretation":("Provisional direction hypothesis; unvalidated mutant affinity" if strict else
                                  "Exploratory only; fails parent and/or direction criteria"),
            })
    pairs.sort(key=lambda r:(r["provisional_candidate_all_models"],
                             r["both_pKas_lowered_all_models"],
                             r["min_predicted_Kd6_over_Kd7_4"] or 0),reverse=True)
    args.out.mkdir(parents=True,exist_ok=True)
    if pairs:
        write_csv(args.out/"all_histidine_pairs.csv",pairs)
    qualified=[r for r in pairs if r["provisional_candidate_all_models"]]
    (args.out/"provisional_switches.fasta").write_text(
        "".join(f">{r['variant']}\n{r['sequence']}\n" for r in qualified))
    exploratory=[r for r in pairs if r["both_pKas_lowered_all_models"]
                 and not r["provisional_candidate_all_models"]]
    (args.out/"exploratory_pka_lowering.fasta").write_text(
        "".join(f">{r['variant']} evidence=EXPLORATORY parent_quality_pass={r['parent_quality_pass']}\n"
                f"{r['sequence']}\n" for r in exploratory))
    summary={
        "pilot_parents_analyzed":parent_count,"histidine_pairs":len(pairs),
        "model_counts":dict(counts),
        "pairs_with_both_pKas_lowered_in_all_models":sum(r["both_pKas_lowered_all_models"] for r in pairs),
        "pairs_passing_direction_criteria":sum(r["direction_criteria_pass_all_models"] for r in pairs),
        "provisional_candidates":len(qualified),
        "exploratory_pka_lowering_variants":len(exploratory),
        "strict_criteria":"Parent ipTM >=0.5, backbone refold RMSD <2.5 A, no receptor clashes; two successful seeds, both shifts <= -0.5, at least one free pKa 6-7, His-only Kd ratio >=2, no PROPKA incomplete-group warnings",
        "limitations":["No experimental binding evidence","Mutant affinity and folding not validated",
                       "His-only independent-site linkage model",
                       "PROPKA uses trimmed target with artificial termini; omitted receptor/glycans may affect pKa",
                       "Incomplete conformational sampling"],
    }
    (args.out/"direction_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
