"""Run direction checks on up to three clash-free parents per protocol."""
import argparse
import json
import math
from pathlib import Path
from direction_check import scan
from prepare_target import ROOT
from read_pilot import parent_quality


def metric(row,*names):
    for name in names:
        try:
            value=float(row["metrics"][name])
            if math.isfinite(value):
                return value
        except (KeyError,ValueError,TypeError):
            pass
    return -1.0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest",type=Path)
    p.add_argument("--out",type=Path,default=ROOT/"outputs/direction")
    p.add_argument("--parents-per-protocol",type=int,default=3)
    p.add_argument("--pairs",type=int,default=6)
    p.add_argument("--protocol",choices=["protein-anything","nanobody-anything"],
                   help="Process a completed modality while the other GPU batch continues.")
    p.add_argument("--exploratory",action="store_true",
                   help="Allow low-confidence parents for diagnostics, never qualified switch candidates.")
    args=p.parse_args()
    parents=json.loads(args.manifest.read_text())
    selected=[]
    for protocol in ("protein-anything","nanobody-anything"):
        if args.protocol and protocol!=args.protocol:
            continue
        group=[r for r in parents if r["protocol"]==protocol and r["full_dimer_clashes"]==0
               and (args.exploratory or parent_quality(r["metrics"],r["full_dimer_clashes"]))]
        group.sort(key=lambda r:metric(r,"design_to_target_iptm","design_iptm","iptm","iPTM","complex_iptm"),reverse=True)
        selected.extend(group[:args.parents_per_protocol])
    args.out.mkdir(parents=True,exist_ok=True)
    suffix=f"_{args.protocol}" if args.protocol else ""
    (args.out/f"selected_parents{suffix}.json").write_text(json.dumps(selected,indent=2))
    if not selected:
        (args.out/f"pilot_direction_summary{suffix}.json").write_text(json.dumps({
            "status":"no_qualified_parents",
            "parents_evaluated":sum(not args.protocol or r["protocol"]==args.protocol for r in parents),
            "criteria":"zero receptor clashes, design_to_target_iptm >=0.5, bb_rmsd <2.5 A",
            "direction_models_run":0},indent=2))
        print("No qualified parents; a valid negative pilot result. No pKa candidates fabricated.")
        return
    results=[]
    for row in selected:
        models=scan(Path(row["complex"]),args.out/row["design"],args.pairs,
                    allowed=set(row["allowed_histidine_positions"]))
        results.append({"design":row["design"],"models":len(models),
                        "provisional_models":sum(m["provisional_candidate"] for m in models)})
    (args.out/f"pilot_direction_summary{suffix}.json").write_text(json.dumps(results,indent=2))


if __name__=="__main__":
    main()
