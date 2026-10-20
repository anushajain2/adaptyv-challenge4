"""Plot measured pipeline outputs; no synthetic design data."""
import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1]/"outputs/.matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--designs",type=Path,required=True)
    p.add_argument("--pairs",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    designs=list(csv.DictReader(args.designs.open()))
    pairs=list(csv.DictReader(args.pairs.open())) if args.pairs.exists() else []
    complete=[r for r in pairs if int(r["models_successful"])>=2
              and r["models_successful"]==r["models_attempted"]]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(1,2,figsize=(12,5.8),layout="constrained")
    ax=axes[0]
    ax.add_patch(Rectangle((.5,0),.5,2.5,color="#DDF2E2",zorder=0))
    groups=[
        ("Initial miniproteins","#376DA5",lambda r:"protein-anything"==r["protocol"] and "conserved_326_347" in r["design"]),
        ("Initial nanobodies","#B06428",lambda r:"nanobody-anything"==r["protocol"]),
        ("Apical-tip miniproteins","#358063",lambda r:"apical_tip_basic_direction_probe" in r["design"]),
    ]
    for label,color,select in groups:
        rows=[r for r in designs if select(r)]
        if rows:
            ax.scatter([float(r["design_to_target_iptm"]) for r in rows],
                       [float(r["backbone_refolding_rmsd_A"]) for r in rows],
                       s=34,color=color,alpha=.8,label=f"{label} (n={len(rows)})")
    ax.axvline(.5,color="#666666",ls="--",lw=1)
    ax.axhline(2.5,color="#666666",ls="--",lw=1)
    ax.set(xlim=(0,1),ylim=(0,max(25,max(float(r["backbone_refolding_rmsd_A"]) for r in designs)*1.08)),
           xlabel="Interface ipTM (higher is better)",
           ylabel="Design-to-refold backbone RMSD (Å)",
           title="Parent pose quality")
    ax.legend(loc="upper right",frameon=False,fontsize=8)
    ax.text(.53,1.1,"Confidence gate*",fontsize=8,color="#266B40")
    ax=axes[1]
    for protocol,color in [("protein-anything","#376DA5"),("nanobody-anything","#B06428")]:
        rows=[r for r in complete if r["protocol"]==protocol]
        if rows:
            ax.scatter([float(r["least_negative_shift_1"]) for r in rows],
                       [float(r["least_negative_shift_2"]) for r in rows],
                       s=65,facecolors="none",edgecolors=color,linewidths=1.5,
                       label=f"{protocol.split('-')[0]} pairs (n={len(rows)})")
    xs=[float(r["least_negative_shift_1"]) for r in complete]
    ys=[float(r["least_negative_shift_2"]) for r in complete]
    low=min([-1.8]+xs+ys)-.2
    high=max([.4]+xs+ys)+.2
    ax.add_patch(Rectangle((low,low),-.5-low,-.5-low,color="#DDF2E2",zorder=0))
    for value,style in [(0,"-"),(-.5,"--")]:
        ax.axvline(value,color="#999999",ls=style,lw=.8)
        ax.axhline(value,color="#999999",ls=style,lw=.8)
    ax.set(xlim=(low,high),ylim=(low,high),
           xlabel="His 1: bound − free pKa",
           ylabel="His 2: bound − free pKa",
           title="Direction check: least negative shift across seeds")
    if complete:
        ax.legend(loc="upper right",frameon=False,fontsize=8)
    else:
        ax.text(.5,.5,"No complete two-model pairs",transform=ax.transAxes,ha="center")
    ax.text(.02,.02,"Negative shifts favor acidic release\nin the modeled independent-site approximation.",
            transform=ax.transAxes,fontsize=8)
    fig.suptitle("TfR1 days 1–3: computational pilot",fontsize=16)
    fig.supxlabel(
        "*Parent also requires zero full-receptor clashes. Shaded regions are partial gates, not proof of binding.\n"
        "pKa points require two successful models; PROPKA uses the trimmed target. Mutant affinity is unvalidated.",
        fontsize=9)
    args.out.mkdir(parents=True,exist_ok=True)
    for extension in ["png","pdf"]:
        fig.savefig(args.out/f"pilot_overview.{extension}",dpi=180)
    print(args.out/"pilot_overview.png")


if __name__=="__main__":
    main()
