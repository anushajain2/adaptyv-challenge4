"""Compact scientific figure for the prospective 6WRW comparison."""
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR",str(ROOT/"outputs/matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    folder = ROOT/"reports/generation_comparison"
    result = json.loads((folder/"results.json").read_text())
    rows = [r for r in result["rows"] if r["confidence_selected"]]
    fig,axes = plt.subplots(1,3,figsize=(14,4.5),layout="constrained")
    groups = [("unconstrained", "Unconstrained", "#2166ac", "o"),
              ("sheet", "Sheet constraint", "#b35806", "s"),
              ("positive", "2DS25.5 control", "#1b7837", "^")]
    for group,label,color,marker in groups:
        subset = [r for r in rows if (r["arm"]==group if group!="positive" else r["kind"]=="positive")]
        axes[0].scatter([r["metrics"]["design_to_target_iptm"] for r in subset],
                        [r["pose"]["binder_CA_fold_rmsd_A"] for r in subset],
                        color=color,marker=marker,label=label,s=35,alpha=.8)
        axes[1].scatter([r["pose"]["native_contact_recall"] for r in subset],
                        [r["pose"]["binder_interface_CA_rmsd_after_target_fit_A"] for r in subset],
                        color=color,marker=marker,s=35,alpha=.8)
        axes[2].scatter([r["scores"]["ipsae_min"] for r in subset],
                        [r["scores"]["contact_pae_both_directions_A"] for r in subset],
                        color=color,marker=marker,s=35,alpha=.8)
    negatives = [r for r in rows if r["kind"]=="negative_hypothesis"]
    axes[2].scatter([r["scores"]["ipsae_min"] for r in negatives],
                    [r["scores"]["contact_pae_both_directions_A"] for r in negatives],
                    color="#777777",marker="x",label="Presumed negatives",s=40)
    axes[0].axvline(.5,color="gray",ls="--",lw=1)
    axes[0].axhline(2,color="gray",ls="--",lw=1)
    axes[0].set(xlabel="Interface ipTM",ylabel="Binder self-fit Cα RMSD (Å)",title="Fold + confidence")
    axes[0].legend(fontsize=8,loc="best")
    axes[1].axvline(.5,color="gray",ls="--",lw=1)
    axes[1].axhline(2,color="gray",ls="--",lw=1)
    axes[1].set(xlabel="Reference contact recall",ylabel="Interface Cα RMSD after receptor fit (Å)",
                title="Interface placement",xlim=(-.03,1.03))
    axes[2].set(xlabel="Minimum directional ipSAE (10 Å / 10 Å)",
                ylabel="Contact-pair PAE (Å)",title="Supporting scores (no pass threshold)")
    axes[2].legend(fontsize=8,loc="best")
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(alpha=.15)
    fig.suptitle("6WRW: 20 novel designs, two refolding seeds each",fontsize=15)
    fig.supxlabel("Confidence-selected model per seed. Reference = generated pose for novel designs; crystal pose for the positive control.\nThese are computational screens, not experimental binding measurements.",fontsize=9)
    fig.savefig(folder/"comparison.png",dpi=180)
    fig.savefig(folder/"comparison.pdf")


if __name__ == "__main__": main()
