"""Static scientific figures for matched calibration and the RMSD frame audit."""
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "outputs/matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    out = ROOT / "reports/controls_matched"
    data = json.loads((out / "results.json").read_text())
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    binders = ["3ds18", "2ds25_5", "scrambled", "caplacizumab"]
    labels = ["3DS18\npositive", "2DS25.5\nheld-out positive",
              "Scrambled\npresumed negative", "Unrelated VHH\npresumed negative"]
    colors = {"6wrv": "#2874a6", "6wrw": "#c76c25"}
    for ax, key, title in [
        (axes[0, 0], "ipsae_min", "Interface confidence: minimum directional ipSAE"),
        (axes[0, 1], "contact_pae_both_directions_A", "PAE over contacting residue pairs (Å)"),
    ]:
        for j, binder in enumerate(binders):
            for k, receptor in enumerate(colors):
                rows = [r for r in data["rows"] if r["binder_id"] == binder
                        and r["receptor_id"] == receptor]
                rows = [r for r in rows if r["scores"][key] is not None]
                x = j + (k-.5)*.28 + np.linspace(-.045, .045, len(rows))
                ax.scatter(x, [r["scores"][key] for r in rows], color=colors[receptor],
                           s=26, alpha=.8, label=receptor.upper() if j == 0 else None)
        ax.set_xticks(range(4), labels, fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.set_ylim(bottom=0)
        if key == "ipsae_min":
            ax.set_ylim(bottom=-.008)
        ax.legend(frameon=False, fontsize=9)
    ax = axes[1, 0]
    for receptor in colors:
        for binder, marker in [("3ds18", "o"), ("2ds25_5", "^")]:
            rows = [r for r in data["rows"] if r["binder_id"] == binder
                    and r["receptor_id"] == receptor]
            ax.scatter([r["scores"]["ipsae_min"] for r in rows],
                       [r["pose"]["native_contact_recall"] for r in rows],
                       c=colors[receptor], marker=marker, s=32, alpha=.8,
                       label=f"{binder} / {receptor.upper()}")
    ax.set(xlabel="Minimum directional ipSAE", ylabel="Native contact recall",
           title="Positive controls: score versus interface recovery", ylim=(-.03, 1.03))
    ax.legend(frameon=False, fontsize=8)
    old = json.loads((ROOT / "reports/controls/all_samples_audit.json").read_text())[0]
    ax = axes[1, 1]
    for row in old["rows"]:
        if row["id"] == "3ds18_crystal" and row["confidence_selected"]:
            ds = row["residue_displacements"]
            ax.plot([d["binder_position"] for d in ds],
                    [d["CA_error_after_target_fit_A"] for d in ds],
                    label=f"Baseline seed {row['seed']}", linewidth=1.6)
    ax.set(xlabel="3DS18 residue index", ylabel="Cα displacement (Å)",
           title="Same receptor alignment: discrepancy at N terminus")
    ax.legend(frameon=False, fontsize=9)
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=.15)
    fig.suptitle("TfR1 matched calibration — all samples, fixed scoring settings", fontsize=15)
    fig.supxlabel("Ten predictions per binder/receptor. Presumed negatives are not measured nonbinders.",
                  fontsize=10)
    fig.savefig(out / "matched_calibration.png", dpi=180)
    fig.savefig(out / "matched_calibration.pdf")


if __name__ == "__main__":
    main()
