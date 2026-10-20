"""Plot control confidence and experimental contact recovery from saved results."""
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
    conditions = [
        ("Pilot settings", ROOT / "reports/controls/results.json", "#4169a5", -.12),
        ("Target MSA added", ROOT / "reports/controls_msa/results.json", "#c06b27", .12),
    ]
    names = ["3ds18_current", "3ds18_crystal", "atv_current", "atv_permuted",
             "3ds18_scrambled", "caplacizumab_current"]
    labels = ["3DS18\n11QC crop", "3DS18\n6WRV crop", "ATV Fc\n11QC crop",
              "ATV Fc\npermuted crop", "Scrambled\n11QC crop", "Unrelated VHH\n11QC crop"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for title, path, color, offset in conditions:
        if not path.exists():
            continue
        result = json.loads(path.read_text())
        for i, name in enumerate(names):
            rows = [r for r in result["rows"] if r["id"] == name]
            for j, row in enumerate(rows):
                x = i+offset+(j-.5)*.035
                axes[0].scatter(x, row["selected_interface_iptm"], color=color, s=45,
                                marker="o" if j == 0 else "D",
                                label=title if i == 0 and j == 0 else None)
                if "pose" in row:
                    axes[1].scatter(x, row["pose"]["native_contact_recall"], color=color,
                                    s=45, marker="o" if j == 0 else "D")
    axes[0].axhline(.6, color="#999999", linestyle="--", linewidth=1)
    axes[0].text(5.45, .615, "Proposed 0.6 cutoff", ha="right", color="#666666", fontsize=9)
    axes[0].set_title("Confidence-selected model: interface ipTM")
    axes[1].set_title("Fraction of experimental residue contacts recovered")
    for ax in axes:
        ax.set_ylim(-.04, 1.04)
        ax.set_xticks(range(6), labels, fontsize=9)
        ax.grid(axis="y", alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    axes[1].text(4.5, .45, "No experimental pose\nfor negative controls",
                 ha="center", color="#666666", fontsize=9)
    fig.suptitle("TfR1 calibration: confidence and experimental-interface recovery", fontsize=15)
    fig.supxlabel("Two seeds per condition; five diffusion samples per seed. Published controls are not new designs.", fontsize=10)
    out = ROOT / "reports/controls"
    fig.savefig(out / "calibration_overview.png", dpi=180)
    fig.savefig(out / "calibration_overview.pdf")


if __name__ == "__main__":
    main()
