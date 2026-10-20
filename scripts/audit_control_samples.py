"""Audit every saved diffusion sample; verify atom mapping against the exported CIF."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
from Bio.PDB import Atom, Superimposer

from analyze_controls import contact_pairs, pose_metrics
from prepare_controls import ROOT, geometry, read_structure


def restore_predicted_atoms(structure, coords, atom_to_token):
    """Restore predicted (not imputed) sidechains dropped by input resolved masks.

    The pinned upstream parser orders standard amino-acid atoms by const.ref_atoms.
    Verify the coordinates of EVERY retained atom before adding any missing atom.
    """
    order = json.loads((ROOT / "configs/boltzgen_atom_order.json").read_text())["ref_atoms"]
    token = 0
    additions = []
    for chain in structure[0]:
        for residue in chain:
            names = order[residue.resname]
            indices = np.flatnonzero(atom_to_token[:, token])
            assert len(indices) == len(names), (residue.resname, len(indices), names)
            positions = {name: coords[index] for name, index in zip(names, indices)}
            for atom in residue:
                assert atom.name in positions
                assert np.max(np.abs(atom.coord - positions[atom.name])) <= .002, (
                    chain.id, residue.id, atom.name, "canonical atom order mismatch")
            additions.extend((residue, name, positions[name])
                             for name in names if name not in residue)
            token += 1
    assert token == atom_to_token.shape[1], "Only complete, standard-protein sequences supported"
    for residue, name, xyz in additions:
        residue.add(Atom.Atom(name, xyz.copy(), 0., 1., " ", f"{name:>4}", 0, element=name[0]))
    return len(additions)


def sample_structures(folder, name, restore_sidechains=False):
    """Use an independently checked selected-CIF atom order, never guess an order."""
    with np.load(folder / "fold_out_npz" / f"{name}.npz", allow_pickle=False) as data:
        coords = data["coords"].copy()
        mask = data["atom_resolved_mask"][0].astype(bool)
        atom_to_token = data["atom_to_token"][0].copy() if restore_sidechains else None
        selected = int(np.argmax(.8 * data["iptm"] + .2 * data["ptm"]))
        metrics = {k: data[k].tolist() for k in
                   ("iptm", "ptm", "design_to_target_iptm", "interaction_pae",
                    "min_design_to_target_pae", "design_ptm")}
    structure = read_structure(folder / "refold_cif" / f"{name}.cif")
    atoms = list(structure.get_atoms())
    exported = np.array([a.coord for a in atoms])
    assert exported.shape == coords[selected, mask].shape, name
    error = float(np.max(np.abs(exported - coords[selected, mask])))
    assert error <= .002, (name, "atom order mismatch", error)
    for i in range(len(coords)):
        sample = copy.deepcopy(structure)
        for a, xyz in zip(sample.get_atoms(), coords[i, mask]):
            a.coord = xyz.copy()
        if restore_sidechains:
            restore_predicted_atoms(sample, coords[i], atom_to_token)
        yield i, sample, {k: v[i] for k, v in metrics.items()}, i == selected


def residue_displacements(predicted, reference):
    indices = [r.id[1] for r in reference[0]["A"] if "CA" in r
               and (" ", r.id[1], " ") in predicted[0]["A"]]
    sup = Superimposer()
    sup.set_atoms([reference[0]["A"][i]["CA"] for i in indices],
                  [predicted[0]["A"][i]["CA"] for i in indices])
    rot, shift = sup.rotran
    interface = {b for _, b in contact_pairs(reference)}
    return [{"binder_position": r.id[1], "native_interface": r.id[1] in interface,
             "CA_error_after_target_fit_A": float(np.linalg.norm(
                 predicted[0]["B"][r.id]["CA"].coord @ rot + shift - r["CA"].coord))}
            for r in reference[0]["B"] if "CA" in r
            and r.id in predicted[0]["B"]]


def audit(run, restore_sidechains=False):
    manifest = json.loads((run / "inputs/manifest.json").read_text())
    cases = {c["id"]: c for c in manifest["cases"]}
    full = read_structure(ROOT / "data/raw/11QC.cif")
    rows = []
    for seed in manifest["refolding"]["seeds"]:
        for case in cases.values():
            ref_id = case.get("reference_id", {"3ds18_current": "3ds18_crystal"}.get(
                case["id"], case["id"]))
            reference = read_structure(ROOT / cases[ref_id]["pdb"])
            for index, sample, metrics, selected in sample_structures(
                    run / f"seed-{seed}", case["id"], restore_sidechains):
                mapping = {int(k): v for k, v in case["target_native_map"].items()
                           if (" ", int(k), " ") in sample[0]["A"]}
                g = geometry(sample, mapping, full)
                row = {"id": case["id"], "kind": case["kind"], "seed": seed,
                       "sample": index, "confidence_selected": selected,
                       "metrics": metrics, "geometry": g,
                       "full_context_fit_valid": g["target_fit_rmsd_A"] <= 2}
                if case["kind"].startswith("positive"):
                    row["pose"] = pose_metrics(sample, reference)
                    row["interface_recovered"] = (
                        row["pose"]["binder_interface_CA_rmsd_after_target_fit_A"] <= 2
                        and row["pose"]["native_contact_recall"] >= .5)
                    row["residue_displacements"] = residue_displacements(sample, reference)
                rows.append(row)
    summary = {}
    for name, case in cases.items():
        subset = [r for r in rows if r["id"] == name]
        summary[name] = {
            "samples": len(subset),
            "interface_recovered": sum(r.get("interface_recovered", False) for r in subset),
            "selected_interface_recovered": sum(
                r.get("interface_recovered", False) and r["confidence_selected"] for r in subset),
            "interface_pae_range_A": [min(r["metrics"]["interaction_pae"] for r in subset),
                                      max(r["metrics"]["interaction_pae"] for r in subset)],
        }
    return {"run": run.name, "mapping_check": "Every selected CIF matches masked NPZ coordinates",
            "summary": summary, "rows": rows}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("runs", nargs="+", type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    result = [audit(run) for run in args.runs]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2))
    for run in result:
        print(run["run"], json.dumps(run["summary"], indent=2))


if __name__ == "__main__":
    main()
