"""Freeze a narrow, prospective 6WRW generation experiment before inference."""
import copy
import csv
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import yaml
from Bio.PDB import MMCIFIO
from Bio.PDB.SASA import ShrakeRupley

from analyze_controls import contact_pairs
from prepare_controls import ROOT, fit_mapped, read_structure, sequence
from prepare_target import add_polymer_metadata, write_csv

OUT = ROOT / "data/generation_comparison"


def prepare(out=OUT):
    out.mkdir(parents=True, exist_ok=False)
    ref = read_structure(ROOT / "data/controls_matched/2ds25_5_6wrw.pdb")
    full = read_structure(ROOT / "data/raw/11QC.cif")
    mapping = {i: i+189 for i in range(1, 195)}
    fitted, rms = fit_mapped(ref, mapping, full)
    assert rms < 2
    # Bound local coordinates in full-dimer context, for accessibility only.
    # This hybrid is not used for folding, energetic scoring, or PROPKA.
    hybrid = copy.deepcopy(full)
    for i, n in mapping.items():
        hybrid[0]["A"].detach_child((" ", n, " "))
        residue = fitted[0]["A"][i].copy()
        residue.id = (" ", n, " ")
        hybrid[0]["A"].add(residue)
    crop = copy.deepcopy(fitted)
    crop[0].detach_child("B")
    sasa = ShrakeRupley(n_points=200)
    sasa.compute(hybrid, level="R")
    sasa.compute(crop, level="R")
    boundaries = np.array([crop[0]["A"][1]["N"].coord, crop[0]["A"][194]["C"].coord])
    anchors = np.array([c[n]["ND2"].coord for c in full[0] for n in (251,317,727)])
    old = {int(r["human_uniprot"]): r for r in csv.DictReader(
        (ROOT/"reports/controls/spatial_target_map.csv").open())}
    rows = []
    for i, n in mapping.items():
        r = crop[0]["A"][i]
        xyz = np.array([a.coord for a in r if a.element != "H"])
        cut = float(np.linalg.norm(xyz[:, None]-boundaries, axis=2).min())
        glycan = float(np.linalg.norm(xyz[:, None]-anchors, axis=2).min())
        exposed = float(hybrid[0]["A"][n].sasa)
        excess = float(r.sasa-exposed)
        reasons = []
        if cut < 10: reasons.append("crop_cut_buffer")
        if glycan < 10: reasons.append("glycan_anchor_buffer")
        if exposed < 25: reasons.append("low_dimer_context_SASA")
        if excess > 20: reasons.append("crop_exposed_surface")
        rows.append({"native_position": n, "crop_position": i,
                     "human_aa": old[n]["human_aa"], "mouse_aa": old[n]["mouse_aa"],
                     "mouse_position": old[n]["mouse_uniprot"],
                     "conserved": old[n]["conserved"] == "True",
                     "bound_context_sasa_A2": exposed, "crop_sasa_A2": float(r.sasa),
                     "crop_excess_sasa_A2": excess, "crop_cut_distance_A": cut,
                     "glycan_anchor_distance_A": glycan,
                     "apo_spatial_pass": old[n]["spatial_screen_pass"] == "True",
                     "bound_spatial_pass": not reasons, "exclusions": ";".join(reasons)})
    hotspots = [208, 210, 211]
    lookup = {r["native_position"]: r for r in rows}
    assert all(lookup[n]["bound_spatial_pass"] and lookup[n]["apo_spatial_pass"] for n in hotspots)
    reference_region = sorted({a for a,b in contact_pairs(ref)})
    # Native sheet-edge contacts can have low total sidechain SASA. Do not mask
    # an experimentally observed contact solely because its total SASA is low.
    excluded = [r["crop_position"] for r in rows if
                any(reason != "low_dimer_context_SASA" for reason in r["exclusions"].split(";") if reason)
                or (r["bound_context_sasa_A2"] < 25 and r["crop_position"] not in reference_region)]
    target = copy.deepcopy(ref)
    target[0].detach_child("B")
    io = MMCIFIO()
    io.set_structure(target)
    io.save(str(out/"target_6wrw.cif"))
    add_polymer_metadata(out/"target_6wrw.cif", list(target[0]["A"]))
    for arm in ("unconstrained", "sheet"):
        protein = {"id": "B", "sequence": "88"}
        if arm == "sheet":
            protein["secondary_structure"] = {"sheet": "33..37"}
        spec = {"entities": [
            {"file": {"path": "target_6wrw.cif", "include": [{"chain": {"id": "A"}}],
                      "binding_types": [
                          {"chain": {"id": "A", "not_binding": ",".join(map(str, excluded))}},
                          {"chain": {"id": "A", "binding": ",".join(str(n-189) for n in hotspots)}}]}},
            {"protein": protein}]}
        (out/f"{arm}.yaml").write_text(yaml.safe_dump(spec, sort_keys=False))
    controls = []
    old_manifest = json.loads((ROOT/"data/controls_matched/manifest.json").read_text())
    for case in old_manifest["cases"]:
        if case["id"] not in ("2ds25_5_6wrw", "scrambled_6wrw", "caplacizumab_6wrw"):
            continue
        controls.append(case)
        for ext in (".pdb", ".cif", ".npz"):
            source = ROOT/"data/controls_matched"/(case["id"]+ext)
            (out/source.name).write_bytes(source.read_bytes())
    manifest = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Prospective exploratory generation; previous calibration gate remains failed",
        "boltzgen_commit": "247b9bbd8b68a60aba854c2968d6a0cddd21ad6d",
        "target_sequence": sequence(target[0]["A"]), "target_native_map": mapping,
        "receptor": "6WRW author chain A, native residues 190-383",
        "reference_binder": "2DS25.5: control only; no published binder sequence or coordinates condition generation",
        "hotspots_native": hotspots, "hotspots_crop": [n-189 for n in hotspots],
        "reference_contact_region_crop": reference_region,
        "sheet_constraint": {"binder_positions": list(range(33,38)),
            "rationale": "Reference backbone N/O contacts: target208-binder37, target210-binder35, target212-binder33. Only secondary structure is conditioned, not sheet pairing or strand register."},
        "arms": [{"id": arm, "num_designs": 10, "binder_length": 88,
                  "generation_seed": 20261021, "inverse_fold_seed": 20261022}
                 for arm in ("unconstrained", "sheet")],
        "refolding": {**old_manifest["refolding"], "seeds": [101,211]},
        "controls": controls,
        "screen": {
            "selected_sample": "max(0.8*iptm+0.2*ptm), fixed before geometric evaluation",
            "min_design_to_target_iptm": 0.50, "max_binder_CA_self_fit_rmsd_A": 2.0,
            "max_interface_CA_rmsd_after_target_fit_A": 2.0,
            "min_design_contact_recall": 0.5, "min_hotspots_contacted": 2,
            "min_contact_pairs_in_reference_region_fraction": 0.5,
            "max_full_context_fit_rmsd_A": 2.0, "max_clash_pairs_lt2A": 0,
            "replicates": "Confidence-selected model must pass in BOTH fresh seeds; all 5 samples retained per seed",
            "controls": "2DS25.5 selected model must recover fold<=2A, native-contact recall>=0.5, interface pose<=2A, ipTM>=0.5 in both seeds; presumed negatives must not reach 0.5 interface ipTM in selected models",
            "ranking_only": ["ipSAE at PAE10/distance10; both directions", "contact-pair PAE", "pLDDT", "hydropathy"],
            "families": "Exploratory complete-link clustering: >=30% identity OR equal-length self-fit CA RMSD<=3A joins a family; inspect before any expansion",
            "decision": "No automatic expansion. >=2 distinct repeatable families supports modest expansion; one requires focused checking; zero requires strategy review, not lower cutoffs.",
            "limits": "Self-consistency is not binding validation. Low-score exclusions are not demonstrated nonbinders. No affinity or pH claims. Mouse compatibility remains untested."},
        "max_prediction_seconds": 3600, "experiment_cap_usd": 5, "total_cap_usd": 30,
        "spatial_context": {"fit_rmsd_A": rms, "sasa_points": 200,
            "note": "Full 11QC dimer with bound crop graft used ONLY for accessibility; observed glycans and 10A anchor/cut buffers are incomplete glycan coverage. Native-contact positions are exempt from low-total-SASA exclusion alone, never from crop/glycan/excess-exposure exclusions."},
    }
    write_csv(out/"spatial_map.csv", rows)
    manifest["input_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in out.iterdir() if p.is_file()}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"designs":20,"length":88,"hotspots":hotspots,
                      "excluded_target_residues":len(excluded),"output":str(out)}, indent=2))


if __name__ == "__main__":
    prepare()
