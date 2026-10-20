"""Add explicit spatial crop-boundary and glycosylation-site checks to target maps."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from prepare_controls import ROOT, read_structure
from prepare_target import write_csv


def main():
    full = read_structure(ROOT / "data/raw/11QC.cif")
    chain = full[0]["A"]
    boundaries = np.array([chain[190]["N"].coord, chain[383]["C"].coord])
    anchors = []
    for c in full[0]:
        for p in (251, 317, 727):
            if (" ", p, " ") in c and c[p].resname == "ASN" and "ND2" in c[p]:
                anchors.append(c[p]["ND2"].coord)
    assert len(anchors) == 6
    partner = cKDTree(np.array([a.coord for c in full[0] if c.id != "A"
                               for r in c if r.id[0] == " " for a in r if a.element != "H"]))
    prepared = json.loads((ROOT / "data/controls/manifest.json").read_text())
    contacts = set(next(c for c in prepared["cases"] if c["id"] == "3ds18_crystal")
                   ["geometry"]["native_target_contacts_lt5A"])
    rows = []
    with (ROOT / "data/prepared/conservation_map.csv").open() as f:
        for original in csv.DictReader(f):
            n = int(original["human_uniprot"])
            if not 190 <= n <= 383:
                continue
            xyz = np.array([a.coord for a in chain[n] if a.element != "H"])
            boundary = float(np.linalg.norm(xyz[:, None]-boundaries[None], axis=2).min())
            glycan_anchor = float(np.linalg.norm(xyz[:, None]-np.array(anchors)[None], axis=2).min())
            reasons = []
            if boundary < 10:
                reasons.append("within_10A_of_crop_terminus")
            if glycan_anchor < 10:
                reasons.append("within_10A_of_glycan_anchor")
            if float(original["sasa_dimer_A2"]) < 25:
                reasons.append("low_full_dimer_SASA")
            if float(original["sasa_crop_A2"])-float(original["sasa_dimer_A2"]) > 20:
                reasons.append("crop_exposes_buried_surface")
            original.update({
                "distance_to_crop_cut_atoms_A": round(boundary, 3),
                "distance_to_glycan_anchor_ND2_A": round(glycan_anchor, 3),
                "nearest_partner_protein_atom_A": round(float(partner.query(xyz)[0].min()), 3),
                "published_3ds18_contact": n in contacts,
                "spatial_screen_pass": not reasons, "exclusion_reasons": ";".join(reasons)})
            rows.append(original)
    out = ROOT / "reports/controls"
    write_csv(out / "spatial_target_map.csv", rows)
    lookup = {int(r["human_uniprot"]): r for r in rows}
    sets = json.loads((ROOT / "data/prepared/hotspot_sets.json").read_text())
    result = {
        "geometry_source": "Full observed 11QC dimer; existing dimer SASA and observed-glycan distances retained",
        "glycan_anchors": ["A251", "A317", "A727", "B251", "B317", "B727"],
        "limits": "10 A cut/anchor buffers and 25 A2 SASA are heuristics, not validated glycan envelopes or binding guarantees.",
        "original_hotspot_sets": [
            {"name": s["name"], "positions": s["human_positions"],
             "passing_positions": [n for n in s["human_positions"] if lookup[n]["spatial_screen_pass"]],
             "excluded": {n: lookup[n]["exclusion_reasons"] for n in s["human_positions"]
                          if not lookup[n]["spatial_screen_pass"]}} for s in sets],
        "published_interface_passing_spatial_screen": [
            n for n in sorted(contacts) if lookup[n]["spatial_screen_pass"]],
        "total_residues_passing_spatial_screen": sum(r["spatial_screen_pass"] for r in rows),
    }
    (out / "spatial_target_audit.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
