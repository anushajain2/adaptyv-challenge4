"""Prepare canonical target structures, conservation/SASA map and BoltzGen specs."""
import copy
import csv
import json
import shutil
from pathlib import Path

import numpy as np
import yaml
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFIO, MMCIFParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
OUT = ROOT / "data/prepared"
# A contiguous apical-bearing crop avoids joining distant sequence segments.
START, END = 190, 383
DOCUMENTED = {204, 208, 210, 211, 348, 371}


def align_map(human, mouse):
    aligner = Align.PairwiseAligner(
        mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-10, extend_gap_score=-0.5)
    alignment = aligner.align(human, mouse)[0]
    mapping = {int(h) + 1: int(m) + 1 for h, m in zip(*alignment.indices)
               if h >= 0 and m >= 0}
    return alignment, mapping


class ProteinSelect(Select):
    def __init__(self, chains=None, start=0, end=10000):
        self.chains, self.start, self.end = chains, start, end

    def accept_chain(self, chain):
        return self.chains is None or chain.id in self.chains

    def accept_residue(self, residue):
        return residue.id[0] == " " and self.start <= residue.id[1] <= self.end

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def save_pdb(structure, path, select=None):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), select or ProteinSelect())


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def add_polymer_metadata(path, residues):
    """Bio.PDB writes coordinates alone; BoltzGen requires polymer sequence tables."""
    d = MMCIF2Dict(str(path))
    d["_atom_site.label_entity_id"] = ["1"] * len(d["_atom_site.id"])
    sequence = seq1("".join(r.resname for r in residues))
    d.update({
        "_entity.id": ["1"], "_entity.type": ["polymer"],
        "_entity_poly.entity_id": ["1"],
        "_entity_poly.type": ["polypeptide(L)"],
        "_entity_poly.pdbx_strand_id": ["A"],
        "_entity_poly.pdbx_seq_one_letter_code": [sequence],
        "_entity_poly.pdbx_seq_one_letter_code_can": [sequence],
        "_struct_asym.id": ["A"], "_struct_asym.entity_id": ["1"],
        "_entity_poly_seq.entity_id": ["1"] * len(residues),
        "_entity_poly_seq.num": [str(i) for i in range(1, len(residues)+1)],
        "_entity_poly_seq.mon_id": [r.resname for r in residues],
        "_entity_poly_seq.hetero": ["n"] * len(residues),
    })
    writer = MMCIFIO()
    writer.set_dict(d)
    writer.save(str(path))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    human, mouse = [str(SeqIO.read(RAW / f"{p}.fasta", "fasta").seq)
                    for p in ("P02786", "Q62351")]
    alignment, hm = align_map(human, mouse)
    (OUT / "human_mouse_alignment.txt").write_text(str(alignment))
    parser = MMCIFParser(QUIET=True)
    full = parser.get_structure("11QC", RAW / "11QC.cif")
    for chain in full[0]:
        for r in chain:
            if r.id[0] == " ":
                assert seq1(r.resname) == human[r.id[1]-1], (chain.id, r.id)
    assert full[0]["A"][142].resname == "GLY"
    save_pdb(full, OUT / "human_11QC_ectodomain_dimer.pdb")
    save_pdb(full, OUT / "human_11QC_apical_uniprot.pdb",
             ProteinSelect({"A"}, START, END))
    # Sequential author AND label numbering eliminates parser-dependent indices.
    crop = copy.deepcopy(full)
    for c in list(crop[0]):
        if c.id != "A":
            crop[0].detach_child(c.id)
    chain = crop[0]["A"]
    for r in list(chain):
        if r.id[0] != " " or not START <= r.id[1] <= END:
            chain.detach_child(r.id)
    for r in list(chain):
        chain.detach_child(r.id)
        r.id = (" ", r.id[1] - START + 1, " ")
        chain.add(r)
    save_pdb(crop, OUT / "human_apical.pdb")
    io = MMCIFIO()
    io.set_structure(crop)
    io.save(str(OUT / "human_apical.cif"))
    add_polymer_metadata(OUT / "human_apical.cif", list(chain))
    old = parser.get_structure("6WRV", RAW / "6WRV.cif")
    corrections = []
    for c in old[0]:
        if len([r for r in c if r.id[0] == " "]) < 500:
            continue
        r = c[142]
        assert r.resname in {"GLY", "SER"}
        if r.resname == "SER":
            corrections.append({"structure": "6WRV", "chain": c.id,
                                "uniprot": 142, "from": "SER", "to": "GLY"})
            for a in list(r):
                if a.name not in {"N", "CA", "C", "O", "OXT"}:
                    r.detach_child(a.id)
            r.resname = "GLY"
    save_pdb(old, OUT / "human_6WRV_G142_ectodomain.pdb", ProteinSelect({"A", "B"}))
    save_pdb(old, OUT / "human_6WRV_apical_uniprot.pdb", ProteinSelect({"A"}, START, END))
    # Accessibility is evaluated in the full observed dimer, including observed glycans.
    ShrakeRupley(n_points=200).compute(full[0], level="R")
    ShrakeRupley(n_points=200).compute(crop[0], level="R")
    glycans = np.array([a.coord for c in full[0] for r in c
                       if r.id[0] != " " and r.resname != "HOH" for a in r])
    rows = []
    for p, aa in enumerate(human, 1):
        m = hm.get(p)
        r = full[0]["A"].child_dict.get((" ", p, " "))
        glydist = (float(np.min(np.linalg.norm(
            np.array([a.coord for a in r])[:, None] - glycans[None], axis=2)))
            if r is not None and len(glycans) else None)
        in_crop = START <= p <= END
        rows.append({
            "human_uniprot": p, "human_aa": aa, "mouse_uniprot": m,
            "mouse_aa": mouse[m-1] if m else "-", "conserved": bool(m and aa == mouse[m-1]),
            "documented_six": p in DOCUMENTED,
            "published_mouse_position": {204:207,208:211,210:213,211:214,348:350,371:373}.get(p),
            "alignment_caution": "gap placement differs from published structural correspondence"
            if p in {210, 211, 212} else "",
            "recommended_tip": 202 <= p <= 215 or p in {348, 371},
            "in_assay": 121 <= p <= 760, "in_crop": in_crop,
            "boltzgen_index": p-START+1 if in_crop else None,
            "sasa_dimer_A2": round(r.sasa, 2) if r is not None else None,
            "sasa_crop_A2": round(chain[p-START+1].sasa, 2) if in_crop else None,
            "observed_glycan_distance_A": round(glydist, 2) if glydist is not None else None,
            "near_N317_sequon": abs(p-317) <= 5,
            "crop_boundary": in_crop and min(p-START, END-p) < 5,
        })
    write_csv(OUT / "conservation_map.csv", rows)
    write_csv(OUT / "epitope_divergence.csv",
              [r for r in rows if r["recommended_tip"] and not r["conserved"]])
    # Candidate sets are spatially clustered within stated conserved sequence regions.
    sets = []
    for name, lo, hi in [
        ("conserved_326_347", 326, 347),
        ("conserved_230_244", 230, 244),
        ("conserved_298_324_glycan_caution", 298, 324),
    ]:
        candidates = [r for r in rows if lo <= r["human_uniprot"] <= hi
                      and r["conserved"] and r["sasa_dimer_A2"] is not None
                      and r["sasa_dimer_A2"] >= 25 and not r["near_N317_sequon"]
                      and (r["observed_glycan_distance_A"] is None
                           or r["observed_glycan_distance_A"] > 8)]
        candidates.sort(key=lambda r: r["sasa_dimer_A2"], reverse=True)
        if not candidates:
            raise ValueError(f"No exposed residues for {name}")
        anchor = candidates[0]["human_uniprot"]
        selected = [r for r in candidates
                    if np.linalg.norm(full[0]["A"][r["human_uniprot"]]["CA"].coord
                                      - full[0]["A"][anchor]["CA"].coord) <= 12][:3]
        sets.append({"name": name, "human_positions": [r["human_uniprot"] for r in selected],
                     "cross_species": "sequence-conserved; binding not established",
                     "caution": "N317 glycan proximity" if lo == 298 else
                     "screen full-dimer clashes after generation"})
    sets.append({"name": "apical_tip_basic_direction_probe", "human_positions": [205, 208, 371],
                 "cross_species": "divergent: human-first direction probe",
                 "caution": "basic patch is a hypothesis; PROPKA must check sign"})
    for s in sets:
        s["boltzgen_indices"] = [p-START+1 for p in s["human_positions"]]
    (OUT / "hotspot_sets.json").write_text(json.dumps(sets, indent=2) + "\n")
    (OUT / "preparation_audit.json").write_text(json.dumps({
        "primary": "11QC chain A", "crop_uniprot": [START, END],
        "crop_length": END-START+1, "G142_primary": "already GLY; outside crop",
        "corrections": corrections, "missing_primary_assay_residues": [121],
        "6W3H_caution": "Circularly permuted apical fragment; NOT simple native numbering. Not used as design target.",
        "sasa": "Shrake-Rupley, 200 points; observed full dimer and glycans",
        "glycan_limit": "Observed glycans are incomplete; no assurance of glycan clearance.",
        "mouse_model": "No modeled mouse structure generated; conservation map only.",
        "alignment": "Global BLOSUM62, gap open -10, extend -0.5. Sequence gap at V210 differs from published gap at L212; 210-212 excluded from conserved hotspots.",
        "hotspots": "Provisional exposed sets, not experimentally established epitopes.",
    }, indent=2) + "\n")
    configs = ROOT / "configs"
    configs.mkdir(exist_ok=True)
    shutil.copy(OUT / "human_apical.cif", configs / "human_apical.cif")
    scaffold = yaml.safe_load((ROOT / "vendor/boltzgen/example/nanobody_scaffolds/7eow.yaml").read_text())
    shutil.copy(ROOT / "vendor/boltzgen/example/nanobody_scaffolds/7eow.cif", configs / "7eow.cif")
    # Inline scaffold YAML so the upstream wrapper can bundle every path.
    forbidden = [r["boltzgen_index"] for r in rows if r["in_crop"] and
                 (r["crop_boundary"] or not r["conserved"] or r["near_N317_sequon"]
                  or (r["sasa_crop_A2"]-r["sasa_dimer_A2"] > 20))]
    for s in sets:
        for protocol in ("protein-anything", "nanobody-anything"):
            target = {"path": "human_apical.cif", "include": [{"chain": {"id": "A"}}],
                      "binding_types": [{"chain": {
                          "id": "A", "binding": ",".join(map(str, s["boltzgen_indices"])),
                          "not_binding": ",".join(str(p) for p in forbidden
                                                  if p not in s["boltzgen_indices"]),
                      }}]}
            binder = ({"protein": {"id": "B", "sequence": "80..110"}}
                      if protocol == "protein-anything" else {"file": scaffold})
            spec = {"entities": [{"file": target}, binder]}
            (configs / f"{s['name']}_{protocol}.yaml").write_text(yaml.safe_dump(spec, sort_keys=False))
    print(json.dumps({"crop_length": END-START+1, "hotspots": sets}, indent=2))


if __name__ == "__main__":
    main()
