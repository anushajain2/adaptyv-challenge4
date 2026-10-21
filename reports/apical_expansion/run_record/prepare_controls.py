"""Prepare published-complex controls, explicit residue maps, and decoys.

These are calibration inputs, never newly designed competition candidates.
"""
import copy
import hashlib
import json
import random
from pathlib import Path

import numpy as np
from Bio.PDB import Chain, Model, Structure, Superimposer, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1, seq3
from scipy.spatial import cKDTree

from direction_check import read_structure, full_context
from prepare_target import ROOT, save_pdb

OUT = ROOT / "data/controls"
REPORT = ROOT / "reports/controls"


def residues(chain):
    return [r for r in chain if r.id[0] == " "]


def sequence(chain):
    return seq1("".join(r.resname for r in residues(chain)))


def make_complex(target, binder):
    s = Structure.Structure("control")
    m = Model.Model(0)
    s.add(m)
    for name, source in [("A", target), ("B", binder)]:
        c = Chain.Chain(name)
        for i, r in enumerate(source, 1):
            # Avoid recursively copying the parent structure.
            r = r.copy()
            r.id = (" ", i, " ")
            c.add(r)
        m.add(c)
    return s


def fit_mapped(source, native_map, full):
    """Map target A onto full receptor A, applying the same rigid transform to B."""
    s = copy.deepcopy(source)
    pairs = [(i, n) for i, n in native_map.items()
             if "CA" in s[0]["A"][i] and (" ", n, " ") in full[0]["A"]]
    sup = Superimposer()
    sup.set_atoms([full[0]["A"][n]["CA"] for i, n in pairs],
                  [s[0]["A"][i]["CA"] for i, n in pairs])
    sup.apply(list(s.get_atoms()))
    return s, float(sup.rms)


def geometry(source, native_map, full):
    fitted, rms = fit_mapped(source, native_map, full)
    bxyz = np.array([a.coord for a in fitted[0]["B"].get_atoms() if a.element != "H"])
    native_set = set(native_map.values())
    regions = {
        "modeled_target": [a for a in fitted[0]["A"].get_atoms() if a.element != "H"],
        "omitted_same_subunit": [a for r in full[0]["A"]
                                if r.id[0] == " " and r.id[1] not in native_set
                                for a in r if a.element != "H"],
        "partner_subunit": [a for c in full[0] if c.id != "A"
                            for r in c if r.id[0] == " "
                            for a in r if a.element != "H"],
        "observed_glycans": [a for c in full[0] for r in c
                            if r.id[0] != " " and r.resname != "HOH"
                            for a in r if a.element != "H"],
    }
    result = {"target_fit_rmsd_A": rms, "mapped_target_residues": len(native_map)}
    for name, atoms in regions.items():
        result[f"{name}_clash_pairs_lt2A"] = int(cKDTree(
            np.array([a.coord for a in atoms])).query_ball_point(
                bxyz, np.nextafter(2., 0.), return_length=True).sum()) if atoms else 0
    # These are contacts in experimental coordinates, not inferred hotspots.
    tree = cKDTree(bxyz)
    result["native_target_contacts_lt5A"] = [
        native_map[r.id[1]] for r in residues(fitted[0]["A"])
        if tree.query(np.array([a.coord for a in r if a.element != "H"]))[0].min() < 5
    ]
    return result


def save_input_cif(s, path, target_sequence, binder_sequence):
    """Retain explicit full polymer sequences, including unresolved residues."""
    io = MMCIFIO()
    io.set_structure(s)
    io.save(str(path))
    d = MMCIF2Dict(str(path))
    chains = ["A", "B"]
    seqs = [target_sequence, binder_sequence]
    d["_atom_site.label_entity_id"] = [
        str(chains.index(c)+1) for c in d["_atom_site.auth_asym_id"]]
    # Author indices were explicitly normalized to full-sequence indices.
    d["_atom_site.label_seq_id"] = list(d["_atom_site.auth_seq_id"])
    d["_entity.id"] = ["1", "2"]
    d["_entity.type"] = ["polymer", "polymer"]
    d["_entity_poly.entity_id"] = ["1", "2"]
    d["_entity_poly.type"] = ["polypeptide(L)"] * 2
    d["_entity_poly.pdbx_strand_id"] = chains
    for key in ["pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can"]:
        d["_entity_poly."+key] = seqs
    d["_struct_asym.id"] = chains
    d["_struct_asym.entity_id"] = ["1", "2"]
    d["_entity_poly_seq.entity_id"] = [str(i+1) for i, seq in enumerate(seqs) for _ in seq]
    d["_entity_poly_seq.num"] = [str(j+1) for seq in seqs for j in range(len(seq))]
    d["_entity_poly_seq.mon_id"] = [seq3(aa).upper() for seq in seqs for aa in seq]
    d["_entity_poly_seq.hetero"] = ["n"] * sum(map(len, seqs))
    io.set_dict(d)
    io.save(str(path))


def main():
    OUT.mkdir(exist_ok=True)
    REPORT.mkdir(exist_ok=True)
    raw = read_structure(ROOT / "data/raw/6WRV.cif")
    full = read_structure(ROOT / "data/raw/11QC.cif")
    fixture = read_structure(ROOT / "tests/fixtures/reference_6WRV_not_pilot.pdb")
    # Identify which experimentally observed binder copy contacts receptor A.
    target_xyz = np.array([a.coord for r in raw[0]["A"] if r.id[0] == " " for a in r])
    tree = cKDTree(target_xyz)
    choices = {}
    for c in raw[0]:
        if len(residues(c)) == 93:
            choices[c.id] = int(tree.query_ball_point(
                np.array([a.coord for r in residues(c) for a in r]), 5,
                return_length=True).sum())
    binder_id = max(choices, key=choices.get)
    xtal = make_complex([raw[0]["A"][i] for i in range(190, 384)],
                        residues(raw[0][binder_id]))
    assert sequence(xtal[0]["B"]) == sequence(fixture[0]["B"])
    native_map = {i: i + 189 for i in range(1, 195)}
    fitted, _ = fit_mapped(xtal, native_map, full)
    current = make_complex([full[0]["A"][i] for i in range(190, 384)],
                           residues(fitted[0]["B"]))
    shuffled = copy.deepcopy(current)
    letters = list(sequence(shuffled[0]["B"]))
    random.Random(20261020).shuffle(letters)
    for r, aa in zip(shuffled[0]["B"], letters):
        r.resname = seq3(aa).upper()
        for a in list(r):
            if a.name not in {"N", "CA", "C", "O"}:
                r.detach_child(a.id)
    scaffold = read_structure(ROOT / "configs/7eow.cif")
    # caplacizumab is a VWF binder; use the same observed single-domain framework
    # as the nanobody pilot, without redesigning its CDRs.
    unrelated_chain = next(c for c in scaffold[0] if 110 <= len(residues(c)) <= 150)
    unrelated = make_complex(residues(current[0]["A"]), residues(unrelated_chain))
    # Binder coordinates are excluded from all target templates below.
    for a in unrelated[0]["B"].get_atoms():
        a.coord = a.coord + np.array([100., 100., 100.])
    atv_raw = read_structure(ROOT / "data/raw/6W3H.cif")
    atv_target_residues = [r for r in residues(atv_raw[0]["D"]) if r.id[1] <= 169]
    atv = make_complex(atv_target_residues, residues(atv_raw[0]["A"]))
    # Author D binds author A. PDB label A/C correspond to author A/D.
    # Restore gaps in sequence indices; leave the two unobserved residues unmodeled.
    d = MMCIF2Dict(str(ROOT / "data/raw/6W3H.cif"))
    label_by_auth = {int(a): int(l) for a, l, c, group in zip(
        d["_atom_site.auth_seq_id"], d["_atom_site.label_seq_id"],
        d["_atom_site.auth_asym_id"], d["_atom_site.group_PDB"])
        if c == "D" and group == "ATOM"}
    saved = list(atv[0]["A"])
    for r in saved:
        atv[0]["A"].detach_child(r.id)
    for r, original in zip(saved, atv_target_residues):
        r.id = (" ", label_by_auth[original.id[1]] - 5, " ")
        atv[0]["A"].add(r)
    # Use full canonical sequence, since observed 11QC starts at native residue 122.
    from Bio import SeqIO
    human = str(SeqIO.read(ROOT / "data/raw/P02786.fasta", "fasta").seq)
    atv_seq = d["_entity_poly.pdbx_seq_one_letter_code_can"][1].replace("\n", "")[5:149]
    split = atv_seq.index("DSAQNS")
    first, second = atv_seq[:split], atv_seq[split:]
    start1, start2 = human.index(first) + 1, human.index(second) + 1
    atv_map = {i+1: (start1+i if i < split else start2+i-split)
               for i in range(len(atv_seq))}
    observed_atv_map = {i: n for i, n in atv_map.items() if (" ", i, " ") in atv[0]["A"]}
    assert all(seq1(atv[0]["A"][i].resname) == human[n-1] for i, n in observed_atv_map.items())
    atv_fitted, _ = fit_mapped(atv, observed_atv_map, full)
    atv_current = make_complex(residues(current[0]["A"]), residues(atv_fitted[0]["B"]))
    cases = [
        ("3ds18_current", current, native_map, "positive", "11QC native 190–383 crop"),
        ("3ds18_crystal", xtal, native_map, "positive", "6WRV native 190–383 crop"),
        ("3ds18_scrambled", shuffled, native_map, "negative_hypothesis",
         "Composition-preserving shuffle; not an experimentally established nonbinder"),
        ("caplacizumab_current", unrelated, native_map, "negative_hypothesis",
         "Unrelated folded VWF-binding nanobody from 7EOW; not a tested TfR1 nonbinder"),
        ("atv_permuted", atv, atv_map, "positive_construct",
         "6W3H author receptor D and Fc A; target construct positions 6–149 (shared observed span of C/D), with two unresolved residues retained in sequence; terminal tags and displaced D-only tail omitted"),
        ("atv_current", atv_current, native_map, "positive_context_transfer",
         "Observed Fc A mapped to 11QC native crop; construct-transfer comparison"),
    ]
    records = []
    for name, s, mapping, kind, note in cases:
        path = OUT / f"{name}.pdb"
        save_pdb(s, path)
        target_seq = atv_seq if name == "atv_permuted" else sequence(s[0]["A"])
        binder_seq = sequence(s[0]["B"])
        save_input_cif(s, path.with_suffix(".cif"), target_seq, binder_seq)
        mask = np.array([False] * len(target_seq) + [True] * len(binder_seq))
        np.savez_compressed(path.with_suffix(".npz"), design_mask=mask)
        records.append({
            "id": name, "kind": kind, "note": note,
            "target_sequence": target_seq, "binder_sequence": binder_seq,
            "target_native_map": mapping, "target_chain": "A", "binder_chain": "B",
            "pdb": str(path.relative_to(ROOT)),
            "cif": str(path.with_suffix(".cif").relative_to(ROOT)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "geometry": geometry(s, {i:n for i,n in mapping.items()
                                    if (" ",i," ") in s[0]["A"]}, full)
                        if kind.startswith("positive") else None,
        })
    # Round-trip the existing exact pipeline, plus numbering/chain-independent mapping.
    raw_out = ROOT / "outputs/controls_geometry"
    raw_out.mkdir(exist_ok=True)
    fit, clashes = full_context(xtal, raw_out / "3ds18_full_context.pdb")
    manifest = {
        "scope": "Published controls and hypothesized negatives, not competition candidates",
        "3ds18_source_receptor_chain": "A", "3ds18_source_binder_chain": binder_id,
        "3ds18_chain_pair_contact_counts": choices,
        "legacy_full_context_roundtrip": {"rmsd_A": fit, "clashes": clashes},
        "atv_source_chains": {"target": "D", "binder": "A"},
        "atv_mapping_segments": [[1, split, start1, start1+split-1],
                                 [split+1, len(atv_seq), start2, start2+len(second)-1]],
        "refolding": {"target_templates": True, "binder_templates": False,
                      "msa": "empty, matching pinned pilot loader",
                      "seeds": [17, 29], "diffusion_samples_per_seed": 5,
                      "recycling_steps": 3, "sampling_steps": 200},
        "cases": records,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    (REPORT / "preparation.json").write_text(json.dumps(manifest, indent=2))
    for r in records:
        print(r["id"], len(r["target_sequence"]), len(r["binder_sequence"]), r["geometry"])


if __name__ == "__main__":
    main()
