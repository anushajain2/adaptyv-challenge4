import csv
import sys
from pathlib import Path

import pytest
import yaml
import json
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from direction_check import release_ratio, histidines, read_structure, enumerate_pairs
from prepare_target import align_map
from read_pilot import normalize, cdr_positions
from Bio.PDB import MMCIFIO, MMCIFParser
from Bio.SeqUtils import seq1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def test_linkage_sign_and_no_shift():
    assert release_ratio([6.5,6.5],[5.5,5.5])>1
    assert release_ratio([5.5,5.5],[6.5,6.5])<1
    assert release_ratio([6.5,6.2],[6.5,6.2])==pytest.approx(1)


def test_linkage_two_site_theoretical_bound():
    assert 1<release_ratio([20,20],[-20,-20])<=10**2.8*(1+1e-10)


def test_missing_histidine_is_error():
    with pytest.raises(ValueError):
        histidines([],"B",[1,2])


def test_alignment_insertion_and_numbering():
    _,mapping=align_map("ACDEFGHIK","ACDWEFGHIK")
    assert mapping[4]==5
    assert mapping[9]==10


def test_target_sequence_indices_and_G142():
    source=read_structure(ROOT/"data/prepared/human_11QC_ectodomain_dimer.pdb")
    old=read_structure(ROOT/"data/prepared/human_6WRV_G142_ectodomain.pdb")
    assert source[0]["A"][142].resname=="GLY"
    assert old[0]["A"][142].resname=="GLY"
    assert set(old[0]["A"][142].child_dict)=={"N","CA","C","O"}
    crop=read_structure(ROOT/"data/prepared/human_apical.pdb")
    assert len(crop[0]["A"])==194
    with (ROOT/"data/prepared/conservation_map.csv").open() as f:
        rows=list(csv.DictReader(f))
    for row in rows:
        if row["in_crop"]=="True":
            i=int(row["boltzgen_index"])
            p=int(row["human_uniprot"])
            assert crop[0]["A"][i].resname==source[0]["A"][p].resname


def test_flat_nanobody_spec_has_real_framework_and_design():
    spec=yaml.safe_load((ROOT/"configs/conserved_326_347_nanobody-anything.yaml").read_text())
    scaffold=spec["entities"][1]["file"]
    assert scaffold["design"]
    assert (ROOT/"configs"/scaffold["path"]).exists()
    target={k:v for e in spec["entities"][0]["file"]["binding_types"] for k,v in e["chain"].items()}
    assert not set(target["binding"].split(",")) & set(target["not_binding"].split(","))


def test_boltzgen_target_has_required_polymer_metadata():
    cif=MMCIF2Dict(str(ROOT/"configs/human_apical.cif"))
    assert len(cif["_entity_poly_seq.num"])==194
    assert cif["_entity_poly_seq.num"]==list(map(str,range(1,195)))
    assert set(cif["_atom_site.label_entity_id"])=={"1"}
    assert cif["_struct_asym.id"]==["A"]


def test_pair_selection_restricts_allowed_positions():
    fixture=read_structure(ROOT/"tests/fixtures/reference_6WRV_not_pilot.pdb")
    pairs=enumerate_pairs(fixture,6,allowed={27,42})
    assert len(pairs)==1
    assert pairs[0]["positions"]==[27,42]


def test_invalid_complex_not_silently_imported(tmp_path):
    with pytest.raises(ValueError):
        normalize(ROOT/"configs/human_apical.cif",tmp_path/"bad.pdb")


def test_import_reference_restores_dimer_and_chain_mapping(tmp_path):
    reference=read_structure(ROOT/"tests/fixtures/reference_6WRV_not_pilot.pdb")
    reference[0]["A"].id="X"
    reference[0]["B"].id="Y"
    io=MMCIFIO()
    io.set_structure(reference)
    io.save(str(tmp_path/"reference.cif"))
    sequence,rms,clashes=normalize(tmp_path/"reference.cif",tmp_path/"normalized.pdb")
    assert len(sequence)==93
    assert rms<2
    assert clashes>=0
    normalized=read_structure(tmp_path/"normalized.pdb")
    assert {c.id for c in normalized[0]}=={"A","B"}
    full=read_structure(tmp_path/"normalized_full_context.pdb")
    target=np.array([a.coord for c in full[0] if c.id!="C" for a in c.get_atoms()])
    binder=np.array([a.coord for a in full[0]["C"].get_atoms()])
    # Independent brute-force count checks the spatial-index implementation.
    brute=sum(int((np.linalg.norm(binder[i:i+32,None]-target[None],axis=2)<2).sum())
              for i in range(0,len(binder),32))
    assert clashes==brute


def test_cdr_mapping_on_upstream_framework():
    scaffold=MMCIFParser(QUIET=True,auth_chains=False).get_structure("scaffold",ROOT/"configs/7eow.cif")
    sequence=seq1("".join(r.resname for r in scaffold[0]["B"] if r.id[0]==" "))
    positions=cdr_positions(sequence)
    assert positions==list(range(26,35))+list(range(52,60))+list(range(98,119))


@pytest.mark.parametrize("states,expected",[
    (["ok","ok"],1), (["ok"],0), (["ok","failed"],0),
])
def test_direction_report_requires_two_successful_models(tmp_path,states,expected):
    # Synthetic unit-test data, never used as experimental or pilot evidence.
    parent="test_parent"
    manifest=tmp_path/"manifest.json"
    manifest.write_text(json.dumps([{
        "design":parent,"sequence":"AAAA","protocol":"protein-anything",
        "full_dimer_clashes":0,
        "metrics":{"design_to_target_iptm":0.8,"bb_rmsd":1.0},
    }]))
    directory=tmp_path/"direction"/parent
    directory.mkdir(parents=True)
    rows=[]
    for seed,status in enumerate(states):
        rows.append({"positions":"1,3","seed":seed,"status":status,
                     "direction":"release" if status=="ok" else "unresolved",
                     "shift_1":-0.8 if status=="ok" else "",
                     "shift_2":-0.9 if status=="ok" else "",
                     "his_only_Kd6_over_Kd7_4":3 if status=="ok" else "",
                     "provisional_candidate":status=="ok"})
    with (directory/"direction_results.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    out=tmp_path/"report"
    subprocess.run([sys.executable,str(ROOT/"scripts/summarize_direction.py"),
                    "--direction",str(directory.parent),"--manifest",str(manifest),
                    "--out",str(out)],check=True,capture_output=True,text=True)
    summary=json.loads((out/"direction_summary.json").read_text())
    assert summary["provisional_candidates"]==expected
    if expected:
        assert "HAHA" in (out/"provisional_switches.fasta").read_text()
    else:
        assert (out/"provisional_switches.fasta").read_text()==""
