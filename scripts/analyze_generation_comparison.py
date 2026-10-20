"""Apply the frozen exploratory screen to all samples; never select by geometry."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.PDB import Superimposer
from Bio.Align import PairwiseAligner

from analyze_controls import contact_pairs, pose_metrics
from audit_control_samples import sample_structures
from prepare_controls import ROOT, geometry, read_structure, sequence
from prepare_target import write_csv
from read_pilot import KD
from score_matched_controls import run_ipsae, interface_pae


def design_pose_metrics(predicted, reference):
    """IF sidechains are zero-coordinate placeholders, not a designed interface.

    Compare binder backbone contacts to target heavy atoms in BOTH structures.
    Full predicted sidechains remain available to hotspot, PAE and clash checks.
    """
    structures = [copy.deepcopy(predicted),copy.deepcopy(reference)]
    for s in structures:
        for residue in s[0]["B"]:
            for atom in list(residue):
                if atom.name not in {"N","CA","C","O"}:
                    residue.detach_child(atom.id)
    return pose_metrics(*structures)


def screen_sample(row, plan):
    """Return every failure; missing/invalid geometry must never count as clear."""
    pose, g, contacts = row["pose"], row["geometry"], row["contacts"]
    checks = {
        "interface_confidence": row["metrics"]["design_to_target_iptm"] >= plan["min_design_to_target_iptm"],
        "binder_fold": pose["binder_CA_fold_rmsd_A"] <= plan["max_binder_CA_self_fit_rmsd_A"],
        "interface_pose": pose["binder_interface_CA_rmsd_after_target_fit_A"] is not None and
            pose["binder_interface_CA_rmsd_after_target_fit_A"] <= plan["max_interface_CA_rmsd_after_target_fit_A"],
        "design_contact_recall": pose["native_contact_recall"] is not None and
            pose["native_contact_recall"] >= plan["min_design_contact_recall"],
        "hotspot_contacts": len(contacts["hotspots_contacted_crop"]) >= plan["min_hotspots_contacted"],
        "intended_region": contacts["fraction_pairs_in_reference_region"] >= plan["min_contact_pairs_in_reference_region_fraction"],
        "context_fit": g["target_fit_rmsd_A"] <= plan["max_full_context_fit_rmsd_A"],
        "clash_clearance": all(g[k] <= plan["max_clash_pairs_lt2A"] for k in (
            "modeled_target_clash_pairs_lt2A", "omitted_same_subunit_clash_pairs_lt2A",
            "partner_subunit_clash_pairs_lt2A", "observed_glycans_clash_pairs_lt2A")),
    }
    return {"pass": all(checks.values()), "checks": checks,
            "failed": [k for k,v in checks.items() if not v]}


def contacts_summary(s, manifest, conserved):
    pairs = contact_pairs(s)
    positions = {a for a,b in pairs}
    region = set(manifest["reference_contact_region_crop"])
    backbone = []
    for a in s[0]["A"]:
        for b in s[0]["B"]:
            for x,y in (("N","O"),("O","N")):
                if x in a and y in b:
                    distance = float(np.linalg.norm(a[x].coord-b[y].coord))
                    if distance < 3.5:
                        backbone.append({"target_native":a.id[1]+189,"target_atom":x,
                                         "binder_position":b.id[1],"binder_atom":y,"distance_A":distance})
    return {"pairs":len(pairs),"target_native_positions":sorted(p+189 for p in positions),
            "hotspots_contacted_crop":sorted(positions & set(manifest["hotspots_crop"])),
            "fraction_pairs_in_reference_region":sum(a in region for a,b in pairs)/len(pairs) if pairs else 0,
            "fraction_pairs_at_conserved_target_positions":sum(conserved[a+189] for a,b in pairs)/len(pairs) if pairs else 0,
            "backbone_NO_contacts_lt3_5A":backbone,
            "backbone_contact_note":"Distance-only observations, not hydrogen-bond or beta-sheet assignments"}


def aligned_identity(a,b):
    aligner = PairwiseAligner(mode="global",match_score=2,mismatch_score=-1,
                              open_gap_score=-4,extend_gap_score=-.5)
    alignment = aligner.align(a,b)[0]
    return sum(a[i] == b[j] for i,j in zip(*alignment.indices) if i>=0 and j>=0)/max(len(a),len(b))


def families(names, references, sequences):
    """Conservative grouping to avoid counting close variants as independent hits."""
    groups, comparisons = [], []
    related = {}
    for i,a in enumerate(names):
        for b in names[i+1:]:
            ca = [r["CA"] for r in references[a][0]["B"]]
            cb = [r["CA"] for r in references[b][0]["B"]]
            sup = Superimposer()
            sup.set_atoms(ca,cb)
            identity = aligned_identity(sequences[a],sequences[b])
            related[frozenset((a,b))] = identity >= .3 or sup.rms <= 3
            comparisons.append({"a":a,"b":b,"sequence_identity":identity,
                                "binder_CA_self_fit_rmsd_A":float(sup.rms)})
    for name in names:
        for group in groups:
            if all(related[frozenset((name,other))] for other in group):
                group.append(name)
                break
        else: groups.append([name])
    return groups, comparisons


def analyze(run, output):
    output.mkdir(parents=True,exist_ok=True)
    manifest = json.loads((run/"inputs/manifest.json").read_text())
    cases = json.loads((run/"cases.json").read_text())
    feature_audit = []
    for arm in ("unconstrained","sheet"):
        paths = sorted((run/arm/"intermediate_designs").glob("*.npz"))
        assert len(paths) == 10
        for path in paths:
            with np.load(path,allow_pickle=False) as data:
                hotspots = (np.flatnonzero(data["binding_type"][:194] == 1)+1).tolist()
                sheet = (np.flatnonzero(data["ss_type"][194:] == 3)+1).tolist()
                assert hotspots == manifest["hotspots_crop"]
                assert sheet == (list(range(33,38)) if arm=="sheet" else [])
                assert data["design_mask"][194:].all() and not data["design_mask"][:194].any()
            feature_audit.append({"source":str(path.relative_to(run)),"hotspots_crop":hotspots,
                                  "sheet_binder_positions":sheet,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
    full = read_structure(ROOT/"data/raw/11QC.cif")
    conserved = {int(r["native_position"]):r["conserved"] == "True"
                 for r in csv.DictReader((run/"inputs/spatial_map.csv").open())}
    references = {c["id"]:read_structure(run/"references"/(c["id"]+".pdb")) for c in cases}
    known = {name:sequence(read_structure(ROOT/path)[0]["B"]) for name,path in (
        ("2DS25.5","data/controls_matched/2ds25_5_6wrw.pdb"),
        ("3DS18","data/controls_matched/3ds18_6wrv.pdb"))}
    rows = []
    for seed in manifest["refolding"]["seeds"]:
        folder = run/f"seed-{seed}"
        preflight = json.loads((folder/"preflight.json").read_text())
        assert len(preflight) == len(cases)
        assert all(p["binder_template_excluded"] and p["msa_depth"] == 1 for p in preflight)
        for case in cases:
            name = case["id"]
            reference = references[name]
            with np.load(folder/"fold_out_npz"/f"{name}.npz",allow_pickle=False) as data:
                paes, plddts = data["pae"].copy(),data["plddt"].copy()
                assert np.all(data["asym_id"][0,:194] == 0)
                assert np.all(data["asym_id"][0,194:] == 1)
                assert data["token_pad_mask"].all()
                assert paes.shape == (5,len(case["binder_sequence"])+194,len(case["binder_sequence"])+194)
            for index,s,metrics,selected in sample_structures(folder,name,restore_sidechains=True):
                assert sequence(s[0]["A"]) == manifest["target_sequence"]
                assert sequence(s[0]["B"]) == case["binder_sequence"]
                mean_pae = float((paes[index,:194,194:].mean()+paes[index,194:,:194].mean())/2)
                assert abs(mean_pae-metrics["interaction_pae"]) < .03
                dest = ROOT/"outputs/generation_scoring"/run.name/f"seed-{seed}"/name/str(index)
                scores = run_ipsae(s,paes[index],plddts[index],dest)
                scores.update(interface_pae(s,paes[index]))
                scores["binder_mean_plddt"] = float(plddts[index,194:].mean())
                row = {"id":name,"arm":case["arm"],"kind":case["kind"],"seed":seed,
                       "sample":index,"confidence_selected":selected,"metrics":metrics,"scores":scores,
                       "geometry":geometry(s,{i:i+189 for i in range(1,195)},full),
                       "pose":(design_pose_metrics(s,reference) if case["kind"]=="novel_design"
                               else pose_metrics(s,reference) if case["kind"]=="positive" else None),
                       "pose_reference_kind":("generated_backbone" if case["kind"]=="novel_design"
                                              else "experimental_complex" if case["kind"]=="positive"
                                              else "no_experimental_pose"),
                       "contacts":contacts_summary(s,manifest,conserved),
                       "raw_scoring_output":str(dest.relative_to(ROOT))}
                if case["kind"] == "novel_design":
                    row["screen"] = screen_sample(row,manifest["screen"])
                rows.append(row)
            print(f"Scored seed {seed}: {name}",flush=True)
    assert len(rows) == len(cases)*len(manifest["refolding"]["seeds"])*5 == 230
    selected = [r for r in rows if r["confidence_selected"]]
    assert len(selected) == 46
    control_checks = []
    for r in selected:
        if r["arm"] != "control": continue
        if r["kind"] == "positive":
            passed = (r["metrics"]["design_to_target_iptm"] >= .5 and
                      r["pose"]["binder_CA_fold_rmsd_A"] <= 2 and
                      r["pose"]["binder_interface_CA_rmsd_after_target_fit_A"] <= 2 and
                      r["pose"]["native_contact_recall"] >= .5)
        else: passed = r["metrics"]["design_to_target_iptm"] < .5
        control_checks.append({"id":r["id"],"seed":r["seed"],"pass":bool(passed)})
    controls_pass = len(control_checks) == 6 and all(c["pass"] for c in control_checks)
    summaries = []
    for case in cases:
        if case["kind"] != "novel_design": continue
        name = case["id"]
        subset = [r for r in selected if r["id"] == name]
        assert len(subset) == 2
        seq = case["binder_sequence"]
        structural_pass = all(r["screen"]["pass"] for r in subset)
        summary = {"id":name,"arm":case["arm"],"sequence":seq,"length":len(seq),
                   "mean_kyte_doolittle":sum(KD[a] for a in seq)/len(seq),
                   "identity_to_published_controls":{k:aligned_identity(seq,v) for k,v in known.items()},
                   "both_seeds_structural_pass":structural_pass,
                   "advances":controls_pass and structural_pass and seq not in known.values(),
                   "selected_models":[{"seed":r["seed"],"sample":r["sample"],
                       "interface_iptm":r["metrics"]["design_to_target_iptm"],
                       "binder_CA_fold_rmsd_A":r["pose"]["binder_CA_fold_rmsd_A"],
                       "interface_CA_pose_rmsd_A":r["pose"]["binder_interface_CA_rmsd_after_target_fit_A"],
                       "design_contact_recall":r["pose"]["native_contact_recall"],
                       "contact_pae_A":r["scores"]["contact_pae_both_directions_A"],
                       "ipsae_min":r["scores"]["ipsae_min"],"ipsae_max":r["scores"]["ipsae_max"],
                       "failed_checks":r["screen"]["failed"]} for r in subset],
                   "passing_samples_out_of_10":sum(r["screen"]["pass"] for r in rows if r["id"]==name)}
        summaries.append(summary)
    advances = [s["id"] for s in summaries if s["advances"]]
    groups, comparisons = families(advances,references,{c["id"]:c["binder_sequence"] for c in cases})
    arms = {arm:{"designs":10,"advance_count":sum(s["advances"] for s in summaries if s["arm"]==arm),
                 "selected_interface_iptm_range":[min(r["metrics"]["design_to_target_iptm"] for r in selected if r["arm"]==arm),
                                                  max(r["metrics"]["design_to_target_iptm"] for r in selected if r["arm"]==arm)]}
            for arm in ("unconstrained","sheet")}
    result = {"run":run.name,"manifest_sha256":hashlib.sha256((run/"inputs/manifest.json").read_bytes()).hexdigest(),
              "actual_generated_feature_audit":feature_audit,
              "screen":manifest["screen"],"controls_pass":controls_pass,"control_checks":control_checks,
              "arms":arms,"candidate_summaries":summaries,"advancing_families":groups,
              "family_comparisons":comparisons,"rows":rows,
              "interpretation":"Prospective structural triage only. In novel rows, native_contact_recall means recall of design binder-backbone/target-heavy-atom contacts, not experimental contacts. Inverse-folded reference sidechains are placeholders and excluded from this metric; predicted full atoms are used for hotspot, PAE and clash checks. No binding, affinity, mouse cross-reactivity or pH switch is established."}
    (output/"results.json").write_text(json.dumps(result,indent=2))
    write_csv(output/"designs.csv",[{"id":s["id"],"arm":s["arm"],"sequence":s["sequence"],
        "advances":s["advances"],"mean_kyte_doolittle":s["mean_kyte_doolittle"],
        "passing_samples_out_of_10":s["passing_samples_out_of_10"],
        "selected_models_json":json.dumps(s["selected_models"])} for s in summaries])
    for filename,subset in (("all_designs.fasta",summaries),("advancing_designs.fasta",[s for s in summaries if s["advances"]])):
        (output/filename).write_text("".join(f">{s['id']}\n{s['sequence']}\n" for s in subset))
    print(json.dumps({"controls_pass":controls_pass,"arms":arms,"families":groups},indent=2))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run",type=Path)
    p.add_argument("--output",type=Path,default=ROOT/"reports/generation_comparison")
    args = p.parse_args()
    analyze(args.run,args.output)
