"""Histidine-pair triage. Structural predictions, never experimental binding claims.

Input: a two-chain complex with target A and binder B, fully modeled backbone.
Outputs: rebuilt/relaxed variants, free/complex structures, raw PROPKA tables,
per-site shifts, and the two-site independent proton-linkage estimate.
"""
import argparse
import copy
import itertools
import json
import logging
import math
import os
import time
import multiprocessing
import random
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

from prepare_target import ROOT, ProteinSelect, save_pdb, write_csv

HIS_ATOMS = {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"}
os.environ.setdefault("OPENMM_CPU_THREADS", "2")

class TargetFitError(ValueError):
    def __init__(self,rms):
        self.rms=float(rms)
        super().__init__(f"Target changed too far for full-context fit: {rms:.2f} A")


def read_structure(path):
    parser = MMCIFParser(QUIET=True) if Path(path).suffix == ".cif" else PDBParser(QUIET=True)
    return parser.get_structure("complex", str(path))


def log_ka_factor(ph, free, bound):
    """Independent sites: Ka(pH)/Ka(deprotonated reference)."""
    return sum(math.log10(1 + 10**(b-ph)) - math.log10(1 + 10**(f-ph))
               for f, b in zip(free, bound))


def release_ratio(free, bound):
    """Kd(6.0)/Kd(7.4); >1 indicates weaker binding at acidic pH."""
    return 10**(log_ka_factor(7.4, free, bound) - log_ka_factor(6.0, free, bound))


def enumerate_pairs(structure, max_pairs=6, allowed=None):
    complex_model = copy.deepcopy(structure[0])
    ShrakeRupley(n_points=100).compute(complex_model, level="R")
    free = copy.deepcopy(structure[0])
    free.detach_child("A")
    ShrakeRupley(n_points=100).compute(free, level="R")
    target = [a for a in complex_model["A"].get_atoms() if a.element != "H"]
    xyz = np.array([a.coord for a in target])
    candidates = []
    for r in complex_model["B"]:
        p = r.id[1]
        if r.id[0] != " " or r.resname in {"GLY", "PRO", "CYS", "HIS"} or "CA" not in r:
            continue
        if allowed is not None and p not in allowed:
            continue
        delta_sasa = free["B"][r.id].sasa - r.sasa
        distances = np.linalg.norm(np.array([a.coord for a in r if a.element != "H"])[:,None]-xyz[None], axis=2)
        minimum = float(distances.min())
        # Periphery heuristic: contacts, modest buried area, some exposed free area.
        if not (minimum <= 5.0 and 2 <= delta_sasa <= 65 and free["B"][r.id].sasa >= 25):
            continue
        nearby = [target[i] for i in np.where(distances.min(axis=0) < 5)[0]]
        basic = sum(a.name in {"NZ", "NH1", "NH2", "NE"} and
                    a.parent.resname in {"LYS", "ARG"} for a in nearby)
        acceptor = sum(a.element == "O" for a in nearby)
        candidates.append({"position": p, "parent": seq1(r.resname),
                           "buried_sasa_A2": round(delta_sasa,2),
                           "min_target_distance_A": round(minimum,2),
                           "nearby_basic_atoms": basic, "nearby_oxygen_atoms": acceptor,
                           "heuristic": 3*basic-acceptor})
    pairs = []
    for a,b in itertools.combinations(candidates, 2):
        d = float(np.linalg.norm(complex_model["B"][a["position"]]["CA"].coord -
                                 complex_model["B"][b["position"]]["CA"].coord))
        if 5 <= d <= 16:
            pairs.append({"positions":[a["position"],b["position"]],
                          "sites":[a,b], "ca_separation_A":round(d,2),
                          "heuristic":a["heuristic"]+b["heuristic"]})
    return sorted(pairs,key=lambda r:r["heuristic"],reverse=True)[:max_pairs]


def rebuild_relax(source, out, positions, seed):
    """Build actual HIS side chains; restrain target and binder backbone."""
    from pdbfixer import PDBFixer
    import openmm as mm
    from openmm import app, unit
    random.seed(seed)
    np.random.seed(seed)
    fixer = PDBFixer(filename=str(source))
    mutations = []
    for r in fixer.topology.residues():
        if r.chain.id == "B" and int(r.id) in positions:
            mutations.append(f"{r.name}-{r.id}-HIS")
    if len(mutations) != len(positions):
        raise ValueError("Mutation positions missing/ambiguous")
    fixer.applyMutations(mutations,"B")
    print(f"  rebuilding {mutations}",flush=True)
    fixer.findMissingResidues()
    # Never invent unobserved loops; source must have continuous backbone.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=seed)
    # Place newly built terminal carboxyl oxygens in the peptide plane.
    # Bad initial O/OXT proximity can be interpreted as an O--O bond by PROPKA.
    coordinates=np.array(fixer.positions.value_in_unit(unit.nanometer))
    for residue in fixer.topology.residues():
        atoms={a.name:a.index for a in residue.atoms()}
        if {"CA","C","O","OXT"}<=set(atoms):
            center=coordinates[atoms["C"]]
            u=coordinates[atoms["CA"]]-center
            v=coordinates[atoms["O"]]-center
            direction=-(u/np.linalg.norm(u)+v/np.linalg.norm(v))
            coordinates[atoms["OXT"]]=center+0.125*direction/np.linalg.norm(direction)
    fixer.positions=coordinates*unit.nanometer
    print("  heavy atoms built; adding hydrogens",flush=True)
    fixer.addMissingHydrogens(7.4)
    print("  relaxing complex",flush=True)
    forcefield = app.ForceField("amber14-all.xml","implicit/gbn2.xml")
    # No dynamics are integrated. Avoid holonomic constraints in minimization:
    # constraint-restraint cycling can exceed maxIterations on rebuilt side chains.
    system = forcefield.createSystem(
        fixer.topology,nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2*unit.nanometer,constraints=None)
    restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k",1000.0)  # kJ/mol/nm^2
    for name in ("x0","y0","z0"):
        restraint.addPerParticleParameter(name)
    coords = fixer.positions.value_in_unit(unit.nanometer)
    for a in fixer.topology.atoms():
        if a.residue.chain.id=="A":
            system.setParticleMass(a.index,0)
        if a.element.symbol != "H" and (a.residue.chain.id=="A" or a.name in {"N","CA","C","O"}):
            restraint.addParticle(a.index,coords[a.index])
    system.addForce(restraint)
    integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
    simulation = app.Simulation(fixer.topology,system,integrator,
                                mm.Platform.getPlatformByName("CPU"),{"Threads":"2"})
    simulation.context.setPositions(fixer.positions)
    class BoundedMinimization(mm.MinimizationReporter):
        def __init__(self):
            super().__init__()
            self.started=time.monotonic()
            self.timed_out=False

        def report(self,iteration,x,grad,args):
            self.timed_out=time.monotonic()-self.started>120
            return self.timed_out

    reporter=BoundedMinimization()
    initial_energy=float(simulation.context.getState(getEnergy=True).getPotentialEnergy()
                         .value_in_unit(unit.kilojoules_per_mole))
    simulation.minimizeEnergy(maxIterations=500,reporter=reporter)
    if reporter.timed_out:
        raise RuntimeError("Relaxation exceeded 120-second diagnostic limit; variant not scored")
    state = simulation.context.getState(getPositions=True,getEnergy=True)
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
    if not math.isfinite(energy):
        raise ValueError("Non-finite minimized energy")
    out.with_suffix(".relaxation.json").write_text(json.dumps({
        "initial_energy_kJ_mol":initial_energy,"final_energy_kJ_mol":energy,
        "method":"Amber14/GBn2, 12 A cutoff, fixed target, restrained binder backbone",
        "max_iterations":500,"seed":seed,
    },indent=2))
    with out.open("w") as f:
        app.PDBFile.writeFile(fixer.topology,state.getPositions(),f,keepIds=True)
    modeled = read_structure(out)
    for p in positions:
        r=modeled[0]["B"][p]
        if r.resname != "HIS" or not HIS_ATOMS <= set(r.child_dict):
            raise ValueError(f"Incomplete histidine at {p}")
    # Strip hydrogens so PROPKA uses one consistent protonation treatment.
    save_pdb(modeled,out)
    return energy


def _relax_worker(source,out,positions,seed,queue):
    try:
        queue.put(("ok",rebuild_relax(source,out,positions,seed)))
    except Exception as e:
        queue.put(("failed",f"{type(e).__name__}: {e}"))


def bounded_relax(source,out,positions,seed):
    """Hard wall-clock guard also covers force evaluations between reporter calls."""
    context=multiprocessing.get_context("spawn")
    queue=context.Queue()
    process=context.Process(target=_relax_worker,args=(source,out,positions,seed,queue))
    process.start()
    process.join(180)
    if process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
            process.join()
        raise RuntimeError("Structure rebuilding/relaxation exceeded 180 seconds; not scored")
    if process.exitcode!=0:
        raise RuntimeError(f"Relaxation process exited {process.exitcode}")
    status,result=queue.get(timeout=5)
    if status!="ok":
        raise RuntimeError(result)
    return result


def full_context(cropped, out):
    """Restore full dimer context using target C-alpha fit; binder becomes chain C."""
    full=read_structure(ROOT/"data/prepared/human_11QC_ectodomain_dimer.pdb")
    mobile=cropped[0]["A"]
    fixed=full[0]["A"]
    residues=[r for r in mobile if r.id[0]==" "]
    if len(residues)!=194:
        raise ValueError("Expected 194-residue prepared target")
    if seq1("".join(r.resname for r in residues)) != seq1("".join(fixed[i].resname for i in range(190,384))):
        raise ValueError("Target sequence does not match prepared crop")
    sup=Superimposer()
    sup.set_atoms([fixed[i]["CA"] for i in range(190,384)],[r["CA"] for r in residues])
    if sup.rms > 2:
        raise TargetFitError(sup.rms)
    binder=copy.deepcopy(cropped[0]["B"])
    sup.apply(list(binder.get_atoms()))
    binder.id="C"
    full[0].add(binder)
    # Preserve the modeled interface when restoring receptor context. Replacing it
    # with the reference conformation can manufacture clashes and spurious pKas.
    mobile_copy=copy.deepcopy(mobile)
    sup.apply(list(mobile_copy.get_atoms()))
    for native,r in zip(range(190,384),[r for r in mobile_copy if r.id[0]==" "]):
        fixed.detach_child((" ",native," "))
        # The chain is already copied. Move each residue out of that copy;
        # deepcopy(residue) would recursively copy its entire parent structure.
        mobile_copy.detach_child(r.id)
        if "OXT" in r:
            r.detach_child("OXT")  # crop terminus becomes internal in the full receptor
        r.id=(" ",native," ")
        fixed.add(r)
    fixed.child_list.sort(key=lambda r:r.id[1])
    # Count all interface clashes, including the modeled crop.
    omitted=np.array([a.coord for c in full[0] if c.id!="C"
                      for a in c.get_atoms() if a.element!="H"])
    bxyz=np.array([a.coord for a in binder.get_atoms() if a.element!="H"])
    # Count atom pairs without allocating a binder x receptor x 3 distance array.
    clash_count=int(cKDTree(omitted).query_ball_point(
        bxyz,np.nextafter(2.0,0.0),return_length=True).sum())
    save_pdb(full,out)
    return sup.rms,clash_count


def pka_table(path):
    from propka.run import single
    messages=[]
    class Capture(logging.Handler):
        def emit(self,record):
            messages.append(record.getMessage())
    handler=Capture(logging.WARNING)
    logger=logging.getLogger("propka")
    logger.addHandler(handler)
    try:
        molecule=single(str(path),write_pka=False)
    finally:
        logger.removeHandler(handler)
    path.with_suffix(".propka_warnings.json").write_text(json.dumps(messages,indent=2))
    groups=molecule.conformations["AVR"].groups
    table=[]
    for group in groups:
        atom=group.atom
        if group.titratable:
            table.append({"chain":atom.chain_id,"residue":atom.res_num,
                          "resname":atom.res_name,"group":group.type,
                          "pka":float(group.pka_value)})
    path.with_suffix(".pka.json").write_text(json.dumps(table,indent=2))
    if any("Missing atoms or failed protonation" in msg for msg in messages):
        raise ValueError("PROPKA reported incomplete/protonation-failed groups; raw tables retained, variant not scored")
    return table


def histidines(table,chain,positions):
    result=[]
    for p in positions:
        hits=[r["pka"] for r in table if r["chain"]==chain and r["residue"]==p
              and r["group"]=="HIS"]
        if len(hits)!=1 or not math.isfinite(hits[0]):
            raise ValueError(f"Missing or duplicate PROPKA histidine: {chain}{p}")
        result.append(hits[0])
    return result


def scan(path,out,max_pairs=6,seeds=(11,29),allowed=None,explicit_pair=None,restore=True):
    out.mkdir(parents=True,exist_ok=True)
    source=read_structure(path)
    chains={c.id for c in source[0]}
    if chains!={"A","B"}:
        raise ValueError(f"Normalize target/binder chains to A/B first, found {chains}")
    save_pdb(source,out/"parent_complex.pdb")
    pairs=([{"positions":list(explicit_pair),"sites":[],"heuristic":None}]
           if explicit_pair else enumerate_pairs(source,max_pairs,allowed))
    (out/"enumerated_pairs.json").write_text(json.dumps(pairs,indent=2))
    rows=[]
    for pair in pairs:
        positions=pair["positions"]
        for seed in seeds:
            name=f"H{positions[0]}_H{positions[1]}_seed{seed}"
            folder=out/name
            folder.mkdir(exist_ok=True)
            row={"variant":name,"positions":",".join(map(str,positions)),"seed":seed,
                 "status":"failed","error":"","full_dimer_clashes":None,
                 "pka_target_context":"intact_194_residue_apical_crop",
                 "free_pka_1":None,"bound_pka_1":None,"shift_1":None,
                 "free_pka_2":None,"bound_pka_2":None,"shift_2":None,
                 "his_only_Kd6_over_Kd7_4":None,"direction":"unresolved",
                 "provisional_candidate":False}
            try:
                bounded_relax(out/"parent_complex.pdb",folder/"complex.pdb",positions,seed)
                modeled=read_structure(folder/"complex.pdb")
                save_pdb(modeled,folder/"free_binder.pdb",ProteinSelect({"B"}))
                bound_path=folder/"complex.pdb"
                if restore:
                    rms,clashes=full_context(modeled,folder/"full_complex.pdb")
                    row["full_dimer_clashes"]=clashes
                # A fitted full-context model is useful for collision checks, but
                # grafting a refolded crop can break bonds at its junctions. Use
                # the intact modeled crop for pKa rather than inventing bonds.
                free=histidines(pka_table(folder/"free_binder.pdb"),"B",positions)
                bound=histidines(pka_table(bound_path),"B",positions)
                shifts=[b-f for f,b in zip(free,bound)]
                ratio=release_ratio(free,bound)
                for i,(f,b,d) in enumerate(zip(free,bound,shifts),1):
                    row.update({f"free_pka_{i}":f,f"bound_pka_{i}":b,f"shift_{i}":d})
                row.update(status="ok",his_only_Kd6_over_Kd7_4=ratio,
                           direction="release" if ratio>1.1 else "trapping" if ratio<1/1.1 else "neutral",
                           provisional_candidate=all(d<=-0.5 for d in shifts)
                           and any(6<=f<=7 for f in free) and ratio>=2
                           and row["full_dimer_clashes"]==0)
            except Exception as e:
                row["error"]=f"{type(e).__name__}: {e}"
                logging.exception("Variant failed: %s",name)
            rows.append(row)
            print(name,row["status"],row["direction"],flush=True)
            write_csv(out/"direction_results.csv",rows)
    summaries=[]
    for pair in pairs:
        label=",".join(map(str,pair["positions"]))
        group=[r for r in rows if r["positions"]==label]
        summaries.append({"positions":label,"successful_models":sum(r["status"]=="ok" for r in group),
                          "all_models_provisional":len(group)>=2 and
                          all(r["provisional_candidate"] for r in group)})
    (out/"summary.json").write_text(json.dumps({
        "pairs":summaries,"tested_models":len(rows),
        "interpretation":"His-only independent-site fixed-conformation triage; not measured affinity or proof of switching.",
        "limits":["PROPKA uses trimmed target; omitted receptor and glycans can affect pKa",
                  "Artificial target termini retained; boundary proximity must be reviewed",
                  "Full-dimer graft is for clash checks only, not a covalently valid pKa model",
                  "No unbound conformational ensemble",
                  "No full proton linkage across all target sites", "Model seeds are not exhaustive rotamer sampling",
                  "Absolute affinity loss at pH 7.4 is not estimated", "No pair passing is a valid result"],
    },indent=2))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("complex",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--max-pairs",type=int,default=6)
    p.add_argument("--seeds",default="11,29")
    p.add_argument("--allowed",help="Comma-separated binder positions; restrict nanobodies to CDRs")
    p.add_argument("--pair",help="Explicit pair for testing, e.g. 12,33")
    p.add_argument("--no-full-context",action="store_true",help="Smoke tests only; cannot qualify candidates")
    args=p.parse_args()
    rows=scan(args.complex,args.out,args.max_pairs,
              tuple(map(int,args.seeds.split(","))),
              set(map(int,args.allowed.split(","))) if args.allowed else None,
              tuple(map(int,args.pair.split(","))) if args.pair else None,
              not args.no_full_context)
    if rows and all(r["status"]=="failed" for r in rows):
        raise SystemExit("All variants failed; inspect direction_results.csv")


if __name__=="__main__":
    main()
