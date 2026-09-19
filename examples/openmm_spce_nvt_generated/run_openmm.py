#!/usr/bin/env python3
"""adit が生成: OpenMM で 1 つの計算を実行する。設定の値は openmm_settings.json にある。
実行する環境に openmm が要る (conda-forge の openmm。ADIT は要らない)。力場・原子型・電荷は
外部のトポロジー (GROMACS の .top か Amber の prmtop) が決める。
結果: results.json (エネルギー [kJ/mol と eV]、温度)、final.pdb (最後の構造)、md.log (MD の記録、CSV)、
      trajectory.dcd (MD の軌跡)、trajectory.extxyz (設定で有効にしたときだけ)。"""
import json

import openmm
from openmm import app, unit

with open("openmm_settings.json", encoding="utf-8") as f:
    S = json.load(f)

print("adit-openmm: openmm " + openmm.version.full_version, flush=True)

NONBONDED = {"NoCutoff": app.NoCutoff, "CutoffNonPeriodic": app.CutoffNonPeriodic, "CutoffPeriodic": app.CutoffPeriodic,
             "Ewald": app.Ewald, "PME": app.PME, "LJPME": app.LJPME}
CONSTRAINTS = {"none": None, "HBonds": app.HBonds, "AllBonds": app.AllBonds, "HAngles": app.HAngles}

if S["input_format"] == "gromacs":
    conf = app.GromacsGroFile(S["coordinates_file"])
    box = conf.getPeriodicBoxVectors()
    top = app.GromacsTopFile(S["topology_file"], periodicBoxVectors=box,
                             includeDir=S["include_dir"] or None, defines=S["defines"] or None)
else:
    top = app.AmberPrmtopFile(S["topology_file"])
    conf = app.AmberInpcrdFile(S["coordinates_file"])
    box = conf.boxVectors
positions = conf.positions

kw = {"nonbondedMethod": NONBONDED[S["nonbonded_method"]], "rigidWater": S["rigid_water"]}
if S["nonbonded_method"] != "NoCutoff":
    kw["nonbondedCutoff"] = S["nonbonded_cutoff_nm"] * unit.nanometer
if CONSTRAINTS[S["constraints"]] is not None:
    kw["constraints"] = CONSTRAINTS[S["constraints"]]
system = top.createSystem(**kw)

task, md = S["task"], S["md"]
temperature = md["temperature_k"] * unit.kelvin
step_size = (md["timestep_fs"] if task == "molecular_dynamics" else 1.0) * unit.femtosecond
collision = 1.0 / (md["coupling_time_fs"] * unit.femtosecond)
ensemble = md["ensemble"] if task == "molecular_dynamics" else "NVE"
thermostat = md["thermostat"]
if ensemble == "NVE" or thermostat == "andersen":
    integrator = openmm.VerletIntegrator(step_size)
    if ensemble != "NVE":
        andersen = openmm.AndersenThermostat(temperature, collision)
        andersen.setRandomNumberSeed(S["seed"])  # 乱数の種を残して、同じ結果を出せるようにする
        system.addForce(andersen)
elif thermostat == "langevin":
    integrator = openmm.LangevinMiddleIntegrator(temperature, collision, step_size)
    integrator.setRandomNumberSeed(S["seed"])
elif thermostat == "nose_hoover":
    integrator = openmm.NoseHooverIntegrator(temperature, collision, step_size)  # 決定論的 (乱数を使わない)
else:
    raise SystemExit(f"unsupported thermostat: {thermostat}")
if S["plumed_file"]:
    # openmm-plumed の PlumedForce に、利用者が書いた PLUMED の入力をそのまま渡す
    from openmmplumed import PlumedForce
    with open(S["plumed_file"], encoding="utf-8") as handle:
        system.addForce(PlumedForce(handle.read()))
if ensemble == "NPT":
    barostat = openmm.MonteCarloBarostat(md["pressure_bar"] * unit.bar, temperature, S["barostat_interval_steps"])
    barostat.setRandomNumberSeed(S["seed"])
    system.addForce(barostat)

platform = openmm.Platform.getPlatformByName(S["platform"]) if S["platform"] else None
simulation = app.Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(positions)
if S["periodic"] and box is not None:
    simulation.context.setPeriodicBoxVectors(*box)

results = {"task": task, "input_format": S["input_format"], "openmm_version": openmm.version.full_version,
           "platform": simulation.context.getPlatform().getName(), "ensemble": ensemble if task == "molecular_dynamics" else None}
KJ_EV = 1.0 / 96.48533212331  # 1 kJ/mol → eV (CODATA 2018)


def record(prefix):
    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    results[prefix + "potential_energy_kj_per_mol"] = pe
    results[prefix + "potential_energy_ev"] = pe * KJ_EV
    if task == "molecular_dynamics":
        results[prefix + "kinetic_energy_kj_per_mol"] = ke
        results[prefix + "total_energy_kj_per_mol"] = pe + ke
        results[prefix + "total_energy_ev"] = (pe + ke) * KJ_EV


class ExtxyzReporter:
    """軌跡を拡張 xyz (テキスト) で書く。報告の約束は古い形の組
    (ステップ数, 座標, 速度, 力, エネルギー, 周期境界で折り返すか) で書く (新しい OpenMM は自分で dict に直す)。
    元素はトポロジーから取る。仮想サイト (element が None) がある系では、黙って書かずに止める。"""

    def __init__(self, path, interval, topology):
        self._path, self._interval = path, interval
        self._symbols = [a.element.symbol if a.element is not None else None for a in topology.atoms()]
        if any(s is None for s in self._symbols):  # 仮想サイト (元素の無い粒子)。黙って falsify しない
            raise SystemExit("adit-openmm: this topology has particles without an element (virtual sites); "
                             "set write_xyz_trajectory to false and use trajectory.dcd")

    def describeNextReport(self, simulation):
        return (self._interval - simulation.currentStep % self._interval, True, False, False, False, None)

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=False)
        vectors = state.getPeriodicBoxVectors()
        cell = " ".join(f"{c.value_in_unit(unit.angstrom):.6f}" for v in vectors for c in v)
        pbc = "T T T" if simulation.system.usesPeriodicBoundaryConditions() else "F F F"
        with open(self._path, "a", encoding="utf-8") as out:
            out.write(f"{len(self._symbols)}\n")
            out.write(f'Lattice="{cell}" Properties=species:S:1:pos:R:3 pbc="{pbc}" step={simulation.currentStep}\n')
            for symbol, xyz in zip(self._symbols, pos):
                x, y, z = (c.value_in_unit(unit.angstrom) for c in xyz)
                out.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")


if task == "single_point":
    record("")
elif task == "geometry_optimization":
    record("initial_")
    simulation.minimizeEnergy(tolerance=S["force_tolerance_kj_per_mol_nm"] * unit.kilojoule_per_mole / unit.nanometer,
                              maxIterations=S["max_steps"])
    record("")
    results["max_iterations"] = S["max_steps"]
elif task == "molecular_dynamics":
    if S["velocities_ang_per_fs"] is not None:
        vectors = [openmm.Vec3(*[c * 100.0 for c in v]) for v in S["velocities_ang_per_fs"]]  # Å/fs → nm/ps
        simulation.context.setVelocities(unit.Quantity(vectors, unit.nanometer / unit.picosecond))
    else:
        simulation.context.setVelocitiesToTemperature(temperature, S["seed"])
    simulation.reporters.append(app.StateDataReporter(
        "md.log", md["dump_interval"], step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, density=True))
    simulation.reporters.append(app.DCDReporter("trajectory.dcd", md["dump_interval"]))
    if S["write_xyz_trajectory"]:
        simulation.reporters.append(ExtxyzReporter("trajectory.extxyz", md["dump_interval"], top.topology))
    simulation.step(md["steps"])
    record("")
    results["steps"] = md["steps"]
    state = simulation.context.getState(getEnergy=True)
    dof = sum(3 for i in range(system.getNumParticles()) if system.getParticleMass(i) > 0 * unit.dalton)
    dof -= system.getNumConstraints()
    if any(isinstance(system.getForce(i), openmm.CMMotionRemover) for i in range(system.getNumForces())):
        dof -= 3
    kt = 2 * state.getKineticEnergy() / (dof * unit.MOLAR_GAS_CONSTANT_R)
    results["final_temperature_k"] = kt.value_in_unit(unit.kelvin)
    results["degrees_of_freedom"] = dof
else:
    raise SystemExit(f"unsupported task: {task}")

final_state = simulation.context.getState(getPositions=True)
final = final_state.getPositions()
if S["periodic"]:  # NPT では箱が変わるので、最後の箱を CRYST1 に書く
    simulation.topology.setPeriodicBoxVectors(final_state.getPeriodicBoxVectors())
with open("final.pdb", "w", encoding="utf-8") as f:
    app.PDBFile.writeFile(simulation.topology, final, f)
with open("results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)
print("adit-openmm: done " + json.dumps(results), flush=True)
