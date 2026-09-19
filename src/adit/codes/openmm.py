
from __future__ import annotations

import json
from pathlib import Path

from adit.codes.amber import _prmtop_natoms, _rst7_natoms
from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.gromacs import _gro_atom_count, gmx_top_dir, scan_includes
from adit.codes.plumed import PLUMED_FILE, output_lines as plumed_outputs, plumed_text, readme_lines as plumed_prepare
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, OpenmmMethod
from adit.validate_types import ValidationError

SCRIPT_FILE = "run_openmm.py"
SETTINGS_FILE = "openmm_settings.json"
TOP_NAME = {"gromacs": "topol.top", "amber": "system.prmtop"}
CONF_NAME = {"gromacs": "conf.gro", "amber": "system.rst7"}
DEFAULT_COMMAND = "python3"
KJ_PER_MOL_NM_PER_EV_ANG = 964.8533212
PERIODIC_METHODS = ("CutoffPeriodic", "Ewald", "PME", "LJPME")
THERMOSTATS = ("langevin", "nose_hoover", "andersen")
DOC = "https://docs.openmm.org/latest/userguide/application/02_running_sims.html"

SCRIPT = r'''#!/usr/bin/env python3
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
'''


class OpenmmGenerator(InputGenerator):

    code = "openmm"
    uses_kpoints = False
    writes_velocities = True

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, OpenmmMethod):
            return [ValidationError("method.code", L(f"OpenMM の生成器に {m.code!r} の手法が渡されました",
                                                     f"the OpenMM generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        if not m.input_format:
            errs.append(ValidationError("method.input_format", L(
                "外部トポロジーの形式 (gromacs か amber) を選んでください",
                "choose the external topology format (gromacs or amber)")))
        top = Path(m.topology_file).expanduser() if m.topology_file.strip() else None
        conf = Path(m.coordinates_file).expanduser() if m.coordinates_file.strip() else None
        if top is None:
            errs.append(ValidationError("method.topology_file", L(
                "トポロジー (GROMACS の .top か Amber の prmtop) を指定してください",
                "give the topology (GROMACS .top or Amber prmtop)")))
        elif not top.is_file():
            errs.append(ValidationError("method.topology_file", L(f"ファイルがありません: {top}", f"file not found: {top}")))
        if conf is None:
            errs.append(ValidationError("method.coordinates_file", L(
                "座標のファイル (GROMACS の .gro か Amber の rst7 / inpcrd) を指定してください",
                "give the coordinate file (GROMACS .gro or Amber rst7 / inpcrd)")))
        elif not conf.is_file():
            errs.append(ValidationError("method.coordinates_file", L(f"ファイルがありません: {conf}", f"file not found: {conf}")))
        if m.input_format == "gromacs":
            if conf is not None and conf.is_file() and conf.suffix.lower() != ".gro":
                errs.append(ValidationError("method.coordinates_file", L(
                    "OpenMM の GROMACS の入口は .gro だけを読みます (GromacsGroFile)",
                    "the OpenMM GROMACS entry point reads .gro only (GromacsGroFile)")))
            if top is not None and top.is_file():
                copies, problems, from_lib = scan_includes(top, "", gmx_top_dir())
                errs += [ValidationError("method.topology_file", p) for p in problems]
                inc = Path(m.include_dir).expanduser() if m.include_dir.strip() else None
                if from_lib and inc is None:
                    errs.append(ValidationError("method.include_dir", L(
                        f"トポロジーが GROMACS の力場 ({', '.join(sorted(set(from_lib)))}) を #include しています。"
                        "その置き場所 (share/gromacs/top) を include_dir に指定してください。OpenMM は GMXLIB を見ません",
                        f"the topology includes GROMACS force-field files ({', '.join(sorted(set(from_lib)))}); "
                        "set include_dir to that folder (share/gromacs/top). OpenMM does not read GMXLIB")))
                elif inc is not None and not inc.is_dir():
                    errs.append(ValidationError("method.include_dir", L(
                        f"ディレクトリがありません: {inc}", f"directory not found: {inc}")))
            if conf is not None and conf.is_file():
                n = _gro_atom_count(conf)
                if n is not None and n != len(st.atoms.symbols):
                    errs.append(ValidationError("structure.atoms", L(
                        f"構造の原子数 ({len(st.atoms.symbols)}) と {conf.name} の原子数 ({n}) が違います。構造にも同じファイルを読み込んでください",
                        f"the structure has {len(st.atoms.symbols)} atoms but {conf.name} has {n}; load the same file as the structure")))
        elif m.input_format == "amber":
            if m.include_dir.strip() or m.defines:
                errs.append(ValidationError("method.include_dir", L(
                    "include_dir と defines は GROMACS のトポロジーだけの設定です",
                    "include_dir and defines apply to GROMACS topologies only")))
            for key, path, reader in (("topology_file", top, _prmtop_natoms), ("coordinates_file", conf, _rst7_natoms)):
                if path is None or not path.is_file():
                    continue
                n = reader(path)
                if n is None:
                    errs.append(ValidationError(f"method.{key}", L(
                        f"{path.name} から原子数を読めません (テキスト形式の prmtop / rst7 を指定してください)",
                        f"cannot read the atom count from {path.name} (give a text prmtop / rst7)")))
                elif n != len(st.atoms.symbols):
                    errs.append(ValidationError("structure.atoms", L(
                        f"構造の原子数 ({len(st.atoms.symbols)}) と {path.name} の原子数 ({n}) が違います",
                        f"the structure has {len(st.atoms.symbols)} atoms but {path.name} has {n}")))
        if not m.nonbonded_method:
            errs.append(ValidationError("method.nonbonded_method", L(
                "非結合相互作用の扱い (NoCutoff / CutoffNonPeriodic / CutoffPeriodic / Ewald / PME / LJPME) を選んでください。OpenMM の既定 (NoCutoff) は使いません",
                "choose the nonbonded method (NoCutoff / CutoffNonPeriodic / CutoffPeriodic / Ewald / PME / LJPME); the OpenMM default (NoCutoff) is not assumed")))
        periodic = bool(m.nonbonded_method in PERIODIC_METHODS)
        if m.nonbonded_method and m.nonbonded_method != "NoCutoff" and m.nonbonded_cutoff_nm <= 0:
            errs.append(ValidationError("method.nonbonded_cutoff_nm", L(
                "カットオフ [nm] は 0 より大きい値にしてください", "the cut-off [nm] must be greater than 0")))
        if m.nonbonded_method == "NoCutoff" and m.nonbonded_cutoff_nm != 0:
            errs.append(ValidationError("method.nonbonded_cutoff_nm", L(
                "NoCutoff ではカットオフを 0 にしてください", "set the cut-off to 0 with NoCutoff")))
        if periodic and not all(st.atoms.pbc):
            errs.append(ValidationError("method.nonbonded_method", L(
                f"{m.nonbonded_method} は周期系の扱いです。構造が 3 方向とも周期でないので選べません",
                f"{m.nonbonded_method} is for periodic systems, but the structure is not periodic in all three directions")))
        if m.nonbonded_method in ("NoCutoff", "CutoffNonPeriodic") and all(st.atoms.pbc) and st.periodic:
            errs.append(ValidationError("method.nonbonded_method", L(
                f"{m.nonbonded_method} は周期境界を使いません。周期系では CutoffPeriodic / Ewald / PME / LJPME から選んでください",
                f"{m.nonbonded_method} ignores periodic boundaries; for a periodic system choose CutoffPeriodic / Ewald / PME / LJPME")))
        if st.charge != 0 or st.multiplicity != 1:
            errs.append(ValidationError("structure.charge", L(
                "OpenMM の力場では電荷はトポロジーが決めます。全電荷と多重度は 0 と 1 のままにしてください",
                "with an OpenMM force field the charges come from the topology; leave the total charge and multiplicity at 0 and 1")))
        if st.fixed_atoms or st.fixed_axes:
            errs.append(ValidationError("structure.fixed_atoms", L(
                "この生成器は固定原子・固定軸を入力に書けません (OpenMM では質量 0 や拘束力で表しますが、拘束との組み合わせを確かめていません)",
                "this generator cannot write fixed atoms or axes (OpenMM expresses them with zero masses or restraint forces, whose interaction with constraints was not verified)")))
        if spec.kpoints is not None:
            errs.append(ValidationError("kpoints", L(
                "OpenMM は k 点を使いません。指定を外してください", "OpenMM does not use k-points; remove this setting")))
        if spec.handoff is not None:
            errs.append(ValidationError("handoff", L(
                "この OpenMM 生成器は続きの計算 (チェックポイント) に対応しません",
                "this OpenMM generator does not support restarts (checkpoints)")))
        if t.type in ("vibrations", "band_structure"):
            errs.append(ValidationError("task.type", L(
                "この OpenMM 生成器は一点計算・エネルギー最小化・分子動力学だけです",
                "this OpenMM generator supports single point, energy minimization and molecular dynamics only")))
        if t.type == "geometry_optimization":
            if t.relax_cell != "no":
                errs.append(ValidationError("task.relax_cell", L(
                    "OpenMM のエネルギー最小化 (LocalEnergyMinimizer) は箱を動かしません",
                    "OpenMM energy minimization (LocalEnergyMinimizer) does not change the box")))
            if t.max_steps < 1:
                errs.append(ValidationError("task.max_steps", L(
                    "最大反復回数は 1 以上にしてください (minimizeEnergy の maxIterations)",
                    "the maximum number of iterations must be at least 1 (maxIterations of minimizeEnergy)")))
        if t.type == "molecular_dynamics":
            md = t.md
            if md.ensemble != "NVE" and md.thermostat not in THERMOSTATS:
                errs.append(ValidationError("task.md.thermostat", L(
                    f"OpenMM にある熱浴は {' / '.join(THERMOSTATS)} です ({md.thermostat} は OpenMM にありません)",
                    f"OpenMM provides {' / '.join(THERMOSTATS)} ({md.thermostat} is not in OpenMM)")))
            if md.ensemble != "NVE" and md.coupling_time_fs <= 0:
                errs.append(ValidationError("task.md.coupling_time_fs", L(
                    "熱浴の時定数は 0 より大きい値にしてください (OpenMM には衝突頻度 = 1 / 時定数 として渡します)",
                    "the thermostat coupling time must be greater than 0 (it is passed to OpenMM as a collision frequency = 1 / coupling time)")))
            if md.ensemble == "NPT":
                if not periodic:
                    errs.append(ValidationError("method.nonbonded_method", L(
                        "NPT (MonteCarloBarostat) は周期系の扱い (CutoffPeriodic / Ewald / PME / LJPME) が要ります",
                        "NPT (MonteCarloBarostat) needs a periodic nonbonded method (CutoffPeriodic / Ewald / PME / LJPME)")))
                if m.barostat_interval_steps < 1:
                    errs.append(ValidationError("method.barostat_interval_steps", L(
                        "圧力浴を試す間隔は 1 ステップ以上にしてください (MonteCarloBarostat の frequency)",
                        "the barostat interval must be at least 1 step (frequency of MonteCarloBarostat)")))
            if md.dump_interval < 1:
                errs.append(ValidationError("task.md.dump_interval", L(
                    "書き出しの間隔は 1 以上にしてください", "the dump interval must be at least 1")))
        if m.seed < 0:
            errs.append(ValidationError("method.seed", L("乱数の種は 0 以上の整数にしてください", "the seed must be a non-negative integer")))
        if any(c in m.platform for c in " \n\r"):
            errs.append(ValidationError("method.platform", L(
                "プラットフォームの名前は 1 語で書いてください (CPU / CUDA / OpenCL / Reference)",
                "the platform name must be a single word (CPU / CUDA / OpenCL / Reference)")))
        return errs

    def version_probe(self, spec):
        return ("output.log", "adit-openmm: openmm")

    def settings(self, spec: CalculationSpec) -> dict:
        m, t, st = spec.method, spec.task, spec.structure
        return {"input_format": m.input_format,
                "topology_file": TOP_NAME.get(m.input_format, "topology"),
                "coordinates_file": CONF_NAME.get(m.input_format, "coordinates"),
                "include_dir": m.include_dir.strip(), "defines": dict(m.defines),
                "nonbonded_method": m.nonbonded_method, "nonbonded_cutoff_nm": m.nonbonded_cutoff_nm,
                "constraints": m.constraints, "rigid_water": m.rigid_water,
                "periodic": bool(m.nonbonded_method in PERIODIC_METHODS),
                "platform": m.platform.strip(), "seed": m.seed, "write_xyz_trajectory": m.write_xyz_trajectory,
                "plumed_file": PLUMED_FILE if spec.plumed is not None else "",
                "barostat_interval_steps": m.barostat_interval_steps,
                "task": t.type, "max_steps": t.max_steps,
                "force_tolerance_kj_per_mol_nm": t.force_tolerance_ev_per_ang * KJ_PER_MOL_NM_PER_EV_ANG,
                "md": t.md.model_dump(mode="json"),
                "velocities_ang_per_fs": [list(v) for v in st.velocities] if st.velocities is not None else None}

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        out = {SCRIPT_FILE: SCRIPT, SETTINGS_FILE: json.dumps(self.settings(spec), indent=2, ensure_ascii=False) + "\n"}
        text = plumed_text(spec)
        if text is not None:
            out[PLUMED_FILE] = text
        return out

    def files_to_copy(self, spec: CalculationSpec, res) -> dict[str, Path]:
        m = spec.method
        if not m.input_format:
            raise GenerationError(L("外部トポロジーの形式 (gromacs か amber) を選んでください",
                                    "choose the external topology format (gromacs or amber)"))
        top = Path(m.topology_file).expanduser()
        out = {TOP_NAME[m.input_format]: top, CONF_NAME[m.input_format]: Path(m.coordinates_file).expanduser()}
        if m.input_format == "gromacs":
            copies, problems, _ = scan_includes(top, "", gmx_top_dir())
            if problems:
                raise GenerationError("; ".join(problems))
            out.update(copies)
        return out

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, DEFAULT_COMMAND).format(
            mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} {SCRIPT_FILE} > output.log 2>&1"

    # ---- README ----
    def readme_notes(self, spec: CalculationSpec, res, copies: dict[str, Path]) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        top_name, conf_name = TOP_NAME[m.input_format] if m.input_format else "topology", CONF_NAME.get(m.input_format, "coordinates")
        files = [
            L("  run_openmm.py       計算を実行する Python のスクリプト (openmm を使う。中身は計算の種類によらず同じ文)",
              "  run_openmm.py       the Python script that runs the calculation (uses openmm; the same text for every task)"),
            L("  openmm_settings.json  このディレクトリの設定の値 (非結合の扱い、カットオフ、拘束、ステップ数、温度など)",
              "  openmm_settings.json  the settings of this directory (nonbonded method, cut-off, constraints, steps, temperature, ...)"),
            L(f"  {top_name} / {conf_name}   利用者が用意したトポロジーと座標 (写したもの)",
              f"  {top_name} / {conf_name}   topology and coordinates supplied by the user (copied)"),
        ]
        extra = sorted(n for n in copies if n not in (top_name, conf_name))
        if extra:
            files.append(L(f"  {', '.join(extra)}   トポロジーが #include で読むファイル",
                           f"  {', '.join(extra)}   files read by #include from the topology"))
        prepare = [
            L("  実行する環境に openmm を入れます (conda-forge。ADIT 自身は openmm を使いません):",
              "  Install openmm where this directory is run (conda-forge; ADIT itself does not use openmm):"),
            "       conda install -c conda-forge openmm",
            L("  力場・原子型・電荷はトポロジーが決めます。ADIT は作りません。トポロジーと座標の原子の並びが同じかを確かめてください。",
              "  The topology sets the force field, atom types and charges; ADIT does not build them. Check that topology and coordinates use the same atom order."),
        ]
        if m.input_format == "gromacs" and m.include_dir.strip():
            prepare.append(L(f"  GROMACS の力場は {m.include_dir} から読みます (openmm_settings.json の include_dir)。別の環境で実行するときは、その環境のパスに直してください。",
                             f"  GROMACS force-field files are read from {m.include_dir} (include_dir in openmm_settings.json); change it when running elsewhere."))
        if spec.task.type == "molecular_dynamics":
            md = spec.task.md
            if md.ensemble != "NVE":
                prepare.append(L(f"  熱浴の時定数 {md.coupling_time_fs:g} fs は、衝突頻度 1 / ({md.coupling_time_fs:g} fs) = {1.0 / md.coupling_time_fs:.4g} fs⁻¹ として OpenMM に渡します。",
                                 f"  The coupling time {md.coupling_time_fs:g} fs is passed to OpenMM as a collision frequency of "
                                 f"1/({md.coupling_time_fs:g} fs) = {1.0 / md.coupling_time_fs:.4g} fs^-1."))
            if md.ensemble == "NPT":
                prepare.append(L(f"  圧力浴 (MonteCarloBarostat) は {m.barostat_interval_steps} ステップごとに体積を試します。共通の設定の圧力浴の時定数 ({md.barostat_time_fs:g} fs) は、OpenMM のこの間隔に機械的に換算できないので入力に使っていません。",
                                 f"  The barostat (MonteCarloBarostat) attempts a volume move every {m.barostat_interval_steps} steps. The common barostat time ({md.barostat_time_fs:g} fs) has no mechanical conversion to this interval and is not used."))
        outputs = [L("  output.log     スクリプトの画面の出力 (先頭の adit-openmm: の行に openmm のバージョン)",
                     "  output.log     screen output of the script (the first adit-openmm: line has the openmm version)"),
                   L("  results.json   エネルギー [kJ/mol と eV]。MD では最後の温度も",
                     "  results.json   energies [kJ/mol and eV]; for MD also the final temperature"),
                   L("  final.pdb      最後の構造", "  final.pdb      the final structure")]
        if t == "molecular_dynamics":
            if m.write_xyz_trajectory:
                outputs.append(L("  trajectory.extxyz   軌跡をテキストの拡張 xyz でも書いたもの (ADIT の解析の RDF・MSD はこちらを読みます。DCD より大きくなります)",
                                 "  trajectory.extxyz   the trajectory also as a text extended xyz (ADIT's RDF and MSD read this; larger than the DCD)"))
            outputs.append(L("  md.log / trajectory.dcd   時刻・エネルギー・温度・体積・密度の記録 (CSV。OpenMM の StateDataReporter) と軌跡 (DCD)",
                             "  md.log / trajectory.dcd   time, energies, temperature, volume and density (CSV from the OpenMM StateDataReporter) and the trajectory (DCD)"))
        if spec.plumed is not None:
            files.append(L(f"  {PLUMED_FILE}   PLUMED の入力 (利用者が書いたもの)", f"  {PLUMED_FILE}   the PLUMED input (written by you)"))
            prepare += plumed_prepare(spec)
            outputs += plumed_outputs(spec)
        return ReadmeNotes(program="python3", files=files, prepare=prepare, outputs=outputs)


register(OpenmmGenerator())


# The reference the current OpenMM user guide asks for ("Referencing OpenMM")
CITATIONS = (
    Citation("openmm_eastman2024", r"""@article{openmm_eastman2024,
  author  = {Eastman, Peter and Galvelis, Raimondas and Pel{\'a}ez, Ra{\'u}l P. and Abreu, Charlles R. A. and Farr, Stephen E. and Gallicchio, Emilio and Gorenko, Anton and Henry, Michael M. and Hu, Frank and Huang, Jing and Kr{\"a}mer, Andreas and Michel, Julien and Mitchell, Joshua A. and Pande, Vijay S. and Rodrigues, Jo{\~a}o P. G. L. M. and Rodriguez-Guerra, Jaime and Simmonett, Andrew C. and Singh, Sukrit and Swails, Jason and Turner, Philip and Wang, Yuanqing and Zhang, Ivy and Chodera, John D. and De Fabritiis, Gianni and Markland, Thomas E.},
  title   = {{OpenMM} 8: Molecular Dynamics Simulation with Machine Learning Potentials},
  journal = {The Journal of Physical Chemistry B},
  volume  = {128},
  number  = {1},
  pages   = {109--116},
  year    = {2024},
  doi     = {10.1021/acs.jpcb.3c06662}
}""", doi="10.1021/acs.jpcb.3c06662", source="https://docs.openmm.org/latest/userguide/introduction.html"),
)
