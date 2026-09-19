
from __future__ import annotations

import io
import json
from pathlib import Path

from ase.io import write

from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, MlipMethod
from adit.validate_types import ValidationError

SCRIPT_FILE = "run_mlip.py"
SETTINGS_FILE = "mlip_settings.json"
STRUCTURE_FILE = "structure.extxyz"
DEFAULT_COMMAND = "python3"
OPTIMIZERS = ("LBFGS", "FIRE")
MACE_DOC = "https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html"
MACE_LAMMPS_DOC = "https://mace-docs.readthedocs.io/en/latest/guide/lammps.html"
CHGNET_URL = "https://github.com/CederGroupHub/chgnet"

SCRIPT = r'''#!/usr/bin/env python3
"""ADIT が生成: 機械学習ポテンシャル (ASE の calculator) で 1 つの計算を実行する。設定の値は mlip_settings.json にある。
実行する環境に ase と、選んだモデルのパッケージ (mace-torch か chgnet。どちらも PyTorch を使う) が要る。ADIT は要らない。
結果: results.json (エネルギー [eV]、力の最大値 [eV/Å]、応力 [eV/Å³]、振動数 [cm⁻¹])、final.extxyz (最後の構造)、
      trajectory.extxyz (最適化の各ステップ / MD の dump の間隔ごとの構造)、opt.log / md.log。"""
import json
from importlib import metadata

import numpy as np
from ase import units
from ase.io import read, write

with open("mlip_settings.json", encoding="utf-8") as f:
    S = json.load(f)


def version(pkg):
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return "not installed"


VERSIONS = {p: version(p) for p in ("ase", "mace-torch", "chgnet", "torch")}
print("adit-mlip: " + ", ".join(f"{k} {v}" for k, v in VERSIONS.items()), flush=True)


def make_calculator():
    fam, model = S["model_family"], S["model"]
    kw = {}
    if fam in ("mace_mp", "mace_off"):
        from mace.calculators import mace_mp, mace_off
        if model:
            kw["model"] = model
        if S["device"]:
            kw["device"] = S["device"]
        if S["dtype"]:
            kw["default_dtype"] = S["dtype"]
        if fam == "mace_mp":
            return mace_mp(dispersion=S["dispersion"], **kw)
        return mace_off(**kw)
    if fam == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator
        return CHGNetCalculator(use_device=S["device"] or None)
    raise SystemExit(f"unknown model_family: {fam}")


atoms = read("structure.extxyz")
cons = []
if S["fixed_atoms"]:
    from ase.constraints import FixAtoms
    cons.append(FixAtoms(indices=S["fixed_atoms"]))
for k, move in S["fixed_axes"].items():
    if int(k) in S["fixed_atoms"]:
        continue
    mask = [not m for m in move]  # FixCartesian の mask は、固定する成分が True
    if any(mask):
        from ase.constraints import FixCartesian
        cons.append(FixCartesian(int(k), mask=mask))
if cons:
    atoms.set_constraint(cons)
atoms.calc = make_calculator()
results = {"task": S["task"], "model_family": S["model_family"], "model": S["model"] or "(package default)", "versions": VERSIONS}


def summary(a):
    f = a.get_forces()
    out = {"energy_ev": float(a.get_potential_energy()), "fmax_ev_per_ang": float(np.linalg.norm(f, axis=1).max())}
    if a.pbc.all():
        try:
            out["stress_voigt_ev_per_ang3"] = [float(x) for x in a.get_stress(voigt=True)]
        except Exception as ex:  # 応力を出さない calculator もある
            out["stress_note"] = str(ex)
    return out


def dump():
    write("trajectory.extxyz", atoms, append=True)


task = S["task"]
if task == "single_point":
    results.update(summary(atoms))
elif task == "geometry_optimization":
    from ase.optimize import FIRE, LBFGS
    target = atoms
    if S["relax_cell"] != "no":
        from ase.filters import FrechetCellFilter
        target = FrechetCellFilter(atoms, hydrostatic_strain=(S["relax_cell"] == "volume_only"))
    opt = {"LBFGS": LBFGS, "FIRE": FIRE}[S["optimizer"]](target, logfile="opt.log")
    opt.attach(dump, interval=1)
    converged = opt.run(fmax=S["fmax_ev_per_ang"], steps=S["max_steps"])
    results.update(summary(atoms), converged=bool(converged), steps=int(opt.nsteps))
elif task == "molecular_dynamics":
    from ase.md import MDLogger
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    md = S["md"]
    rng = np.random.default_rng(S["seed"])
    dt, T, tau = md["timestep_fs"] * units.fs, md["temperature_k"], md["coupling_time_fs"] * units.fs
    if S["velocities_ang_per_fs"] is not None:
        atoms.set_velocities(np.array(S["velocities_ang_per_fs"]) / units.fs)  # Å/fs → ASE の速度の単位 (Å / ASE の時間)
    else:
        MaxwellBoltzmannDistribution(atoms, temperature_K=T, rng=rng)
        Stationary(atoms)
    ens, th = md["ensemble"], md["thermostat"]
    if ens == "NVE":
        from ase.md.verlet import VelocityVerlet
        dyn = VelocityVerlet(atoms, timestep=dt)
    elif ens == "NVT" and th == "berendsen":
        from ase.md.nvtberendsen import NVTBerendsen
        dyn = NVTBerendsen(atoms, timestep=dt, temperature_K=T, taut=tau)
    elif ens == "NVT" and th == "andersen":
        from ase.md.andersen import Andersen
        dyn = Andersen(atoms, timestep=dt, temperature_K=T, andersen_prob=min(1.0, md["timestep_fs"] / md["coupling_time_fs"]), rng=rng)
    elif ens == "NVT" and th == "langevin":
        from ase.md.langevin import Langevin
        dyn = Langevin(atoms, timestep=dt, temperature_K=T, friction=1.0 / tau, rng=rng)
    elif ens == "NVT" and th == "nose_hoover":
        from ase.md.nose_hoover_chain import NoseHooverChainNVT
        dyn = NoseHooverChainNVT(atoms, timestep=dt, temperature_K=T, tdamp=tau)
    elif ens == "NVT" and th == "csvr":
        from ase.md.bussi import Bussi
        dyn = Bussi(atoms, timestep=dt, temperature_K=T, taut=tau, rng=rng)
    elif ens == "NPT" and th == "nose_hoover":
        from ase.md.nose_hoover_chain import IsotropicMTKNPT
        dyn = IsotropicMTKNPT(atoms, timestep=dt, temperature_K=T, pressure_au=md["pressure_bar"] * units.bar,
                              tdamp=tau, pdamp=md["barostat_time_fs"] * units.fs)
    else:
        raise SystemExit(f"unsupported ensemble / thermostat: {ens} / {th}")
    dyn.attach(MDLogger(dyn, atoms, "md.log", header=True, stress=False, peratom=False), interval=md["dump_interval"])
    dyn.attach(dump, interval=md["dump_interval"])
    dyn.run(md["steps"])
    results.update(summary(atoms), steps=int(md["steps"]))
elif task == "vibrations":
    from ase.vibrations import Vibrations
    moved = [i for i in range(len(atoms)) if i not in set(S["fixed_atoms"])]
    vib = Vibrations(atoms, indices=moved, name="vib")  # 変位の大きさは ASE の既定 (0.01 Å)
    vib.run()
    vib.summary(log="vib_summary.txt")
    freqs = vib.get_frequencies()  # cm^-1。虚数は負の数にして書く
    results.update(summary(atoms), frequencies_cm1=[float(x.real) if abs(x.imag) < 1e-9 else -float(x.imag) for x in freqs])
else:
    raise SystemExit(f"unsupported task: {task}")
write("final.extxyz", atoms)
with open("results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)
print("adit-mlip: done", json.dumps({k: v for k, v in results.items() if k != "frequencies_cm1"}), flush=True)
'''


class MlipGenerator(InputGenerator):
    code = "mlip"
    uses_kpoints = False
    writes_velocities = True

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, MlipMethod):
            return [ValidationError("method.code", L(f"機械学習ポテンシャルの生成器に {m.code!r} の手法が渡されました", f"the MLIP generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        if not m.model_family:
            errs.append(ValidationError("method.model_family", L("機械学習ポテンシャルの種類 (mace_mp / mace_off / chgnet) を選んでください",
                                                                 "choose the machine-learning potential (mace_mp / mace_off / chgnet)")))
        if st.charge != 0 or st.multiplicity != 1:
            errs.append(ValidationError("structure.charge", L("この生成器のモデル (MACE・CHGNet) は全電荷と多重度を入力に取りません。0 と 1 のままにしてください",
                                                              "the models of this generator (MACE, CHGNet) take no total charge or multiplicity; leave them at 0 and 1")))
        if m.model_family == "chgnet":
            if m.model.strip():
                errs.append(ValidationError("method.model", L("CHGNet のモデルの指定の書き方は確かめていないので、空にしてください (パッケージの既定のモデルを使います)",
                                                              "the way to choose a CHGNet model was not verified; leave it empty (the package default is used)")))
            if m.dtype:
                errs.append(ValidationError("method.dtype", L("数値の精度 (default_dtype) は MACE の引数です。CHGNet では空にしてください", "default_dtype is a MACE argument; leave it empty for CHGNet")))
        if m.dispersion and m.model_family != "mace_mp":
            errs.append(ValidationError("method.dispersion", L("D3 の分散補正 (dispersion) は mace_mp の引数です", "the D3 dispersion option is an argument of mace_mp")))
        if m.seed < 0:
            errs.append(ValidationError("method.seed", L("乱数の種は 0 以上の整数にしてください", "the seed must be a non-negative integer")))
        if any(c in m.model for c in "\n\r"):
            errs.append(ValidationError("method.model", L("モデルは 1 行で書いてください", "the model must be a single line")))
        if t.type == "band_structure":
            errs.append(ValidationError("task.type", L("機械学習ポテンシャルにバンド計算はありません", "a machine-learning potential has no band structure")))
        if t.type == "geometry_optimization":
            if t.optimizer not in OPTIMIZERS:
                errs.append(ValidationError("task.optimizer", L(f"ASE にある最適化の方法から選んでください ({' / '.join(OPTIMIZERS)})。{t.optimizer} は ASE にありません",
                                                                f"choose an optimizer available in ASE ({' / '.join(OPTIMIZERS)}); {t.optimizer} is not in ASE")))
            if t.max_steps < 1:
                errs.append(ValidationError("task.max_steps", L("ASE の最適化には最大ステップ数が要ります (1 以上)", "the ASE optimizer needs a maximum number of steps (at least 1)")))
        if t.type == "molecular_dynamics":
            md = t.md
            if md.ensemble == "NPT" and md.thermostat != "nose_hoover":
                errs.append(ValidationError("task.md.thermostat", L("この生成器の NPT は nose_hoover (ASE の IsotropicMTKNPT) だけです (NPTBerendsen は圧縮率が要り、欄がありません)",
                                                                    "NPT in this generator is nose_hoover only (ASE IsotropicMTKNPT); NPTBerendsen needs a compressibility, which has no field")))
            if md.ensemble == "NPT" and not all(st.atoms.pbc):
                errs.append(ValidationError("task.md.ensemble", L("NPT は 3 方向とも周期の系だけです", "NPT needs a system periodic in all three directions")))
        if t.type == "vibrations" and st.fixed_axes:
            errs.append(ValidationError("structure.fixed_axes", L("振動解析 (ASE の Vibrations) は軸ごとの固定を扱いません (原子ごとの固定は動かす原子から外します)",
                                                                  "vibrations (ASE Vibrations) cannot fix individual axes (fixed atoms are left out of the displaced atoms)")))
        if t.relax_cell != "no" and t.type == "geometry_optimization" and not all(st.atoms.pbc):
            errs.append(ValidationError("task.relax_cell", L("格子を動かすのは 3 方向とも周期の系だけです", "the cell can be relaxed only for systems periodic in all three directions")))
        return errs

    def version_probe(self, spec):
        return ("output.log", "adit-mlip: ase")

    def _model_file(self, spec: CalculationSpec) -> Path | None:
        p = Path(spec.method.model.strip()).expanduser() if spec.method.model.strip() else None
        return p if p is not None and p.is_file() else None

    def settings(self, spec: CalculationSpec) -> dict:
        m, t, st = spec.method, spec.task, spec.structure
        mf = self._model_file(spec)
        return {"model_family": m.model_family, "model": mf.name if mf else m.model.strip(), "device": m.device.strip(), "dtype": m.dtype,
                "dispersion": m.dispersion, "seed": m.seed, "task": t.type, "optimizer": t.optimizer, "fmax_ev_per_ang": t.force_tolerance_ev_per_ang,
                "max_steps": t.max_steps, "relax_cell": t.relax_cell, "md": t.md.model_dump(mode="json"),
                "fixed_atoms": sorted(st.fixed_atoms), "fixed_axes": {k: list(v) for k, v in st.fixed_axes.items()},
                "velocities_ang_per_fs": [list(v) for v in st.velocities] if st.velocities is not None else None}

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        buf = io.StringIO()
        write(buf, spec.atoms, format="extxyz")
        return {SCRIPT_FILE: SCRIPT, SETTINGS_FILE: json.dumps(self.settings(spec), indent=2, ensure_ascii=False) + "\n", STRUCTURE_FILE: buf.getvalue()}

    def files_to_copy(self, spec: CalculationSpec, res) -> dict[str, Path]:
        mf = self._model_file(spec)
        if mf is not None and mf.name in (SCRIPT_FILE, SETTINGS_FILE, STRUCTURE_FILE):
            raise GenerationError(L(f"モデルのファイル名 {mf.name} は生成した名前と重なります。名前を変えてください", f"the model file name {mf.name} clashes with a generated file; rename it"))
        return {mf.name: mf} if mf is not None else {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} {SCRIPT_FILE} > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res, copies) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        fam = m.model_family
        files = [
            L("  run_mlip.py        計算を実行する Python のスクリプト (ASE を使う。中身は計算の種類によらず同じ文)",
              "  run_mlip.py        the Python script that runs the calculation (uses ASE; the same text for every task)"),
            L("  mlip_settings.json  このディレクトリの設定の値 (モデル、計算の種類、ステップ数、温度など)。run_mlip.py が読みます",
              "  mlip_settings.json  the settings of this directory (model, task, steps, temperature, ...); read by run_mlip.py"),
            L("  structure.extxyz   構造 (拡張 xyz 形式。セルと周期も入っています。長さの単位は Å)",
              "  structure.extxyz   structure (extended xyz, with the cell and periodicity; lengths in Å)"),
        ]
        if copies:
            files.append(L(f"  {', '.join(sorted(copies))}   モデルのファイル (写したもの)", f"  {', '.join(sorted(copies))}   model file (copied)"))
        pkg = {"mace_mp": "mace-torch", "mace_off": "mace-torch", "chgnet": "chgnet"}.get(fam, "")
        prep = [
            L(f"  実行する環境に ase と {pkg} を入れます (どちらも pip で入ります。{pkg} は PyTorch を一緒に入れるので数 GB になります):",
              f"  Install ase and {pkg} in the environment that runs it (both via pip; {pkg} pulls in PyTorch, several GB):"),
            f"       pip install ase {pkg}",
            L("  ADIT 自身はこれらを使いません (このディレクトリの計算を実行する環境にだけ要ります)。",
              "  ADIT itself does not use them (they are needed only where this directory is run)."),
        ]
        if fam in ("mace_mp", "mace_off"):
            prep += [L(f"  モデル: {m.model or '(mace-torch のバージョンごとの既定。使ったバージョンは output.log の adit-mlip: の行)'}。初めて使うときに mace-torch がモデルをダウンロードします。",
                       f"  Model: {m.model or '(the default of the installed mace-torch version; the version is on the adit-mlip: line of output.log)'}. mace-torch downloads it on first use."),
                     L(f"  モデルごとのライセンス: MACE の文書の表 ({MACE_DOC}) では MACE-MP-0 は MIT、MACE-OFF23・MACE-MATPES-0・MACE-MH は",
                       f"  Licenses per model: the table in the MACE docs ({MACE_DOC}) gives MIT for MACE-MP-0 and the Academic Software License for"),
                     L("  Academic Software License (学術用) です。表にないモデルは、その配布元で確かめてください。",
                       "  MACE-OFF23, MACE-MATPES-0 and MACE-MH; check other models with their distributor.")]
        if fam == "chgnet":
            prep += [L(f"  モデル: chgnet のバージョンごとの既定 (使ったバージョンは output.log の adit-mlip: の行)。chgnet のライセンスは配布元 ({CHGNET_URL}) の LICENSE を見てください。",
                       f"  Model: the default of the installed chgnet version (on the adit-mlip: line of output.log). See the LICENSE at {CHGNET_URL}.")]
        prep += [L("  同じ MACE のモデルを LAMMPS で使うには (ADIT の LAMMPS の生成器の欄で書けます):",
                   "  To use the same MACE model in LAMMPS (the ADIT LAMMPS generator has fields for this):"),
                 L(f"    1. mace の create_lammps_model.py でモデルを LAMMPS 用 (*.model-lammps.pt) に変換します ({MACE_LAMMPS_DOC})",
                   f"    1. convert the model for LAMMPS (*.model-lammps.pt) with mace's create_lammps_model.py ({MACE_LAMMPS_DOC})"),
                 L("    2. LAMMPS の生成器で units = metal、atom_style = atomic、pair_style = mace no_domain_decomposition、",
                   "    2. in the LAMMPS generator set units = metal, atom_style = atomic, pair_style = mace no_domain_decomposition,"),
                 L("       pair_coeff = * * <変換したファイル> <元素を型番号の順に>、写すファイルに変換したファイルを入れます",
                   "       pair_coeff = * * <converted file> <elements in type order>, and add the converted file to the files to copy"),
                 L("    3. LAMMPS は MACE を組み込んで (PKG_ML-MACE と libtorch) ビルドしたものが要ります",
                   "    3. LAMMPS must be built with MACE (PKG_ML-MACE and libtorch)")]
        out = [L("  output.log     スクリプトの画面の出力 (先頭の adit-mlip: の行に、使ったパッケージのバージョン)", "  output.log     screen output of the script (the first adit-mlip: line has the package versions)"),
               L("  results.json   エネルギー [eV]、力の最大値 [eV/Å]、応力 [eV/Å³] (周期系)", "  results.json   energy [eV], maximum force [eV/Å], stress [eV/Å³] (periodic systems)"),
               L("  final.extxyz   最後の構造", "  final.extxyz   final structure")]
        out += {
            "geometry_optimization": [L("  trajectory.extxyz / opt.log   最適化の各ステップの構造と、ASE の最適化の記録", "  trajectory.extxyz / opt.log   structure at each optimization step and the ASE optimizer log")],
            "molecular_dynamics": [L("  trajectory.extxyz / md.log    MD の軌跡 (dump の間隔ごと) と、時刻・エネルギー・温度の記録 (ASE の MDLogger)",
                                     "  trajectory.extxyz / md.log    MD trajectory (every dump interval) and time, energies and temperature (ASE MDLogger)")],
            "vibrations": [L("  vib_summary.txt / vib/   振動数の表 (ASE の Vibrations) と、変位ごとの力。振動数 [cm⁻¹] は results.json にも (虚数は負の数)",
                             "  vib_summary.txt / vib/   frequency table (ASE Vibrations) and forces per displacement; frequencies [cm^-1] also in results.json (imaginary as negative)")],
        }.get(t, [])
        return ReadmeNotes(program="python3", files=files, prepare=prep, outputs=out)


register(MlipGenerator())
