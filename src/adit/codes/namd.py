"""NAMD 3 NVE input with user-supplied PSF, PDB and CHARMM parameter files."""
# Official syntax and sample configurations:
# https://www.ks.uiuc.edu/Research/namd/3.0/ug/node9.html
# https://www.ks.uiuc.edu/Research/namd/3.0/ug/node12.html
# https://www.ks.uiuc.edu/Research/namd/3.0/ug/node91.html

from __future__ import annotations

import math
from pathlib import Path

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, NamdMethod
from adit.validate_types import ValidationError


def _psf_natoms(path: Path) -> int | None:
    try:
        with path.open(encoding="ascii", errors="replace") as stream:
            for line in stream:
                if "!NATOM" in line.upper():
                    return int(line.split()[0])
    except (OSError, ValueError, IndexError):
        pass
    return None


def _pdb_natoms(path: Path) -> int | None:
    try:
        with path.open(encoding="ascii", errors="replace") as stream:
            return sum(line.startswith(("ATOM  ", "HETATM")) for line in stream)
    except OSError:
        return None


def _pdb_positions(path: Path) -> list[tuple[float, float, float]] | None:
    try:
        positions: list[tuple[float, float, float]] = []
        with path.open(encoding="ascii", errors="replace") as stream:
            for line in stream:
                if line.startswith(("ATOM  ", "HETATM")):
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    if not all(math.isfinite(x) for x in xyz):
                        return None
                    positions.append(xyz)
        return positions
    except (OSError, ValueError):
        return None


class NamdGenerator(InputGenerator):
    code = "namd"
    uses_kpoints = False
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m, st, task = spec.method, spec.structure, spec.task
        if not isinstance(m, NamdMethod):
            return [ValidationError("method.code", L("NAMD の条件ではありません", "the settings are not for NAMD"))]
        errors: list[ValidationError] = []
        if task.type != "molecular_dynamics" or task.md.ensemble != "NVE":
            errors.append(ValidationError("task.type", L(
                "この NAMD 生成器は NVE 分子動力学だけに対応します",
                "this NAMD generator supports NVE molecular dynamics only")))
        if st.fixed_atoms or st.fixed_axes:
            errors.append(ValidationError("structure.fixed_atoms", L(
                "この NAMD 生成器は固定原子・固定軸を入力に書けません",
                "this NAMD generator cannot write fixed atoms or axes")))
        if st.periodic and not all(st.atoms.pbc):
            errors.append(ValidationError("structure.atoms", L(
                "この NAMD 生成器は 3 次元周期系か非周期系だけに対応します",
                "this NAMD generator supports only fully periodic or nonperiodic systems")))
        if spec.kpoints is not None:
            errors.append(ValidationError("kpoints", L(
                "NAMD は k 点を使用しません。指定を外してください",
                "NAMD does not use k-points; remove this setting")))
        if spec.handoff is not None:
            errors.append(ValidationError("handoff", L(
                "この NAMD 生成器は続きの計算に対応しません", "this NAMD generator does not support restarts")))
        for key, suffix in (("structure_file", ".psf"), ("coordinates_file", ".pdb")):
            value = getattr(m, key)
            path = Path(value).expanduser() if value.strip() else None
            if path is None or not path.is_file() or path.suffix.lower() != suffix:
                errors.append(ValidationError(f"method.{key}", L(
                    f"{key} に既存の {suffix} ファイルを指定してください",
                    f"set {key} to an existing {suffix} file")))
        psf = Path(m.structure_file).expanduser() if m.structure_file.strip() else None
        pdb = Path(m.coordinates_file).expanduser() if m.coordinates_file.strip() else None
        for key, path, reader in (("structure_file", psf, _psf_natoms), ("coordinates_file", pdb, _pdb_natoms)):
            if path is not None and path.is_file():
                n = reader(path)
                if n is None or n == 0:
                    errors.append(ValidationError(f"method.{key}", L(
                        f"{path.name} から原子数を読めません", f"cannot read the atom count from {path.name}")))
                elif n != len(st.atoms.symbols):
                    errors.append(ValidationError("structure.atoms", L(
                        f"Spec は {len(st.atoms.symbols)} 原子、{path.name} は {n} 原子です",
                        f"the spec has {len(st.atoms.symbols)} atoms but {path.name} has {n}")))
        if pdb is not None and pdb.is_file() and _pdb_natoms(pdb) == len(st.atoms.symbols):
            coords = _pdb_positions(pdb)
            if coords is None:
                errors.append(ValidationError("method.coordinates_file", L(
                    "PDB の原子座標を読めません", "cannot read the PDB atom coordinates")))
            elif any(abs(a - b) > 1e-3 for xyz, ref in zip(coords, st.atoms.positions)
                     for a, b in zip(xyz, ref)):
                errors.append(ValidationError("structure.atoms", L(
                    "Spec と PDB の原子座標または順序が違います。実際に計算する PDB を構造として読み込んでください",
                    "the spec and PDB differ in coordinates or atom order; load the actual simulation PDB as the structure")))
        if not m.parameter_files:
            errors.append(ValidationError("method.parameter_files", L(
                "CHARMM のパラメータファイルを指定してください", "provide CHARMM parameter files")))
        for path_value in m.parameter_files:
            if not Path(path_value).expanduser().is_file():
                errors.append(ValidationError("method.parameter_files", L(
                    f"パラメータファイルがありません: {path_value}", f"parameter file not found: {path_value}")))
        if not m.exclude:
            errors.append(ValidationError("method.exclude", L(
                "exclude の規則を明示してください", "explicitly set the exclude rule")))
        if m.one_four_scaling is None or m.one_four_scaling < 0:
            errors.append(ValidationError("method.one_four_scaling", L(
                "oneFourScaling を 0 以上の値で明示してください", "explicitly set a non-negative oneFourScaling")))
        if m.switching is None:
            errors.append(ValidationError("method.switching", L(
                "switching を使うかどうか明示してください", "explicitly choose whether switching is enabled")))
        if m.cutoff_ang <= 0 or m.pairlistdist_ang <= m.cutoff_ang:
            errors.append(ValidationError("method.cutoff_ang", L(
                "cutoff は正、pairlistdist は cutoff より大きい Å 値にしてください",
                "cutoff must be positive and pairlistdist must be larger than cutoff (both in Å)")))
        if m.switching and not (0 < m.switchdist_ang < m.cutoff_ang):
            errors.append(ValidationError("method.switchdist_ang", L(
                "switching を使うときは 0 < switchdist < cutoff を満たしてください",
                "with switching enabled, require 0 < switchdist < cutoff")))
        if m.switching is False and m.switchdist_ang != 0:
            errors.append(ValidationError("method.switchdist_ang", L(
                "switching を使わないときは switchdist_ang を 0 にしてください",
                "set switchdist_ang to 0 when switching is off")))
        if m.seed < 1:
            errors.append(ValidationError("method.seed", L("seed は 1 以上にしてください", "seed must be positive")))
        if ((spec.runtime.mpiprocs > 1 or spec.runtime.omp_threads > 1)
                and spec.runtime.profile in cfg.profiles
                and not cfg.profile(spec.runtime.profile).commands.get(self.code)):
            errors.append(ValidationError("runtime.profile", L(
                "並列実行では commands.namd に起動コマンドを指定してください",
                "set commands.namd explicitly for parallel execution")))
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st, md = spec.method, spec.structure, spec.task.md
        lines = [f"numsteps {md.steps}", "structure topology.psf", "coordinates coordinates.pdb",
                 f"temperature {md.temperature_k:.10g}", f"seed {m.seed}",
                 "outputName adit", f"DCDfreq {md.dump_interval}",
                 f"outputEnergies {md.dump_interval}", f"timestep {md.timestep_fs:.10g}"]
        lines += [f"parameters parameter_{i:02d}.prm" for i in range(1, len(m.parameter_files) + 1)]
        lines += [f"exclude {m.exclude}", f"oneFourScaling {m.one_four_scaling:.10g}",
                  f"switching {'on' if m.switching else 'off'}"]
        if m.switching:
            lines.append(f"switchdist {m.switchdist_ang:.10g}")
        lines += [f"cutoff {m.cutoff_ang:.10g}", f"pairlistdist {m.pairlistdist_ang:.10g}"]
        if st.periodic:
            lines += [f"cellBasisVector{i} {x:.10f} {y:.10f} {z:.10f}"
                      for i, (x, y, z) in enumerate(st.atoms.cell, 1)]
        lines.append("")
        return {"namd.conf": "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        m = spec.method
        copies = {"topology.psf": Path(m.structure_file).expanduser(),
                  "coordinates.pdb": Path(m.coordinates_file).expanduser()}
        copies.update({f"parameter_{i:02d}.prm": Path(path).expanduser()
                       for i, path in enumerate(m.parameter_files, 1)})
        return copies

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, "namd3").format(
            mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} namd.conf > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="namd3", files=[L(
            "  namd.conf / topology.psf / coordinates.pdb / parameter_*.prm   NAMD の NVE 入力と利用者の力場",
            "  namd.conf / topology.psf / coordinates.pdb / parameter_*.prm   NAMD NVE input and user-supplied force field")],
            prepare=[L("  PSF の原子電荷と、PSF・PDB・パラメータの原子順・力場を確認してください。NAMD は初速度を指定した温度と seed から作ります。ADIT は力場を作りません。",
                       "  Check atom order, charges and force field across PSF, PDB and parameter files. NAMD draws initial velocities from the specified temperature and seed. ADIT does not build force fields."),
                     L("  NVE では共通 Spec の熱浴・結合時定数・圧力・圧力浴時定数は使用しません。",
                       "  NVE does not use the common thermostat, coupling time, pressure, or barostat time.")],
            outputs=[L("  output.log / adit.dcd   NAMD の出力。ADIT はまだ内容を解析しません。",
                       "  output.log / adit.dcd   NAMD outputs; ADIT does not yet parse them.")])


register(NamdGenerator())
