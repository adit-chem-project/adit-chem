"""OpenMX 4.0 neutral, non-spin-polarized single-point inputs."""
# Official keywords and file layout:
# https://openmx-square.org/openmx_man4.0/s8_2_keywords.html
# https://openmx-square.org/openmx_man4.0/s13_1_valence-pseudo.html

from __future__ import annotations

import re
from pathlib import Path

from ase.data import atomic_numbers
import numpy as np

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, OpenmxMethod
from adit.validate_types import ValidationError


_NAME = re.compile(r"[A-Za-z0-9_.+-]+\Z")


class OpenmxGenerator(InputGenerator):
    code = "openmx"
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None  # PAO and VPS stay in the user's own DFT_DATA directory.

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m, st = spec.method, spec.structure
        if not isinstance(m, OpenmxMethod):
            return [ValidationError("method.code", L("OpenMX の条件ではありません", "not an OpenMX method"))]
        errors: list[ValidationError] = []
        if spec.task.type != "single_point":
            errors.append(ValidationError("task.type", L(
                "この OpenMX 生成器は一点計算だけに対応します", "this OpenMX generator supports single-point calculations only")))
        if st.charge != 0 or st.multiplicity != 1:
            errors.append(ValidationError("structure.charge", L(
                "この OpenMX 生成器は中性・一重項だけに対応します",
                "this OpenMX generator supports neutral singlets only")))
        if st.fixed_atoms or st.fixed_axes:
            errors.append(ValidationError("structure.fixed_atoms", L(
                "この OpenMX 生成器は固定原子・固定軸を入力に書けません",
                "this OpenMX generator cannot write fixed atoms or axes")))
        if st.periodic and not all(st.atoms.pbc):
            errors.append(ValidationError("structure.atoms", L(
                "この OpenMX 生成器は 3 次元周期系か非周期系だけに対応します",
                "this OpenMX generator supports only fully periodic or nonperiodic systems")))
        if not st.periodic and np.any(st.atoms.cell) and abs(np.linalg.det(st.atoms.cell)) < 1e-6:
            errors.append(ValidationError("structure.atoms", L(
                "非周期系のセルを指定する場合は、3 本の独立した格子ベクトルが必要です",
                "a specified nonperiodic cell needs three independent lattice vectors")))
        if spec.handoff is not None:
            errors.append(ValidationError("handoff", L("この OpenMX 生成器は続きの計算に対応しません", "this OpenMX generator does not support restarts")))
        if not m.data_path.strip() or any(c.isspace() for c in m.data_path):
            errors.append(ValidationError("method.data_path", L(
                "実行先の DATA.PATH を空白・改行のないパスで指定してください",
                "set DATA.PATH on the execution host to a path without spaces or newlines")))
        if not m.xc:
            errors.append(ValidationError("method.xc", L("OpenMX の XC を明示してください", "explicitly set the OpenMX XC functional")))
        if m.energycutoff_ry <= 0:
            errors.append(ValidationError("method.energycutoff_ry", L(
                "scf.energycutoff を正の Ry 値で明示してください", "explicitly set a positive scf.energycutoff in Ry")))
        symbols = set(st.atoms.symbols)
        for symbol in sorted(symbols):
            if not m.pao.get(symbol) or not _NAME.fullmatch(m.pao[symbol]):
                errors.append(ValidationError("method.pao", L(
                    f"{symbol} の PAO 名を指定してください", f"set the PAO name for {symbol}")))
            if not m.vps.get(symbol) or not _NAME.fullmatch(m.vps[symbol]):
                errors.append(ValidationError("method.vps", L(
                    f"{symbol} の VPS 名を指定してください", f"set the VPS name for {symbol}")))
            value = m.valence.get(symbol, 0)
            if value <= 0 or value > atomic_numbers[symbol]:
                errors.append(ValidationError("method.valence", L(
                    f"{symbol} の VPS に対応する価電子数を明示してください",
                    f"explicitly set the number of valence electrons for the {symbol} VPS")))
        if st.periodic and spec.kpoints is not None and any(spec.kpoints.shift):
            errors.append(ValidationError("kpoints.shift", L(
                "OpenMX のこの生成器では k 点の shift を入力に書けません",
                "this OpenMX generator cannot write a shifted k-point mesh")))
        if ((spec.runtime.mpiprocs > 1 or spec.runtime.omp_threads > 1)
                and spec.runtime.profile in cfg.profiles
                and not cfg.profile(spec.runtime.profile).commands.get(self.code)):
            errors.append(ValidationError("runtime.profile", L(
                "並列実行では commands.openmx に起動コマンドを指定してください",
                "set commands.openmx explicitly for parallel execution")))
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st = spec.method, spec.structure
        symbols = sorted(set(st.atoms.symbols))
        lines = ["System.Name adit", f"DATA.PATH {m.data_path}",
                 f"Species.Number {len(symbols)}", "<Definition.of.Atomic.Species"]
        lines += [f"{s} {m.pao[s]} {m.vps[s]}" for s in symbols]
        lines += ["Definition.of.Atomic.Species>", f"Atoms.Number {len(st.atoms.symbols)}",
                  "Atoms.SpeciesAndCoordinates.Unit Ang", "<Atoms.SpeciesAndCoordinates"]
        lines += [f"{i} {s} {x:.10f} {y:.10f} {z:.10f} {m.valence[s] / 2:g} {m.valence[s] / 2:g}"
                  for i, (s, (x, y, z)) in enumerate(zip(st.atoms.symbols, st.atoms.positions), 1)]
        lines += ["Atoms.SpeciesAndCoordinates>"]
        if st.periodic or np.any(st.atoms.cell):
            lines += ["Atoms.UnitVectors.Unit Ang", "<Atoms.UnitVectors"]
            lines += [f"{x:.10f} {y:.10f} {z:.10f}" for x, y, z in st.atoms.cell]
            lines += ["Atoms.UnitVectors>"]
        if st.periodic:
            lines += ["scf.EigenvalueSolver Band"]
            mesh = spec.kpoints.resolved_mesh(st.atoms.cell)
            lines.append(f"scf.Kgrid {mesh[0]} {mesh[1]} {mesh[2]}")
        else:
            lines.append("scf.EigenvalueSolver Cluster")
        lines += ["scf.SpinPolarization OFF", f"scf.XcType {m.xc}",
                  f"scf.energycutoff {m.energycutoff_ry:.10g}", ""]
        return {"openmx.dat": "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, "openmx").format(
            mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} openmx.dat > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="openmx", files=[L("  openmx.dat   OpenMX の一点計算の入力", "  openmx.dat   OpenMX single-point input")],
                           prepare=[L("  DATA.PATH の PAO/VPS が入力の名前・価電子数と一致するか確認してください。ADIT はこれらを同梱・検証しません。",
                                      "  Check that the PAO/VPS in DATA.PATH match the names and valence counts in the input. ADIT does not package or verify them.")],
                           outputs=[L("  output.log   OpenMX の標準出力。ADIT はまだ内容を解析しません。",
                                      "  output.log   OpenMX output; ADIT does not yet parse it.")])


register(OpenmxGenerator())
