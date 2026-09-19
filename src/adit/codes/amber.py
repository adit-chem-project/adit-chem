"""Amber/sander energy-minimization input with user-supplied topology and restart."""
# Official command and &cntrl example:
# https://ambermd.org/tutorials/basic/tutorial5/index.php
# Reference manuals: https://ambermd.org/Manuals.php

from __future__ import annotations

import math
import re
from pathlib import Path

from ase.geometry import cell_to_cellpar

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import AmberMethod, CalculationSpec
from adit.validate_types import ValidationError


def _prmtop_natoms(path: Path) -> int | None:
    # Read only the first POINTERS value (NATOM), without loading a large topology.
    try:
        with path.open(encoding="ascii", errors="replace") as stream:
            for line in stream:
                if line.strip().upper() != "%FLAG POINTERS":
                    continue
                next(stream)  # %FORMAT
                return int(next(stream).split()[0])
    except (OSError, StopIteration, ValueError, IndexError):
        pass
    return None


def _rst7_natoms(path: Path) -> int | None:
    try:
        with path.open(encoding="ascii", errors="replace") as stream:
            next(stream)  # title
            return int(next(stream).split()[0])
    except (OSError, StopIteration, ValueError, IndexError):
        return None


_NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


def _rst7_positions(path: Path, count: int) -> list[tuple[float, float, float]] | None:
    # Read only the coordinate prefix of a text restart; later velocities are untouched.
    try:
        values: list[float] = []
        with path.open(encoding="ascii", errors="replace") as stream:
            next(stream)
            next(stream)
            for line in stream:
                values.extend(float(token.replace("D", "E").replace("d", "e")) for token in _NUMBER.findall(line))
                if len(values) >= count * 3:
                    break
        if len(values) < count * 3 or not all(math.isfinite(x) for x in values[:count * 3]):
            return None
        return [tuple(values[i:i + 3]) for i in range(0, count * 3, 3)]
    except (OSError, StopIteration, ValueError):
        return None


def _rst7_box(path: Path, count: int) -> tuple[float, ...] | None:
    # A text restart has coordinates, optional velocities, then six cell parameters.
    try:
        remaining: list[float] = []
        seen = 0
        with path.open(encoding="ascii", errors="replace") as stream:
            next(stream)
            next(stream)
            for line in stream:
                for token in _NUMBER.findall(line):
                    value = float(token.replace("D", "E").replace("d", "e"))
                    seen += 1
                    if seen > 3 * count:
                        remaining.append(value)
                    if len(remaining) > 3 * count + 6:
                        return None
        if len(remaining) not in (6, 3 * count + 6):
            return None
        box = tuple(remaining[-6:])
        return box if all(math.isfinite(x) for x in box) else None
    except (OSError, StopIteration, ValueError):
        return None


class AmberGenerator(InputGenerator):
    code = "amber"
    uses_kpoints = False
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m, st = spec.method, spec.structure
        if not isinstance(m, AmberMethod):
            return [ValidationError("method.code", L("Amber の条件ではありません", "not an Amber method"))]
        errors: list[ValidationError] = []
        if spec.task.type != "geometry_optimization":
            errors.append(ValidationError("task.type", L(
                "この Amber 生成器はエネルギー最小化だけに対応します",
                "this Amber generator supports energy minimization only")))
        if spec.task.relax_cell != "no":
            errors.append(ValidationError("task.relax_cell", L(
                "Amber のこの生成器はセルを緩和しません", "this Amber generator cannot relax the cell")))
        if st.fixed_atoms or st.fixed_axes:
            errors.append(ValidationError("structure.fixed_atoms", L(
                "この Amber 生成器は固定原子・固定軸を入力へ写せません",
                "this Amber generator cannot write fixed atoms or axes")))
        if st.periodic and not all(st.atoms.pbc):
            errors.append(ValidationError("structure.atoms", L(
                "この Amber 生成器は 3 次元周期系か非周期系だけに対応します",
                "this Amber generator supports only fully periodic or nonperiodic systems")))
        if spec.kpoints is not None:
            errors.append(ValidationError("kpoints", L(
                "Amber は k 点を使用しません。指定を外してください",
                "Amber does not use k-points; remove this setting")))
        if spec.handoff is not None:
            errors.append(ValidationError("handoff", L(
                "この Amber 生成器は続きの計算に対応しません", "this Amber generator does not support restarts")))
        for key in ("topology_file", "coordinates_file"):
            value = getattr(m, key)
            path = Path(value).expanduser() if value.strip() else None
            if path is None or not path.is_file():
                errors.append(ValidationError(f"method.{key}", L(
                    f"{key} に既存のファイルを指定してください", f"set {key} to an existing file")))
        top = Path(m.topology_file).expanduser() if m.topology_file.strip() else None
        crd = Path(m.coordinates_file).expanduser() if m.coordinates_file.strip() else None
        if top and top.is_file():
            n = _prmtop_natoms(top)
            if n is None:
                errors.append(ValidationError("method.topology_file", L(
                    "prmtop の %FLAG POINTERS から原子数を読めません", "cannot read NATOM from the prmtop %FLAG POINTERS")))
            elif n != len(st.atoms.symbols):
                errors.append(ValidationError("structure.atoms", L(
                    f"Spec は {len(st.atoms.symbols)} 原子、prmtop は {n} 原子です",
                    f"the spec has {len(st.atoms.symbols)} atoms but the prmtop has {n}")))
        if crd and crd.is_file():
            n = _rst7_natoms(crd)
            if n is None:
                errors.append(ValidationError("method.coordinates_file", L(
                    "座標は原子数を読めるテキスト形式の rst7/inpcrd にしてください",
                    "coordinates must be a text rst7/inpcrd with a readable atom count")))
            elif n != len(st.atoms.symbols):
                errors.append(ValidationError("structure.atoms", L(
                    f"Spec は {len(st.atoms.symbols)} 原子、rst7/inpcrd は {n} 原子です",
                    f"the spec has {len(st.atoms.symbols)} atoms but the rst7/inpcrd has {n}")))
            else:
                coords = _rst7_positions(crd, n)
                if coords is None:
                    errors.append(ValidationError("method.coordinates_file", L(
                        "rst7/inpcrd の座標を読めません", "cannot read the rst7/inpcrd coordinates")))
                elif any(abs(a - b) > 1e-3 for xyz, ref in zip(coords, st.atoms.positions)
                         for a, b in zip(xyz, ref)):
                    errors.append(ValidationError("structure.atoms", L(
                        "Spec と rst7/inpcrd の原子座標または順序が違います。実際に計算する座標の構造を読み込んでください",
                        "the spec and rst7/inpcrd differ in coordinates or atom order; load the structure actually used for the run")))
                if st.periodic:
                    box = _rst7_box(crd, n)
                    if box is None:
                        errors.append(ValidationError("method.coordinates_file", L(
                            "周期系の rst7/inpcrd にセルの長さと角度がありません",
                            "the periodic rst7/inpcrd lacks cell lengths and angles")))
                    elif any(abs(actual - expected) > 1e-3 for actual, expected in zip(box, cell_to_cellpar(st.atoms.cell))):
                        errors.append(ValidationError("structure.atoms", L(
                            "Spec と rst7/inpcrd のセルが違います。実際の計算に使うセルを構造へ読み込んでください",
                            "the spec and rst7/inpcrd cells differ; load the cell actually used for the run")))
        if m.cutoff_ang <= 0:
            errors.append(ValidationError("method.cutoff_ang", L(
                "cut を正の Å 値で明示してください", "explicitly set a positive cut in Å")))
        if (not st.periodic and m.igb is None) or (m.igb is not None and (m.igb < 0 or (st.periodic and m.igb != 0))):
            errors.append(ValidationError("method.igb", L(
                "非周期系では igb を明示してください (真空は 0)。周期系では 0 か空欄にしてください",
                "explicitly set igb for a nonperiodic system (0 for vacuum); use 0 or leave it empty for a periodic system")))
        if ((spec.runtime.mpiprocs > 1 or spec.runtime.omp_threads > 1)
                and spec.runtime.profile in cfg.profiles
                and not cfg.profile(spec.runtime.profile).commands.get(self.code)):
            errors.append(ValidationError("runtime.profile", L(
                "並列実行では commands.amber に起動コマンドを指定してください",
                "set commands.amber explicitly for parallel execution")))
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, task, st = spec.method, spec.task, spec.structure
        lines = ["ADIT energy minimization", "&cntrl", "  imin=1,",
                 f"  maxcyc={task.max_steps},",
                 f"  ntb={1 if st.periodic else 0},",
                 f"  cut={m.cutoff_ang:.10g},"]
        if not st.periodic:
            lines.append(f"  igb={m.igb},")
        lines += ["/", ""]
        return {"amber.in": "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        m = spec.method
        return {"topology.prmtop": Path(m.topology_file).expanduser(),
                "coordinates.rst7": Path(m.coordinates_file).expanduser()}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, "sander").format(
            mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} -O -i amber.in -o output.log -p topology.prmtop -c coordinates.rst7 -r final.rst7"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="sander", files=[L(
            "  amber.in / topology.prmtop / coordinates.rst7   Amber の最小化入力と利用者のトポロジー・座標",
            "  amber.in / topology.prmtop / coordinates.rst7   Amber minimization input and user-supplied topology/coordinates")],
            prepare=[L("  prmtop の原子順・力場・電荷、および rst7/inpcrd の座標と周期セルを確認してください。ADIT は原子型や力場を作りません。",
                       "  Check atom order, force field and charges in the prmtop, and coordinates/periodic cell in rst7/inpcrd. ADIT does not build atom types or force fields."),
                     L("  共通 Spec の最適化器と力の閾値は amber.in に適用していません。",
                       "  The common optimizer and force threshold are not applied to amber.in.")],
            outputs=[L("  output.log / final.rst7   Amber の出力。ADIT はまだ内容を解析しません。",
                       "  output.log / final.rst7   Amber outputs; ADIT does not yet parse them.")])


register(AmberGenerator())
