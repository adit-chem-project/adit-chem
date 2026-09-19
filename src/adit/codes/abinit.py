
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
from ase.data import atomic_numbers

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import AbinitMethod, CalculationSpec
from adit.validate_types import ValidationError

INPUT_FILE = "input.abi"
DEFAULT_COMMAND = "abinit"
BOHR_ANG = 0.529177210903
HARTREE_EV = 27.211386245988
EV_ANG_TO_HA_BOHR = BOHR_ANG / HARTREE_EV
_NAME = re.compile(r"[a-zA-Z][a-zA-Z0-9_]*\Z")
DOC = "https://docs.abinit.org/variables/"


def _num(x: float) -> str:
    return f"{x:.12g}"


class AbinitGenerator(InputGenerator):
    code = "abinit"

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, AbinitMethod):
            return [ValidationError("method.code", L(f"ABINIT の生成器に {m.code!r} の手法が渡されました",
                                                     f"the ABINIT generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        if not all(st.atoms.pbc):
            errs.append(ValidationError("structure.atoms", L(
                "この ABINIT 生成器は 3 方向とも周期の系だけに対応します (分子は大きな箱に入れてください)",
                "this ABINIT generator supports systems periodic in all three directions (put molecules in a large box)")))
        if st.charge != 0:
            errs.append(ValidationError("structure.charge", L(
                "この ABINIT 生成器は中性の系だけです (帯電系は補正の指定が要ります)",
                "this ABINIT generator supports neutral systems only (charged cells need a compensating scheme)")))
        if st.multiplicity != 1:
            errs.append(ValidationError("structure.multiplicity", L(
                "この ABINIT 生成器はスピン分極 (nsppol 2) を扱いません。多重度は 1 のままにしてください",
                "this ABINIT generator does not handle spin polarization (nsppol 2); leave the multiplicity at 1")))
        if st.velocities is not None:
            errs.append(ValidationError("structure.velocities", L(
                "この ABINIT 生成器は初速度を入力へ書けません (MD に対応していません)",
                "this ABINIT generator cannot write initial velocities (MD is not supported)")))
        if spec.handoff is not None:
            errs.append(ValidationError("handoff", L(
                "この ABINIT 生成器は続きの計算に対応しません", "this ABINIT generator does not support restarts")))
        if t.type in ("molecular_dynamics", "vibrations", "band_structure"):
            errs.append(ValidationError("task.type", L(
                "この ABINIT 生成器は一点計算と構造最適化だけです (MD・振動 (anaddb)・バンドは別の入力が要ります)",
                "this ABINIT generator supports single point and geometry optimization only (MD, phonons via anaddb and bands need separate inputs)")))
        if t.type == "geometry_optimization":
            if t.relax_cell != "no":
                errs.append(ValidationError("task.relax_cell", L(
                    "セルを動かす最適化 (optcell) は、ecutsm と dilatmx を利用者が決める必要があるので、この生成器では扱いません",
                    "cell optimization (optcell) needs user-chosen ecutsm and dilatmx, so this generator does not handle it")))
            if t.max_steps < 1:
                errs.append(ValidationError("task.max_steps", L(
                    "構造最適化の最大回数 (ntime) は 1 以上にしてください", "the maximum number of steps (ntime) must be at least 1")))
            if t.force_tolerance_ev_per_ang <= 0:
                errs.append(ValidationError("task.force_tolerance_ev_per_ang", L(
                    "力の収束基準 (tolmxf) は 0 より大きい値にしてください", "the force tolerance (tolmxf) must be greater than 0")))
            if m.tolerance == "toldfe":
                errs.append(ValidationError("method.tolerance", L(
                    "ABINIT は構造最適化 (ionmov) で toldfe を受け付けません。toldff・tolrff・tolvrs (または tolwfr) から選んでください",
                    "ABINIT does not accept toldfe with a structural optimization (ionmov); choose toldff, tolrff or tolvrs (or tolwfr)")))
        if not m.tolerance:
            errs.append(ValidationError("method.tolerance", L(
                "SCF の収束の種類 (toldfe / toldff / tolrff / tolvrs / tolwfr) を 1 つ選んでください。ABINIT は 1 つだけ受け取ります",
                "choose one SCF convergence criterion (toldfe / toldff / tolrff / tolvrs / tolwfr); ABINIT accepts exactly one")))
        elif m.tolerance_value <= 0:
            errs.append(ValidationError("method.tolerance_value", L(
                "収束の閾値は 0 より大きい値にしてください", "the convergence threshold must be greater than 0")))
        if m.ecut_ha <= 0:
            errs.append(ValidationError("method.ecut_ha", L(
                "平面波のカットオフ ecut [Hartree] を指定してください (ABINIT に既定はありません)",
                "give the plane-wave cut-off ecut in hartree (ABINIT has no default)")))
        if m.nstep < 1:
            errs.append(ValidationError("method.nstep", L("SCF の最大反復回数は 1 以上にしてください", "the maximum number of SCF steps must be at least 1")))
        if m.tsmear_ha < 0:
            errs.append(ValidationError("method.tsmear_ha", L("広がり tsmear は 0 以上にしてください", "tsmear must be non-negative")))
        if m.nband is not None and m.nband < 1:
            errs.append(ValidationError("method.nband", L("バンド数は 1 以上にしてください", "the number of bands must be at least 1")))
        if m.pawecutdg_ha < 0 or (m.pawecutdg_ha and m.pawecutdg_ha < m.ecut_ha):
            errs.append(ValidationError("method.pawecutdg_ha", L(
                "pawecutdg は ecut 以上にしてください", "pawecutdg must be at least as large as ecut")))
        elements = sorted(set(st.atoms.symbols))
        missing = [e for e in elements if not m.pseudos.get(e, "").strip()]
        if missing:
            errs.append(ValidationError("method.pseudos", L(
                f"擬ポテンシャルのファイルが指定されていない元素があります: {', '.join(missing)}。ADIT は擬ポテンシャルを持ちません",
                f"no pseudopotential file was given for: {', '.join(missing)}; ADIT does not ship pseudopotentials")))
        for element in elements:
            value = m.pseudos.get(element, "").strip()
            if value and not Path(value).expanduser().is_file():
                errs.append(ValidationError("method.pseudos", L(
                    f"{element} の擬ポテンシャルがありません: {value}", f"pseudopotential file for {element} not found: {value}")))
        names = [Path(v).name for v in m.pseudos.values() if v.strip()]
        if len(set(names)) != len(names):
            errs.append(ValidationError("method.pseudos", L(
                "擬ポテンシャルのファイル名が重なっています。生成したファイルでは同じ名前で写すので、名前を分けてください",
                "two pseudopotential files share a name; they are copied under the same names, so rename them")))
        if st.fixed_axes and all(st.atoms.pbc):
            from adit.vasp_constraints import VaspConstraintError

            try:
                self._reduced_axis_masks(spec)
            except VaspConstraintError as ex:
                errs.append(ValidationError("structure.fixed_axes", L(
                    f"ABINIT の iatfixx/y/z は格子ベクトル方向 (還元座標) の固定です。{ex}",
                    f"ABINIT's iatfixx/y/z fix lattice-vector (reduced) directions. {ex}")))
        if spec.kpoints is None:
            errs.append(ValidationError("kpoints", L("周期系では k 点を指定してください", "give the k-points for a periodic system")))
        elif any(x <= 0 for x in spec.kpoints.resolved_mesh(st.atoms.cell)):
            errs.append(ValidationError("kpoints.mesh", L("k 点の分割数は 1 以上にしてください", "each k-point division must be at least 1")))
        bad = [k for k in m.extra if not _NAME.fullmatch(k.strip())]
        if bad:
            errs.append(ValidationError("method.extra", L(
                f"ABINIT の変数名として読めません: {bad}", f"not valid ABINIT variable names: {bad}")))
        overlap = sorted(set(k.strip() for k in m.extra) & self._written_names(spec))
        if overlap:
            errs.append(ValidationError("method.extra", L(
                f"この生成器が書く変数と重なっています: {', '.join(overlap)}。画面の欄で指定してください",
                f"these clash with variables this generator already writes: {', '.join(overlap)}; set them through their own fields")))
        if not all(math.isfinite(x) for x in np.asarray(st.atoms.cell, dtype=float).ravel()):
            errs.append(ValidationError("structure.atoms", L("セルに数でない値があります", "the cell contains non-numeric values")))
        return errs

    def _written_names(self, spec: CalculationSpec) -> set[str]:
        m = spec.method
        names = {"acell", "rprim", "natom", "ntypat", "znucl", "typat", "xred", "pseudos", "pp_dirpath",
                 "ecut", "nstep", "ngkpt", "nshiftk", "shiftk", "chkprim", m.tolerance or "toldfe"}
        if spec.task.type == "geometry_optimization":
            names |= {"ionmov", "ntime", "tolmxf"}
        if spec.structure.fixed_atoms or spec.structure.fixed_axes:
            names |= {"natfix", "iatfix", "natfixx", "iatfixx", "natfixy", "iatfixy", "natfixz", "iatfixz"}
        for key, value in (("pawecutdg", m.pawecutdg_ha), ("tsmear", m.tsmear_ha)):
            if value:
                names.add(key)
        for key, value in (("ixc", m.ixc), ("occopt", m.occopt), ("nband", m.nband)):
            if value is not None:
                names.add(key)
        return names

    def _types(self, spec: CalculationSpec) -> tuple[list[str], list[int]]:
        order: list[str] = []
        for symbol in spec.structure.atoms.symbols:
            if symbol not in order:
                order.append(symbol)
        return order, [order.index(s) + 1 for s in spec.structure.atoms.symbols]

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        m, st, t = spec.method, spec.structure, spec.task
        cell = np.asarray(st.atoms.cell, dtype=float)
        scaled = np.asarray(spec.atoms.get_scaled_positions(), dtype=float)
        order, typat = self._types(spec)
        mesh = spec.kpoints.resolved_mesh(cell)
        shift = spec.kpoints.shift
        lines = ["# ADIT が生成した ABINIT の入力 / ABINIT input generated by ADIT",
                 "# 長さの単位は Bohr (ABINIT の既定) / lengths are in bohr (the ABINIT default)",
                 "",
                 "# ---- 構造 / structure ----",
                 "acell 1.0 1.0 1.0",
                 "rprim"]
        lines += ["  " + " ".join(_num(x / BOHR_ANG) for x in row) for row in cell]
        lines += [f"natom {len(typat)}", f"ntypat {len(order)}",
                  "znucl " + " ".join(str(atomic_numbers[s]) for s in order),
                  "typat " + " ".join(str(i) for i in typat), "xred"]
        lines += ["  " + " ".join(_num(x) for x in row) for row in scaled]
        lines += ["", "# ---- 擬ポテンシャル / pseudopotentials ----", 'pp_dirpath "./"',
                  'pseudos "' + ", ".join(Path(m.pseudos[s]).name for s in order) + '"',
                  "", "# ---- 計算の条件 / calculation ----", f"ecut {_num(m.ecut_ha)}"]
        if m.pawecutdg_ha:
            lines.append(f"pawecutdg {_num(m.pawecutdg_ha)}")
        if m.ixc is not None:
            lines.append(f"ixc {m.ixc}")
        if m.occopt is not None:
            lines.append(f"occopt {m.occopt}")
        if m.tsmear_ha:
            lines.append(f"tsmear {_num(m.tsmear_ha)}")
        if m.nband is not None:
            lines.append(f"nband {m.nband}")
        lines += ["chkprim 0", f"nstep {m.nstep}", f"{m.tolerance} {_num(m.tolerance_value)}", "",
                  "# ---- k 点 / k-points ----",
                  "ngkpt " + " ".join(str(x) for x in mesh), "nshiftk 1",
                  "shiftk " + " ".join(_num(x) for x in shift)]
        if t.type == "geometry_optimization":
            lines += ["", "# ---- 構造最適化 / geometry optimization ----", "ionmov 2", f"ntime {t.max_steps}",
                      f"tolmxf {_num(t.force_tolerance_ev_per_ang * EV_ANG_TO_HA_BOHR)}"]
        fixed = self._fixed_lines(spec)
        if fixed:
            lines += ["", "# ---- 動かさない原子 / fixed atoms ----"] + fixed
        if m.extra:
            lines += ["", "# ---- 画面にない変数 (利用者の指定) / extra variables set by the user ----"]
            lines += [f"{k.strip()} {v}" for k, v in m.extra.items()]
        lines.append("")
        return {INPUT_FILE: "\n".join(lines)}

    @staticmethod
    def _reduced_axis_masks(spec: CalculationSpec) -> dict[int, np.ndarray]:
        # iatfixx/y/z act on reduced coordinates for ionmov /= 1 (docs.abinit.org iatfix), so the
        # Cartesian masks are converted the same way as VASP's Selective dynamics flags.
        from adit.vasp_constraints import cartesian_to_direct_mask

        st = spec.structure
        full = set(st.fixed_atoms)
        masks: dict[int, np.ndarray] = {}
        for key, move in st.fixed_axes.items():
            index = int(key)
            if index in full:
                continue
            masks[index] = cartesian_to_direct_mask(st.atoms.cell, [not c for c in move])
        return masks

    def _fixed_lines(self, spec: CalculationSpec) -> list[str]:
        st = spec.structure
        full = sorted(set(st.fixed_atoms))
        per_axis: dict[int, list[int]] = {0: [], 1: [], 2: []}
        for index, mask in self._reduced_axis_masks(spec).items():
            for axis in range(3):
                if mask[axis]:
                    per_axis[axis].append(index)
        lines = []
        if full:
            lines += [f"natfix {len(full)}", "iatfix " + " ".join(str(i + 1) for i in full)]
        for axis, name in enumerate("xyz"):
            picked = sorted(per_axis[axis])
            if picked:
                lines += [f"natfix{name} {len(picked)}", f"iatfix{name} " + " ".join(str(i + 1) for i in picked)]
        return lines

    def files_to_copy(self, spec: CalculationSpec, res) -> dict[str, Path]:
        m = spec.method
        order, _ = self._types(spec)
        return {Path(m.pseudos[s]).name: Path(m.pseudos[s]).expanduser() for s in order if m.pseudos.get(s, "").strip()}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, DEFAULT_COMMAND).format(
            mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} {INPUT_FILE} > output.log 2>&1"

    def version_probe(self, spec):
        return ("output.log", "Version .* of ABINIT")

    # ---- README ----
    def readme_notes(self, spec: CalculationSpec, res, copies: dict[str, Path]) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        files = [L("  input.abi    ABINIT の入力 (構造・擬ポテンシャルの指定・計算の条件・k 点)",
                   "  input.abi    ABINIT input (structure, pseudopotential names, run settings, k-points)")]
        if copies:
            files.append(L(f"  {', '.join(sorted(copies))}   擬ポテンシャル (利用者が指定したものを写しました)",
                           f"  {', '.join(sorted(copies))}   pseudopotentials (copies of the files you selected)"))
        prepare = [
            L("  擬ポテンシャルは利用者のものです。ADIT は同梱も検証もしません。汎関数は擬ポテンシャルに書かれたものを ABINIT が使います"
              + (f" (この入力は ixc {m.ixc} を明示しています)。" if m.ixc is not None else " (ixc を書いていません)。"),
              "  The pseudopotentials are yours; ADIT neither ships nor verifies them. ABINIT uses the functional recorded in the pseudopotential"
              + (f" (this input sets ixc {m.ixc})." if m.ixc is not None else " (no ixc is written).")),
            L("  ecut と k 点は収束を確かめてから決めてください。ADIT は推奨値を持ちません (1 つの値を振る一括生成が使えます)。",
              "  Converge ecut and the k-point mesh yourself; ADIT has no recommended values (the scan feature can vary one value)."),
            L("  入力には chkprim 0 を書いています (ADIT は利用者が与えたセルをそのまま使うため。スラブや超格子は"
              "プリミティブではないので、既定の chkprim 1 では ABINIT が止まります)。",
              "  The input sets chkprim 0 because ADIT uses the cell you supplied as it is; slabs and supercells are not primitive, "
              "and the ABINIT default chkprim 1 would stop the run."),
        ]
        if m.pawecutdg_ha:
            prepare.append(L("  pawecutdg を書いています (PAW の擬ポテンシャル用)。ノルム保存の擬ポテンシャルでは要りません。",
                             "  pawecutdg is written (for PAW datasets); norm-conserving pseudopotentials do not need it."))
        outputs = [L("  output.log   ABINIT の標準出力。ADIT は etotal (全エネルギー、Hartree) を読みます",
                     "  output.log   ABINIT standard output; ADIT reads etotal (total energy, hartree)"),
                   L("  input.abo    ABINIT の主出力ファイル (同じ内容の表と、最後の構造)",
                     "  input.abo    the main ABINIT output file (the same tables and the final structure)")]
        if t == "geometry_optimization":
            outputs.append(L("  構造最適化では、各段階のエネルギーと力が output.log に並びます (ADIT はエネルギーの推移として読みます)",
                             "  For a geometry optimization, the energy and forces of each step appear in output.log (ADIT reads them as the energy trace)"))
        return ReadmeNotes(program=DEFAULT_COMMAND, files=files, prepare=prepare, outputs=outputs)


register(AbinitGenerator())
