
from __future__ import annotations

from pathlib import Path

import numpy as np

from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.cp2k_data import Cp2kData
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import HARTREE_PER_BOHR_IN_EV_PER_ANG, CalculationSpec, Cp2kMethod
from adit.validate import electron_parity_error
from adit.validate_types import ValidationError

INPUT_FILE = "cp2k.inp"
BASIS_OUT = "BASIS_adit"
POTENTIAL_OUT = "POTENTIAL_adit"
PROJECT = "adit"
DEFAULT_COMMAND = "mpirun -np {mpiprocs} cp2k.psmp"
D3_FILE = "dftd3.dat"
EXTRA_SECTIONS = ("GLOBAL", "FORCE_EVAL", "FORCE_EVAL/DFT", "FORCE_EVAL/DFT/SCF", "FORCE_EVAL/DFT/XC", "FORCE_EVAL/SUBSYS",
                  "MOTION", "MOTION/GEO_OPT", "MOTION/CELL_OPT", "MOTION/MD")
THERMOSTATS = {"nose_hoover": "NOSE", "csvr": "CSVR"}


def _periodic_letters(pbc) -> str:
    s = "".join(a for a, p in zip("XYZ", pbc) if p)
    return s or "NONE"


def _f(x: float) -> str:
    return f"{x:.10g}"


def _data(cfg: Config) -> Cp2kData:
    return Cp2kData.discover(cfg.cp2k_data)


def _pick(m: Cp2kMethod, data: Cp2kData, element: str, *, potential: bool) -> str | None:
    given = (m.potential if potential else m.basis).get(element, "").strip()
    if given:
        return given
    names = data.names_for(m.potential_file if potential else m.basis_file, element, potential=potential)
    return names[0] if len(names) == 1 else None


def _data_problem(cfg: Config, data: Cp2kData) -> str:
    from adit.config import config_path
    if cfg.cp2k_data:
        return L(f"cp2k_data に書かれたフォルダがありません: {data.root}。環境設定ファイル {config_path()} の cp2k_data = \"...\" を直してください",
                 f"the folder given as cp2k_data does not exist: {data.root}. Fix cp2k_data = \"...\" in the settings file {config_path()}")
    return ""


class Cp2kGenerator(InputGenerator):
    code = "cp2k"

    def resolve(self, spec: CalculationSpec, cfg: Config) -> Cp2kData:
        return _data(cfg)

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, Cp2kMethod):
            return [ValidationError("method.code", L(f"CP2K の生成器に {m.code!r} の手法が渡されました", f"the CP2K generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        if not m.xc.strip():
            errs.append(ValidationError("method.xc", L("汎関数の名前を入れてください (CP2K の XC_FUNCTIONAL。例 PBE、BLYP、PBE0)",
                                                       "enter the functional (CP2K XC_FUNCTIONAL, e.g. PBE, BLYP, PBE0)")))
        elif len(m.xc.split()) != 1:
            errs.append(ValidationError("method.xc", L("汎関数の名前は空白を含まない 1 語で書いてください", "the functional must be a single word")))
        if m.cutoff_ry <= 0:
            errs.append(ValidationError("method.cutoff_ry", L("値を入れてください (MGRID の CUTOFF [Ry]。CP2K のマニュアルの「CUTOFF と REL_CUTOFF の収束」の手順で決めます)",
                                                              "enter a value (MGRID CUTOFF in Ry; the CP2K manual page on converging CUTOFF and REL_CUTOFF shows how to choose it)")))
        if m.rel_cutoff_ry <= 0:
            errs.append(ValidationError("method.rel_cutoff_ry", L("値を入れてください (MGRID の REL_CUTOFF [Ry])", "enter a value (MGRID REL_CUTOFF in Ry)")))
        if m.eps_scf <= 0:
            errs.append(ValidationError("method.eps_scf", L("0 より大きい値が必要です", "must be greater than 0")))
        if m.max_scf < 1:
            errs.append(ValidationError("method.max_scf", L("1 以上が必要です", "must be at least 1")))
        if st.multiplicity != 1 and not m.uks:
            errs.append(ValidationError("method.uks", L("多重度が 1 でないなら UKS (スピン分極) が必要です", "UKS (spin polarization) is required when the multiplicity is not 1")))
        bad = [k for k in m.extra_sections if k.strip().upper() not in EXTRA_SECTIONS]
        if bad:
            errs.append(ValidationError("method.extra_sections", L(f"追加の行を足せる節は {', '.join(EXTRA_SECTIONS)} です: {bad}",
                                                                   f"extra lines can be added to {', '.join(EXTRA_SECTIONS)}: {bad}")))
        if not all(st.atoms.pbc):
            if not m.poisson_solver:
                errs.append(ValidationError("method.poisson_solver", L(
                    f"周期でない方向 ({'/'.join(a for a, p in zip('xyz', st.atoms.pbc) if not p)}) があるので、ポアソン方程式の解き方 (POISSON_SOLVER: MT / WAVELET / ANALYTIC / MULTIPOLE) を選んでください",
                    f"there are non-periodic directions ({'/'.join(a for a, p in zip('xyz', st.atoms.pbc) if not p)}); choose a POISSON_SOLVER (MT / WAVELET / ANALYTIC / MULTIPOLE)")))
            if not any(st.atoms.pbc) and abs(np.linalg.det(np.asarray(st.atoms.cell, dtype=float))) < 1e-6 and m.isolated_box_ang <= 0:
                errs.append(ValidationError("method.isolated_box_ang", L(
                    "分子 (非周期) の構造に箱がありません。CP2K は箱の中で計算するので、箱の一辺 [Å] を入れてください (分子は箱の中心に置かれます)",
                    "the molecule (non-periodic) has no box; CP2K computes inside a box, so enter the box edge in Å (the molecule is centered)")))
        if spec.kpoints is not None and spec.kpoints.mode != "gamma" and any(s != 0 for s in spec.kpoints.shift):
            errs.append(ValidationError("kpoints.shift", L("CP2K の生成器は k 点のずらし (shift) に対応していません (0 にしてください。MONKHORST-PACK の格子を書きます)",
                                                           "the CP2K generator does not support k-point shifts (use 0; a MONKHORST-PACK grid is written)")))
        if t.type == "band_structure":
            errs.append(ValidationError("task.type", L("CP2K の生成器にバンド計算はありません", "the CP2K generator has no band structure")))
        if t.type == "geometry_optimization" and t.relax_cell == "volume_only":
            errs.append(ValidationError("task.relax_cell", L("CP2K の CELL_OPT には体積だけを変える指定がありません (形と体積、または格子を動かさない、から選んでください)",
                                                             "CP2K CELL_OPT cannot change only the volume (choose shape and volume, or no cell relaxation)")))
        if t.type == "vibrations" and (st.fixed_atoms or st.fixed_axes):
            errs.append(ValidationError("structure.fixed_atoms", L("CP2K の振動解析 (VIBRATIONAL_ANALYSIS) は固定原子を扱いません", "CP2K VIBRATIONAL_ANALYSIS does not handle fixed atoms")))
        errs += self._check_magnetism(spec)
        if t.type == "molecular_dynamics":
            md = t.md
            if md.ensemble != "NVE" and md.thermostat not in (*THERMOSTATS, "langevin"):
                errs.append(ValidationError("task.md.thermostat", L(f"CP2K にはない熱浴です: {md.thermostat} (nose_hoover / csvr / langevin)",
                                                                    f"thermostat not available in CP2K: {md.thermostat} (nose_hoover / csvr / langevin)")))
            if md.ensemble == "NPT" and md.thermostat == "langevin":
                errs.append(ValidationError("task.md.thermostat", L("CP2K の Langevin は NVT (ENSEMBLE LANGEVIN) だけです。NPT では nose_hoover か csvr にしてください",
                                                                    "CP2K Langevin is NVT only (ENSEMBLE LANGEVIN); use nose_hoover or csvr for NPT")))
        errs += self._check_data(spec, cfg)
        return errs

    @staticmethod
    def _check_magnetism(spec: CalculationSpec) -> list[ValidationError]:
        from adit.validate import check_hubbard

        m, errs = spec.method, []
        absent = [e for e in m.magnetization_by_element if e not in spec.elements]
        if absent:
            errs.append(ValidationError("method.magnetization_by_element", L(f"構造に無い元素です: {absent}", f"elements not in the structure: {absent}")))
        if not m.uks and any(v != 0 for v in m.magnetization_by_element.values()):
            errs.append(ValidationError("method.uks", L("元素ごとの MAGNETIZATION は UKS (スピン分極) のときだけ使われます", "per-element MAGNETIZATION is used only with UKS (spin polarization)")))
        errs += check_hubbard(m.hubbard, spec.elements, "method.hubbard")
        if m.sccs_relative_permittivity < 0:
            errs.append(ValidationError("method.sccs_relative_permittivity", L("負の値は指定できません (0 は溶媒なし)", "negative values are not allowed (0 means no solvent)")))
        return errs

    def version_probe(self, spec):
        return ("output.log", r"CP2K\| version string")

    def parameter_sources(self, spec: CalculationSpec, data: Cp2kData) -> dict[str, Path]:
        m = spec.method
        if not self._copies_entries(spec, data):
            return {}
        return {BASIS_OUT: data.resolve(m.basis_file), POTENTIAL_OUT: data.resolve(m.potential_file)}

    def _check_data(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m, st = spec.method, spec.structure
        data = _data(cfg)
        problem = _data_problem(cfg, data) if cfg.cp2k_data and not data.found else ""
        if problem:
            return [ValidationError("method.basis_file", problem)]
        errs: list[ValidationError] = []
        if not data.found and not (Path(m.basis_file).expanduser().is_absolute() and Path(m.potential_file).expanduser().is_absolute()):
            for kind, table in (("method.basis", m.basis), ("method.potential", m.potential)):
                missing = [e for e in spec.elements if not table.get(e, "").strip()]
                if missing:
                    what = L("基底関数名", "basis-set names") if kind == "method.basis" else L("擬ポテンシャル名", "pseudopotential names")
                    errs.append(ValidationError(kind, L(
                        f"この PC に CP2K の data ディレクトリが見つからないので、元素ごとの{what}を自動では選べません。{missing} の{what}を入れるか、"
                        "環境設定の cp2k_data に data ディレクトリ (BASIS_MOLOPT などのある場所) を書いてください",
                        f"no CP2K data directory was found on this PC, so {what} cannot be chosen automatically; enter {what} for {missing}, "
                        "or write the data directory (where BASIS_MOLOPT etc. are) in cp2k_data in the settings")))
            return errs
        for fkey, fname, potential in (("method.basis_file", m.basis_file, False), ("method.potential_file", m.potential_file, True)):
            if data.resolve(fname) is None:
                pref = "BASIS" if not potential else ""
                cands = [n for n in data.files(pref) if (not potential and "BASIS" in n) or (potential and "POTENTIAL" in n)]
                errs.append(ValidationError(fkey, L(f"{fname!r} が data ディレクトリ {data.root} にありません (ある例: {', '.join(cands[:12])})",
                                                    f"{fname!r} is not in the data directory {data.root} (available e.g.: {', '.join(cands[:12])})")))
        if errs:
            return errs
        zsum = 0
        for e in spec.elements:
            chosen = {}
            for kind, fname, potential in (("method.basis", m.basis_file, False), ("method.potential", m.potential_file, True)):
                name = _pick(m, data, e, potential=potential)
                cands = data.names_for(fname, e, potential=potential)
                if name is None:
                    msg = (L(f"{fname} に元素 {e} の項目がありません", f"{fname} has no entry for element {e}") if not cands else
                           L(f"元素 {e} の名前を選んでください ({fname} にある名前: {', '.join(cands)})", f"choose a name for element {e} (in {fname}: {', '.join(cands)})"))
                    errs.append(ValidationError(kind, msg))
                    continue
                entry = data.find(fname, e, name, potential=potential)
                if entry is None:
                    errs.append(ValidationError(kind, L(f"{fname} に元素 {e} の {name!r} がありません (ある名前: {', '.join(cands) or 'なし'})",
                                                        f"{fname} has no {name!r} for element {e} (available: {', '.join(cands) or 'none'})")))
                    continue
                chosen[potential] = entry
            if True in chosen and chosen[True].valence is None:
                errs.append(ValidationError("method.potential", L(f"{m.potential_file} の {chosen[True].names[0]} から価電子数を読めません (GTH の形ではない)",
                                                                  f"cannot read the valence of {chosen[True].names[0]} in {m.potential_file} (not GTH format)")))
            elif True in chosen:
                zsum += chosen[True].valence * st.atoms.symbols.count(e)
                qb = chosen[False].q if False in chosen else None
                if qb is not None and qb != chosen[True].valence:
                    errs.append(ValidationError("method.basis", L(
                        f"元素 {e}: 基底関数 {chosen[False].names[0]} は価電子 {qb} 個用 (-q{qb}) ですが、擬ポテンシャル {chosen[True].names[0]} の価電子は {chosen[True].valence} 個です",
                        f"element {e}: basis {chosen[False].names[0]} is for {qb} valence electrons (-q{qb}) but potential {chosen[True].names[0]} has {chosen[True].valence}")))
        if m.dispersion != "none" and data.found and not (data.root / D3_FILE).is_file():
            errs.append(ValidationError("method.dispersion", L(f"D3 の係数のファイル {D3_FILE} が {data.root} にありません", f"{D3_FILE} (D3 parameters) is not in {data.root}")))
        if not any(e.location in ("method.basis", "method.potential") for e in errs):
            pe = electron_parity_error(zsum - st.charge, st.charge, st.multiplicity)
            if pe:
                errs.append(pe)
        return errs

    def generate(self, spec: CalculationSpec, data: Cp2kData) -> dict[str, str]:
        out = {}
        if self._copies_entries(spec, data):
            out[BASIS_OUT] = self._subset(spec, data, potential=False)
            out[POTENTIAL_OUT] = self._subset(spec, data, potential=True)
        out[INPUT_FILE] = self.cp2k_inp(spec, data)
        return out

    @staticmethod
    def _copies_entries(spec: CalculationSpec, data: Cp2kData) -> bool:
        m = spec.method
        return data.resolve(m.basis_file) is not None and data.resolve(m.potential_file) is not None

    def _subset(self, spec: CalculationSpec, data: Cp2kData, *, potential: bool) -> str:
        m = spec.method
        fname = m.potential_file if potential else m.basis_file
        src = data.resolve(fname)
        lines = [f"# ADIT: {src} から、この計算で使う項目だけを写したもの (CP2K の配布物。GPL-2.0-or-later)", ""]
        for e in spec.elements:
            name = _pick(m, data, e, potential=potential)
            entry = data.find(fname, e, name, potential=potential) if name else None
            if entry is None:
                raise GenerationError(L(f"{fname} から元素 {e} の項目を決められません", f"cannot decide the entry for element {e} in {fname}"))
            lines += [entry.text, "#"]
        return "\n".join(lines) + "\n"

    def files_to_copy(self, spec: CalculationSpec, data: Cp2kData) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        cmd = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{cmd} -i {INPUT_FILE} > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, data: Cp2kData, copies) -> ReadmeNotes:
        m, t = spec.method, spec.task
        copied = self._copies_entries(spec, data)
        files = [L("  cp2k.inp      CP2K の入力 (&GLOBAL、&FORCE_EVAL (電子状態と構造)、&MOTION (最適化・MD) の節)",
                   "  cp2k.inp      CP2K input (&GLOBAL, &FORCE_EVAL (electronic structure and structure), &MOTION (optimization / MD))")]
        if copied:
            kinds = ", ".join(f"{e}: {_pick(m, data, e, potential=False)} / {_pick(m, data, e, potential=True)}" for e in spec.elements)
            files += [L(f"  BASIS_adit / POTENTIAL_adit   基底関数と擬ポテンシャル (内殻電子の効果をまとめたもの) のうち、この計算で使う項目だけを",
                        f"  BASIS_adit / POTENTIAL_adit   the basis-set and pseudopotential (stands in for the core electrons) entries used here, copied from"),
                      L(f"                {m.basis_file} と {m.potential_file} から写したもの ({kinds})",
                        f"                {m.basis_file} and {m.potential_file} ({kinds})")]
        prep = []
        if not copied:
            prep += [L(f"  基底関数と擬ポテンシャルは、実行する計算機の CP2K の data ディレクトリの {m.basis_file} と {m.potential_file} から読みます。",
                       f"  The basis sets and pseudopotentials are read from {m.basis_file} and {m.potential_file} in the CP2K data directory of the machine that runs it."),
                     L("  ADIT を動かした PC に data ディレクトリが無かったため、名前は未確認です。",
                       "  The PC that ran ADIT had no data directory, so the names were not checked (unverified).")]
        if m.dispersion != "none":
            prep.append(L(f"  D3 の係数は、実行する計算機の CP2K が自分の data ディレクトリの {D3_FILE} から読みます (このディレクトリには入っていません)。",
                          f"  The D3 parameters are read by CP2K from {D3_FILE} in its own data directory (not included here)."))
        out = [L("  output.log    CP2K の出力 ('ENERGY| Total FORCE_EVAL' の行が全エネルギー。単位は Hartree)",
                 "  output.log    CP2K output (the 'ENERGY| Total FORCE_EVAL' lines give the total energy, in hartree)")]
        out += {
            "geometry_optimization": [L(f"  {PROJECT}-pos-1.xyz  最適化の各ステップの構造 (最後が最適化後)。{PROJECT}-1.restart は続きから実行するためのもの",
                                        f"  {PROJECT}-pos-1.xyz  structure at each optimization step (the last is the result); {PROJECT}-1.restart is for continuing")],
            "molecular_dynamics": [L(f"  {PROJECT}-pos-1.xyz  MD の軌跡 (dump の間隔ごとの座標、xyz 形式)", f"  {PROJECT}-pos-1.xyz  MD trajectory (coordinates every dump interval, xyz format)"),
                                   L(f"  {PROJECT}-1.ener   各ステップの時刻・運動エネルギー・温度・ポテンシャルエネルギー・保存量",
                                     f"  {PROJECT}-1.ener   time, kinetic energy, temperature, potential energy and conserved quantity of each step")],
            "vibrations": [L("                振動数は output.log の VIB| Frequency の行 (cm⁻¹)。基準振動は adit-VIBRATIONS-1.mol (Molden 形式)",
                             "                frequencies are the VIB| Frequency lines of output.log (cm⁻¹); normal modes in adit-VIBRATIONS-1.mol (Molden format)")],
        }.get(t.type, [])
        if t.type == "molecular_dynamics" and t.md.ensemble == "NPT":
            out.append(L(f"  {PROJECT}-1.cell   各ステップのセル", f"  {PROJECT}-1.cell   cell at each step"))
        return ReadmeNotes(program=profile_program(), files=files, prepare=prep, outputs=out)

    # ---- cp2k.inp ----
    def cp2k_inp(self, spec: CalculationSpec, data: Cp2kData) -> str:
        m, st, t = spec.method, spec.structure, spec.task
        extra = {k.strip().upper(): v for k, v in m.extra_sections.items()}
        copied = self._copies_entries(spec, data)

        def more(path: str, indent: str) -> list[str]:
            text = extra.get(path, "").strip("\n")
            return [indent + l.strip() for l in text.splitlines() if l.strip()] if text else []

        run_type = {"single_point": "ENERGY_FORCE", "geometry_optimization": "CELL_OPT" if t.relax_cell != "no" else "GEO_OPT",
                    "molecular_dynamics": "MD", "vibrations": "VIBRATIONAL_ANALYSIS"}.get(t.type)
        if run_type is None:
            raise GenerationError(L(f"CP2K の生成器にない計算の種類です: {t.type}", f"calculation type not supported by the CP2K generator: {t.type}"))
        npt = t.type == "molecular_dynamics" and t.md.ensemble == "NPT"
        lines = ["# ADIT が生成した CP2K の入力 (キーワードは CP2K 2026.2 の入力の説明で確かめたもの)",
                 "&GLOBAL", f"  PROJECT {PROJECT}", f"  RUN_TYPE {run_type}", *more("GLOBAL", "  "), "&END GLOBAL"]
        lines += self._ext_restart(spec)
        lines += ["&FORCE_EVAL", "  METHOD QUICKSTEP"]
        if run_type == "CELL_OPT" or npt:
            lines.append("  STRESS_TENSOR ANALYTICAL")
        lines += ["  &DFT",
                  f"    BASIS_SET_FILE_NAME {BASIS_OUT if copied else m.basis_file}",
                  f"    POTENTIAL_FILE_NAME {POTENTIAL_OUT if copied else m.potential_file}",
                  f"    CHARGE {st.charge}", f"    MULTIPLICITY {st.multiplicity}"]
        if m.uks:
            lines.append("    UKS .TRUE.")
        if m.hubbard:
            lines.append(f"    PLUS_U_METHOD {m.plus_u_method}")
        if m.surface_dipole_correction:
            lines.append("    SURFACE_DIPOLE_CORRECTION .TRUE.")
        if m.surf_dip_dir:
            lines.append(f"    SURF_DIP_DIR {m.surf_dip_dir}")
        scf = ["    &SCF", f"      EPS_SCF {_f(m.eps_scf)}", f"      MAX_SCF {m.max_scf}"]
        if m.ot:
            scf += ["      &OT"]
            if m.ot_minimizer: scf.append(f"        MINIMIZER {m.ot_minimizer}")
            if m.ot_preconditioner: scf.append(f"        PRECONDITIONER {m.ot_preconditioner}")
            scf += ["        " + line.strip() for line in m.ot_extra.splitlines() if line.strip()]
            scf.append("      &END OT")
        lines += ["    &MGRID", f"      CUTOFF {_f(m.cutoff_ry)}", f"      REL_CUTOFF {_f(m.rel_cutoff_ry)}", "    &END MGRID",
                  *scf, *more("FORCE_EVAL/DFT/SCF", "      "), "    &END SCF",
                  "    &XC", f"      &XC_FUNCTIONAL {m.xc.strip()}", "      &END XC_FUNCTIONAL"]
        if m.dispersion != "none":
            lines += ["      &VDW_POTENTIAL", "        POTENTIAL_TYPE PAIR_POTENTIAL", "        &PAIR_POTENTIAL",
                      f"          TYPE {'DFTD3(BJ)' if m.dispersion == 'd3bj' else 'DFTD3'}", f"          PARAMETER_FILE_NAME {D3_FILE}",
                      f"          REFERENCE_FUNCTIONAL {m.xc.strip()}", "        &END PAIR_POTENTIAL", "      &END VDW_POTENTIAL"]
        lines += [*more("FORCE_EVAL/DFT/XC", "      "), "    &END XC"]
        kp = spec.kpoints
        if st.periodic and kp is not None and kp.mode != "gamma":
            n1, n2, n3 = kp.resolved_mesh(st.atoms.cell)
            lines += ["    &KPOINTS", f"      SCHEME MONKHORST-PACK {n1} {n2} {n3}", "    &END KPOINTS"]
        periodic = _periodic_letters(st.atoms.pbc)
        if periodic != "XYZ":
            lines += ["    &POISSON", f"      PERIODIC {periodic}", f"      POISSON_SOLVER {m.poisson_solver}", "    &END POISSON"]
        if m.sccs_relative_permittivity > 0:
            lines += ["    &SCCS", f"      RELATIVE_PERMITTIVITY {_f(m.sccs_relative_permittivity)}", "    &END SCCS"]
        lines += [*more("FORCE_EVAL/DFT", "    "), "  &END DFT", "  &SUBSYS", "    &CELL"]
        cell = np.asarray(st.atoms.cell, dtype=float)
        if not any(st.atoms.pbc) and abs(np.linalg.det(cell)) < 1e-6:
            cell = np.eye(3) * m.isolated_box_ang
        for name, v in zip("ABC", cell):
            lines.append(f"      {name} {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
        lines += [f"      PERIODIC {periodic}", "    &END CELL", "    &COORD"]
        lines += [f"      {s:2s} {x:16.10f} {y:16.10f} {z:16.10f}" for s, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        lines.append("    &END COORD")
        if periodic == "NONE":
            lines += ["    &TOPOLOGY", "      &CENTER_COORDINATES", "      &END CENTER_COORDINATES", "    &END TOPOLOGY"]
        for e in spec.elements:
            b = _pick(m, data, e, potential=False) or m.basis.get(e, "")
            p = _pick(m, data, e, potential=True) or m.potential.get(e, "")
            lines += [f"    &KIND {e}", f"      BASIS_SET {b}", f"      POTENTIAL {p}"]
            if m.magnetization_by_element.get(e):
                lines.append(f"      MAGNETIZATION {_f(m.magnetization_by_element[e])}")
            h = m.hubbard.get(e)
            if h is not None:
                from adit.validate import orbital_l
                lines += ["      &DFT_PLUS_U", f"        L {orbital_l(h.orbital)}", f"        U_MINUS_J [eV] {_f(h.u_ev - h.j_ev)}", "      &END DFT_PLUS_U"]
            lines.append("    &END KIND")
        lines += [*more("FORCE_EVAL/SUBSYS", "    "), "  &END SUBSYS",
                  "  &PRINT", "    &FORCES ON", "    &END FORCES", "  &END PRINT",
                  *more("FORCE_EVAL", "  "), "&END FORCE_EVAL"]
        motion = self._motion(spec, run_type, more)
        if motion:
            lines += motion
        if run_type == "VIBRATIONAL_ANALYSIS":
            lines += ["&VIBRATIONAL_ANALYSIS", "&END VIBRATIONAL_ANALYSIS"]
        return "\n".join(lines) + "\n"

    @staticmethod
    def _ext_restart(spec: CalculationSpec) -> list[str]:
        h = spec.handoff
        if h is None or "prev.restart" not in h.files:
            return []
        out = ["&EXT_RESTART", "  RESTART_FILE_NAME prev.restart", "  RESTART_DEFAULT .FALSE."]
        if h.at_run:
            out += ["  RESTART_POS .TRUE.", "  RESTART_CELL .TRUE."]
        if h.velocities:
            out.append("  RESTART_VEL .TRUE.")
        return out + ["&END EXT_RESTART"]

    def _motion(self, spec: CalculationSpec, run_type: str, more) -> list[str]:
        st, t = spec.structure, spec.task
        if run_type in ("ENERGY_FORCE", "VIBRATIONAL_ANALYSIS"):
            return []
        out = ["&MOTION"]
        if run_type in ("GEO_OPT", "CELL_OPT"):
            sec = run_type
            out += [f"  &{sec}", f"    OPTIMIZER {'LBFGS' if t.optimizer == 'LBFGS' else 'BFGS'}"]
            if t.max_steps > 0:
                out.append(f"    MAX_ITER {t.max_steps}")
            out.append(f"    MAX_FORCE {t.force_tolerance_ev_per_ang / HARTREE_PER_BOHR_IN_EV_PER_ANG:.6e}")  # eV/Å → Hartree/Bohr
            if sec == "CELL_OPT":
                out.append("    EXTERNAL_PRESSURE [bar] 0.0")
            out += [*more(f"MOTION/{sec}", "    "), f"  &END {sec}"]
        elif run_type == "MD":
            md = t.md
            ens = {"NVE": "NVE", "NVT": "LANGEVIN" if md.thermostat == "langevin" else "NVT", "NPT": "NPT_I"}[md.ensemble]
            out += ["  &MD", f"    ENSEMBLE {ens}", f"    STEPS {md.steps}", f"    TIMESTEP {_f(md.timestep_fs)}", f"    TEMPERATURE {_f(md.temperature_k)}"]
            if ens == "LANGEVIN":
                out += ["    &LANGEVIN", f"      GAMMA {_f(1.0 / md.coupling_time_fs)}", "    &END LANGEVIN"]
            elif ens in ("NVT", "NPT_I"):
                th = THERMOSTATS.get(md.thermostat)
                if th is None:
                    raise GenerationError(L(f"CP2K にはない熱浴です: {md.thermostat}", f"thermostat not available in CP2K: {md.thermostat}"))
                out += ["    &THERMOSTAT", f"      TYPE {th}", f"      &{th}", f"        TIMECON {_f(md.coupling_time_fs)}", f"      &END {th}", "    &END THERMOSTAT"]
            if ens == "NPT_I":
                out += ["    &BAROSTAT", f"      PRESSURE {_f(md.pressure_bar)}", f"      TIMECON {_f(md.barostat_time_fs)}", "    &END BAROSTAT"]
            out += [*more("MOTION/MD", "    "), "  &END MD"]
        fixed = self._fixed_blocks(spec)
        if fixed:
            out += ["  &CONSTRAINT", *fixed, "  &END CONSTRAINT"]
        every = t.md.dump_interval if run_type == "MD" else 1
        out += ["  &PRINT", "    &TRAJECTORY", "      FORMAT XYZ", "      &EACH", f"        {run_type if run_type != 'CELL_OPT' else 'CELL_OPT'} {every}",
                "      &END EACH", "    &END TRAJECTORY"]
        if run_type == "MD" and t.md.ensemble == "NPT":
            out += ["    &CELL ON", "      &EACH", f"        MD {every}", "      &END EACH", "    &END CELL"]
        out += ["  &END PRINT", *more("MOTION", "  "), "&END MOTION"]
        return out

    @staticmethod
    def _fixed_blocks(spec: CalculationSpec) -> list[str]:
        st = spec.structure
        groups: dict[str, list[int]] = {}
        for i in sorted(st.fixed_atoms):
            groups.setdefault("XYZ", []).append(i + 1)
        for k, move in sorted(st.fixed_axes.items(), key=lambda kv: int(kv[0])):
            if int(k) in st.fixed_atoms:
                continue
            comp = "".join(a for a, mv in zip("XYZ", move) if not mv)
            if comp:
                groups.setdefault(comp, []).append(int(k) + 1)
        out = []
        for comp, idx in groups.items():
            out += ["    &FIXED_ATOMS", f"      LIST {' '.join(str(i) for i in idx)}", f"      COMPONENTS_TO_FIX {comp}", "    &END FIXED_ATOMS"]
        return out


def profile_program() -> str:
    return "cp2k.psmp"


register(Cp2kGenerator())


# The review the CP2K FAQ names when the REFERENCES block of the output is not reproduced
CITATIONS = (
    Citation("cp2k_kuhne2020", r"""@article{cp2k_kuhne2020,
  author  = {K{\"u}hne, Thomas D. and Iannuzzi, Marcella and Del Ben, Mauro and Rybkin, Vladimir V. and Seewald, Patrick and Stein, Frederick and Laino, Teodoro and Khaliullin, Rustam Z. and Sch{\"u}tt, Ole and Schiffmann, Florian and Golze, Dorothea and Wilhelm, Jan and Chulkov, Sergey and Bani-Hashemian, Mohammad Hossein and Weber, Val{\'e}ry and Bor{\v{s}}tnik, Urban and Taillefumier, Mathieu and Jakobovits, Alice Shoshana and Lazzaro, Alfio and Pabst, Hans and M{\"u}ller, Tiziano and Schade, Robert and Guidon, Manuel and Andermatt, Samuel and Holmberg, Nico and Schenter, Gregory K. and Hehn, Anna and Bussy, Augustin and Belleflamme, Fabian and Tabacchi, Gloria and Gl{\"o}{\ss}, Andreas and Lass, Michael and Bethune, Iain and Mundy, Christopher J. and Plessl, Christian and Watkins, Matt and VandeVondele, Joost and Krack, Matthias and Hutter, J{\"u}rg},
  title   = {{CP2K}: An electronic structure and molecular dynamics software package -- {Quickstep}: Efficient and accurate electronic structure calculations},
  journal = {The Journal of Chemical Physics},
  volume  = {152},
  number  = {19},
  pages   = {194103},
  year    = {2020},
  doi     = {10.1063/5.0007045}
}""", doi="10.1063/5.0007045", source="https://www.cp2k.org/faq:cite"),
)
