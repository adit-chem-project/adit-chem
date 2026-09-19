
from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np
from ase.data import atomic_numbers

from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.plumed import PLUMED_FILE, PLUMED_LOG, output_lines as plumed_outputs, plumed_text, readme_lines as plumed_prepare
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, LammpsMethod
from adit.validate_types import ValidationError

INPUT_FILE = "in.lammps"
DATA_FILE = "data.lammps"
DUMP_FILE = "traj.lammpstrj"
FORCES_FILE = "forces.lammpstrj"
FINAL_DATA = "final.data"
DEFAULT_COMMAND = "mpirun -np {mpiprocs} lmp"
EV_PER_ANG_IN_KCAL_PER_MOL_ANG = 23.060547830619
BAR_PER_ATM = 1.01325
WRITABLE_STYLES = ("atomic", "charge")
MIN_STYLE = {"FIRE": "fire", "SteepestDescent": "sd"}
DATA_MARGIN_ANG = 1.0


def _f(x: float) -> str:
    return f"{x:.10g}"


def _data_header(path: Path) -> dict:
    info: dict = {"atoms": None, "types": None, "style": None, "pair_coeffs": False}
    with open(path, encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            s = line.split("#", 1)[0].strip()
            m = re.match(r"^(\d+)\s+atoms$", s)
            if m:
                info["atoms"] = int(m.group(1))
            m = re.match(r"^(\d+)\s+atom types$", s)
            if m:
                info["types"] = int(m.group(1))
            if s.startswith(("Pair Coeffs", "PairIJ Coeffs")):
                info["pair_coeffs"] = True
            if s.startswith("Atoms"):
                c = line.split("#", 1)
                info["style"] = c[1].strip().split()[0] if len(c) > 1 and c[1].strip() else None
                break
    return info


def _pair_coeff_lines(m: LammpsMethod) -> list[str]:
    out = []
    for l in m.pair_coeff.splitlines():
        l = l.strip()
        if l:
            out.append(l if l.split()[0] == "pair_coeff" else "pair_coeff " + l)
    return out


def type_elements(spec: CalculationSpec) -> list[str]:
    m = spec.method
    return list(m.type_elements) if m.type_elements else list(spec.elements)


class LammpsGenerator(InputGenerator):
    code = "lammps"
    uses_kpoints = False

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, LammpsMethod):
            return [ValidationError("method.code", L(f"LAMMPS の生成器に {m.code!r} の手法が渡されました", f"the LAMMPS generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        if not m.units:
            errs.append(ValidationError("method.units", L("単位系 (units) を選んでください。力場のファイルが前提とする単位系です (EAM は多くが metal、ReaxFF や分子の力場は多くが real。ファイルの説明で確かめます)",
                                                          "choose the unit system (units) the force field assumes (most EAM files: metal; ReaxFF and molecular force fields: real; check the file)")))
        if not m.pair_style.strip():
            errs.append(ValidationError("method.pair_style", L("pair_style の行を入れてください (例 eam、eam/alloy、reaxff NULL、mace no_domain_decomposition)",
                                                               "enter the pair_style line (e.g. eam, eam/alloy, reaxff NULL, mace no_domain_decomposition)")))
        elif "\n" in m.pair_style.strip():
            errs.append(ValidationError("method.pair_style", L("pair_style は 1 行で書いてください", "pair_style must be a single line")))
        if m.seed < 1:
            errs.append(ValidationError("method.seed", L("乱数の種は 1 以上の整数にしてください (LAMMPS の velocity create の決まり)", "the seed must be a positive integer (LAMMPS velocity create)")))
        if st.charge != 0 or st.multiplicity != 1:
            errs.append(ValidationError("structure.charge", L("LAMMPS では全電荷と多重度を使いません (電荷は data ファイルの原子ごとの値か、力場が決めます)。0 と 1 のままにしてください",
                                                              "LAMMPS does not use the total charge or multiplicity (charges come from the data file or the force field); leave them at 0 and 1")))
        if st.fixed_axes:
            errs.append(ValidationError("structure.fixed_axes", L("LAMMPS の生成器は軸ごとの固定に対応していません (原子ごとの固定は使えます)", "the LAMMPS generator cannot fix individual axes (whole atoms can be fixed)")))
        if 0 < sum(st.atoms.pbc) < 3:
            errs.append(ValidationError("structure.atoms", L("一部の方向だけ周期の構造は、LAMMPS の生成器では扱いません (3 方向とも周期か、分子にしてください)",
                                                             "structures periodic in only some directions are not handled by the LAMMPS generator (use fully periodic or a molecule)")))
        if t.type in ("vibrations", "band_structure"):
            errs.append(ValidationError("task.type", L("LAMMPS の生成器は一点計算・構造最適化 (最小化)・分子動力学だけです", "the LAMMPS generator supports single point, minimization and MD only")))
        if t.type == "geometry_optimization" and t.max_steps < 1:
            errs.append(ValidationError("task.max_steps", L("LAMMPS の minimize には最大反復回数が要ります (1 以上)", "LAMMPS minimize needs a maximum number of iterations (at least 1)")))
        if t.type == "molecular_dynamics":
            md = t.md
            if md.ensemble != "NVE" and md.thermostat == "andersen":
                errs.append(ValidationError("task.md.thermostat", L("LAMMPS に Andersen 熱浴はありません (nose_hoover / berendsen / csvr / langevin)",
                                                                    "LAMMPS has no Andersen thermostat (nose_hoover / berendsen / csvr / langevin)")))
            if md.ensemble == "NPT" and md.thermostat not in ("nose_hoover", "berendsen"):
                errs.append(ValidationError("task.md.thermostat", L("LAMMPS の NPT は nose_hoover (fix npt) か berendsen (fix temp/berendsen + press/berendsen) にしてください",
                                                                    "LAMMPS NPT: use nose_hoover (fix npt) or berendsen (fix temp/berendsen + press/berendsen)")))
        bad_el = [e for e in m.type_elements if e not in atomic_numbers]
        if bad_el:
            errs.append(ValidationError("method.type_elements", L(f"元素記号として読めません: {bad_el}", f"not chemical symbols: {bad_el}")))
        names = [Path(p).name for p in m.potential_files]
        for p in m.potential_files:
            if not Path(p).expanduser().is_file():
                errs.append(ValidationError("method.potential_files", L(f"ファイルがありません: {p}", f"file not found: {p}")))
        dup = sorted({n for n in names if names.count(n) > 1} | ({DATA_FILE, INPUT_FILE} & set(names)))
        if dup:
            errs.append(ValidationError("method.potential_files", L(f"同じ名前のファイルは 1 つのディレクトリに置けません: {dup}", f"files with the same name cannot share the directory: {dup}")))
        for l in [m.pair_style, *_pair_coeff_lines(m), *m.extra_commands.splitlines(), *m.style_commands.splitlines()]:
            for tok in l.split():
                if ("/" in tok or "\\" in tok) and Path(tok).expanduser().is_file():
                    errs.append(ValidationError("method.pair_coeff", L(f"{tok} はこの PC のパスです。ファイル名だけを書き、そのファイルを「写すファイル」に入れてください",
                                                                      f"{tok} is a path on this PC; write only the file name and add the file to the files to copy")))
        # data
        if m.data_file.strip():
            p = Path(m.data_file).expanduser()
            if not p.is_file():
                errs.append(ValidationError("method.data_file", L(f"data ファイルがありません: {p}", f"data file not found: {p}")))
            else:
                h = _data_header(p)
                if h["types"] is None or h["atoms"] is None:
                    errs.append(ValidationError("method.data_file", L(f"{p.name} の先頭に「N atoms」「N atom types」の行がありません (LAMMPS の data ファイルの形ではありません)",
                                                                      f"{p.name} lacks the 'N atoms' / 'N atom types' lines (not a LAMMPS data file)")))
                else:
                    if len(m.type_elements) != h["types"]:
                        errs.append(ValidationError("method.type_elements", L(
                            f"{p.name} の原子の型は {h['types']} 種類です。型番号 1〜{h['types']} の元素を順に入れてください (いまは {len(m.type_elements)} 個。軌跡に元素名を書くのに使います)",
                            f"{p.name} has {h['types']} atom types; give the elements of types 1..{h['types']} in order (now {len(m.type_elements)}; used for element names in the trajectory)")))
                    if h["atoms"] != len(st.atoms.symbols):
                        errs.append(ValidationError("structure.atoms", L(
                            f"構造の原子数 ({len(st.atoms.symbols)}) と {p.name} の原子数 ({h['atoms']}) が違います。構造にも同じ data ファイルを読み込んでください",
                            f"the structure has {len(st.atoms.symbols)} atoms but {p.name} has {h['atoms']}; load the same data file as the structure")))
                    if h["style"] and h["style"] != m.atom_style.split()[0]:
                        errs.append(ValidationError("method.atom_style", L(f"{p.name} の Atoms の節は {h['style']} 形式ですが、atom_style は {m.atom_style} です",
                                                                           f"the Atoms section of {p.name} is {h['style']} style, but atom_style is {m.atom_style}")))
                    if not _pair_coeff_lines(m) and not h["pair_coeffs"]:
                        errs.append(ValidationError("method.pair_coeff", L(f"pair_coeff の行がなく、{p.name} にも Pair Coeffs の節がありません", f"no pair_coeff lines and no Pair Coeffs section in {p.name}")))
        else:
            if m.atom_style not in WRITABLE_STYLES:
                errs.append(ValidationError("method.atom_style", L(
                    f"構造から data を書けるのは atomic と charge だけです ({m.atom_style} には結合などの情報が要ります。外部で作った data ファイルを指定してください)",
                    f"a data file can be written from the structure only for atomic and charge ({m.atom_style} needs bonds etc.; give a data file made elsewhere)")))
            if m.type_elements and sorted(m.type_elements) != sorted(spec.elements):
                errs.append(ValidationError("method.type_elements", L(f"型番号の元素 {m.type_elements} は、構造の元素 {spec.elements} を 1 回ずつ並べたものにしてください",
                                                                      f"the type elements {m.type_elements} must list each element of the structure {spec.elements} once")))
            if not _pair_coeff_lines(m):
                errs.append(ValidationError("method.pair_coeff", L("pair_coeff の行を入れてください (例 * * Cu_u3.eam)。構造から書く data には力場の係数が入りません",
                                                                   "enter the pair_coeff lines (e.g. * * Cu_u3.eam); the data file written from the structure has no force-field coefficients")))
        errs += self._check_handoff(spec)
        return errs

    @staticmethod
    def _check_handoff(spec: CalculationSpec) -> list[ValidationError]:
        h, m = spec.handoff, spec.method
        if h is None or not h.velocities or h.at_run:
            return []
        if not m.data_file.strip():
            return [ValidationError("method.data_file", L("前の MD の速度を使うには、前の計算の final.data を data ファイルに指定します",
                                                          "to use the velocities of the previous MD, give its final.data as the data file"))]
        p = Path(m.data_file).expanduser()
        if p.is_file():
            with open(p, encoding="utf-8", errors="replace") as f:
                if not any(line.split("#", 1)[0].strip() == "Velocities" for line in f):
                    return [ValidationError("method.data_file", L(f"{p.name} に Velocities の節がありません (速度を引き継げません)", f"{p.name} has no Velocities section (velocities cannot be carried over)"))]
        return []

    def version_probe(self, spec):
        return ("log.lammps", r"^LAMMPS \(")

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        out = {INPUT_FILE: self.in_lammps(spec)}
        if not spec.method.data_file.strip():
            out[DATA_FILE] = self.data(spec)
        text = plumed_text(spec)
        if text is not None:
            out[PLUMED_FILE] = text
        return out

    def files_to_copy(self, spec: CalculationSpec, res) -> dict[str, Path]:
        m = spec.method
        out = {Path(p).name: Path(p).expanduser() for p in m.potential_files}
        if m.data_file.strip():
            out[DATA_FILE] = Path(m.data_file).expanduser()
        return out

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        cmd = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{cmd} -in {INPUT_FILE} -log log.lammps > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res, copies) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        types = ", ".join(f"{i} = {e}" for i, e in enumerate(type_elements(spec), 1))
        u = {"metal": L("metal (時間 ps、エネルギー eV、長さ Å、圧力 bar)", "metal (time ps, energy eV, length Å, pressure bar)"),
             "real": L("real (時間 fs、エネルギー kcal/mol、長さ Å、圧力 atm)", "real (time fs, energy kcal/mol, length Å, pressure atm)")}.get(m.units, m.units)
        files = [L(f"  in.lammps     LAMMPS の入力 (単位系 {u})", f"  in.lammps     LAMMPS input (units {u})")]
        if m.data_file.strip():
            files.append(L(f"  data.lammps   構造と力場の情報 ({Path(m.data_file).name} を写したもの。そのまま使います)", f"  data.lammps   structure and force-field data (a copy of {Path(m.data_file).name}, used as is)"))
        else:
            files.append(L(f"  data.lammps   構造 ({m.atom_style} 形式。ASE で書いたもの。斜めのセルは LAMMPS の向きに回して書かれます)",
                           f"  data.lammps   structure ({m.atom_style} style, written with ASE; a skewed cell is rotated to the LAMMPS orientation)"))
        files.append(L(f"                原子の型番号と元素: {types}", f"                atom types and elements: {types}"))
        others = sorted(n for n in copies if n != DATA_FILE)
        if others:
            files.append(L(f"  {', '.join(others)}   力場・モデルのファイル (pair_coeff などから名前で読みます)", f"  {', '.join(others)}   force-field / model files (read by name from pair_coeff etc.)"))
        if spec.method.thermo_pressure_tensor:
            files.append(L("  (thermo に圧力テンソル pxy pxz pyz を足し、毎ステップ書きます。Green-Kubo で粘度を出すためで、log.lammps が大きくなります)",
                           "  (the thermo table includes pxy pxz pyz and is written every step, for Green-Kubo viscosity; log.lammps becomes large)"))
        out = [L("  output.log / log.lammps   LAMMPS の出力 (thermo の表: ステップ、時刻、温度、ポテンシャルエネルギー pe、運動エネルギー、全エネルギー、圧力、体積)",
                 "  output.log / log.lammps   LAMMPS output (thermo table: step, time, temperature, potential energy pe, kinetic and total energy, pressure, volume)")]
        out += {
            "single_point": [L(f"  {FORCES_FILE}  原子ごとの座標と力 (dump custom 形式)", f"  {FORCES_FILE}  per-atom positions and forces (dump custom format)")],
            "geometry_optimization": [L(f"  {DUMP_FILE}  最小化の途中の構造 (dump custom 形式、id type element x y z)。{FINAL_DATA} は最後の構造 (data 形式)",
                                        f"  {DUMP_FILE}  structures during minimization (dump custom, id type element x y z); {FINAL_DATA} is the final structure (data format)")],
            "molecular_dynamics": [L(f"  {DUMP_FILE}  MD の軌跡 (dump custom 形式、id type element x y z。OVITO や VMD でそのまま開けます)。{FINAL_DATA} は最後の構造",
                                     f"  {DUMP_FILE}  MD trajectory (dump custom, id type element x y z; opens as is in OVITO or VMD); {FINAL_DATA} is the final structure")],
        }.get(t, [])
        if spec.plumed is not None:
            files.append(L(f"  {PLUMED_FILE}   PLUMED の入力 (利用者が書いたもの)", f"  {PLUMED_FILE}   the PLUMED input (written by you)"))
            out += plumed_outputs(spec)
        return ReadmeNotes(program="lmp", files=files, prepare=plumed_prepare(spec), outputs=out)

    # ---- data ----
    def data(self, spec: CalculationSpec) -> str:
        from ase.io.lammpsdata import write_lammps_data
        m = spec.method
        atoms = spec.atoms
        if not any(atoms.pbc):
            pos = atoms.get_positions()
            lo = pos.min(axis=0) - DATA_MARGIN_ANG
            atoms.positions = pos - lo
            atoms.cell = np.diag(pos.max(axis=0) - pos.min(axis=0) + 2 * DATA_MARGIN_ANG)
        buf = io.StringIO()
        write_lammps_data(buf, atoms, specorder=type_elements(spec), masses=True, units=m.units or "metal", atom_style=m.atom_style)
        lines = buf.getvalue().splitlines()
        types = " ".join(f"{i}={e}" for i, e in enumerate(type_elements(spec), 1))
        lines[0] = f"adit: atom types {types} (written by ASE)"
        return "\n".join(lines) + "\n"

    # ---- in.lammps ----
    def in_lammps(self, spec: CalculationSpec) -> str:
        m, st, t = spec.method, spec.structure, spec.task
        tf = (lambda fs: fs / 1000.0) if m.units == "metal" else (lambda fs: fs)
        pf = (lambda bar: bar) if m.units == "metal" else (lambda bar: bar / BAR_PER_ATM)
        ff = 1.0 if m.units == "metal" else EV_PER_ANG_IN_KCAL_PER_MOL_ANG
        els = type_elements(spec)
        boundary = "p p p" if all(st.atoms.pbc) else "s s s"
        lines = ["# ADIT が生成した LAMMPS の入力 (コマンドは https://docs.lammps.org/ で確かめたもの)",
                 f"# 原子の型番号と元素: {', '.join(f'{i}={e}' for i, e in enumerate(els, 1))}",
                 f"units {m.units}", f"atom_style {m.atom_style}", f"boundary {boundary}"]
        lines += [l.strip() for l in m.style_commands.splitlines() if l.strip()]
        lines += [f"pair_style {m.pair_style.strip()}", f"read_data {DATA_FILE}", *_pair_coeff_lines(m)]
        lines += [l.strip() for l in m.extra_commands.splitlines() if l.strip()]
        mobile = "all"
        if st.fixed_atoms:
            lines += [f"group fixed id {' '.join(str(i + 1) for i in sorted(st.fixed_atoms))}", "group mobile subtract all fixed",
                      "fix freeze fixed setforce 0.0 0.0 0.0"]
            mobile = "mobile"
        every = t.md.dump_interval
        columns = "step time temp pe ke etotal press vol" + (" pxy pxz pyz" if m.thermo_pressure_tensor else "")
        interval = 1 if (m.thermo_pressure_tensor or t.type == "single_point") else every
        lines += [f"thermo_style custom {columns}", f"thermo {interval}"]
        dump = [f"dump traj all custom {every} {DUMP_FILE} id type element x y z", f"dump_modify traj element {' '.join(els)} sort id"]
        if t.type == "single_point":
            lines += [f"dump forces all custom 1 {FORCES_FILE} id type element x y z fx fy fz", f"dump_modify forces element {' '.join(els)} sort id", "run 0"]
            return "\n".join(lines) + "\n"
        if t.type == "geometry_optimization":
            if t.relax_cell != "no":
                cell = st.atoms.to_ase().cell
                ortho = np.allclose(cell.angles(), 90.0)
                kind = "iso" if t.relax_cell == "volume_only" else ("aniso" if ortho else "tri")
                lines.append(f"fix relax all box/relax {kind} 0.0")
            ftol = t.force_tolerance_ev_per_ang * ff
            lines += [*dump, f"min_style {MIN_STYLE.get(t.optimizer, 'cg')}",
                      f"minimize 0.0 {_f(ftol)} {t.max_steps} {10 * t.max_steps}", f"write_data {FINAL_DATA}"]
            return "\n".join(lines) + "\n"
        md = t.md
        T, tau, dt = _f(md.temperature_k), _f(tf(md.coupling_time_fs)), _f(tf(md.timestep_fs))
        P, ptau = _f(pf(md.pressure_bar)), _f(tf(md.barostat_time_fs))
        if spec.handoff is not None and spec.handoff.velocities:
            lines.append("# 速度は前の MD の write_data が data の Velocities の節に書いたものを read_data が読む (velocity create はしない)")
        else:
            lines.append(f"velocity {mobile} create {T} {m.seed} mom yes rot no dist gaussian")
        lines += [f"timestep {dt}", *dump]
        if md.ensemble == "NVE":
            lines.append(f"fix integ {mobile} nve")
        elif md.ensemble == "NVT":
            if md.thermostat == "nose_hoover":
                lines.append(f"fix integ {mobile} nvt temp {T} {T} {tau}")
            else:
                th = {"berendsen": f"temp/berendsen {T} {T} {tau}", "csvr": f"temp/csvr {T} {T} {tau} {m.seed}",
                      "langevin": f"langevin {T} {T} {tau} {m.seed}"}.get(md.thermostat)
                if th is None:
                    raise GenerationError(L(f"LAMMPS にない熱浴です: {md.thermostat}", f"thermostat not available in LAMMPS: {md.thermostat}"))
                lines += [f"fix integ {mobile} nve", f"fix tstat {mobile} {th}"]
        else:
            if md.thermostat == "nose_hoover":
                lines.append(f"fix integ {mobile} npt temp {T} {T} {tau} iso {P} {P} {ptau}")
            elif md.thermostat == "berendsen":
                lines += [f"fix integ {mobile} nve", f"fix tstat {mobile} temp/berendsen {T} {T} {tau}", f"fix pstat all press/berendsen iso {P} {P} {ptau}"]
            else:
                raise GenerationError(L(f"LAMMPS の NPT に対応していない熱浴です: {md.thermostat}", f"thermostat not supported for LAMMPS NPT: {md.thermostat}"))
        if spec.plumed is not None:
            lines.append(f"fix adit_plumed all plumed plumedfile {PLUMED_FILE} outfile {PLUMED_LOG}")
        lines += [f"run {md.steps}", f"write_data {FINAL_DATA}"]
        return "\n".join(lines) + "\n"


register(LammpsGenerator())


# The reference the LAMMPS manual asks for ("Citing LAMMPS")
CITATIONS = (
    Citation("lammps_thompson2022", r"""@article{lammps_thompson2022,
  author  = {Thompson, Aidan P. and Aktulga, H. Metin and Berger, Richard and Bolintineanu, Dan S. and Brown, W. Michael and Crozier, Paul S. and in 't Veld, Pieter J. and Kohlmeyer, Axel and Moore, Stan G. and Nguyen, Trung Dac and Shan, Ray and Stevens, Mark J. and Tranchida, Julien and Trott, Christian and Plimpton, Steven J.},
  title   = {{LAMMPS} - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales},
  journal = {Computer Physics Communications},
  volume  = {271},
  pages   = {108171},
  year    = {2022},
  doi     = {10.1016/j.cpc.2021.108171}
}""", doi="10.1016/j.cpc.2021.108171", source="https://docs.lammps.org/Intro_citing.html"),
)
