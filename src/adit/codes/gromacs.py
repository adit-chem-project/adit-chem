
from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.plumed import PLUMED_FILE, output_lines as plumed_outputs, plumed_text, readme_lines as plumed_prepare
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, GromacsMethod
from adit.validate_types import ValidationError

MDP_FILE = "grompp.mdp"
TOP_FILE = "topol.top"
CPT_FILE = "prev.cpt"
DEFFNM = "adit"
DEFAULT_COMMAND = "gmx"
KJ_PER_MOL_NM_PER_EV_ANG = 964.8533212  # 1 eV/Å = 96.485 kJ/mol ÷ 0.1 nm (CODATA 2018)
TCOUPL = {"berendsen": "berendsen", "nose_hoover": "nose-hoover", "andersen": "andersen", "csvr": "V-rescale"}
_INCLUDE = re.compile(r'^\s*#include\s+["<]([^">]+)[">]')
_IFDEF = re.compile(r"^\s*#(ifdef|ifndef)\s+(\S+)")


def _f(x: float) -> str:
    return f"{x:.10g}"


def gmx_top_dir(path: str | None = None) -> Path | None:
    env = os.environ.get("GMXLIB", "")
    if env and Path(env).is_dir():
        return Path(env)
    for exe in ("gmx", "gmx_mpi", "gmx_d", "gmx_mpi_d"):
        w = shutil.which(exe, path=path)
        if w:
            d = Path(w).resolve().parent.parent / "share" / "gromacs" / "top"
            if d.is_dir():
                return d
    return None


def scan_includes(top: Path, define: str = "", gmx_dir: Path | None = None) -> tuple[dict[str, Path], list[str], list[str]]:
    copies: dict[str, Path] = {}
    problems: list[str] = []
    from_lib: list[str] = []
    defined = set(re.findall(r"-D(\w+)", define))
    seen: set[Path] = set()

    def walk(f: Path, rel_dir: Path) -> None:
        if f in seen:
            return
        seen.add(f)
        stack: list[tuple[str, str]] = []
        for line in f.read_text(encoding="utf-8", errors="replace").splitlines():
            m = _IFDEF.match(line)
            if m:
                stack.append((m.group(1), m.group(2)))
                continue
            s = line.strip()
            if s.startswith("#endif") and stack:
                stack.pop()
                continue
            if s.startswith("#else") and stack:
                kind, name = stack[-1]
                stack[-1] = ("ifndef" if kind == "ifdef" else "ifdef", name)
                continue
            m = _INCLUDE.match(line)
            if not m:
                continue
            inc = m.group(1)
            active = all((name in defined) == (kind == "ifdef") for kind, name in stack)
            p = Path(inc)
            local = (f.parent / p) if not p.is_absolute() else p
            if local.is_file():
                if p.is_absolute() or ".." in p.parts:
                    problems.append(L(f"{f.name} の #include \"{inc}\" は、置き場所を変えると読めない書き方です (絶対パスか ..)。トポロジーと同じフォルダかその下に置き、相対の名前で書いてください",
                                      f"#include \"{inc}\" in {f.name} breaks when moved (absolute or ..); put it in or below the topology folder and use a relative name"))
                    continue
                rel = (rel_dir / p).as_posix()
                copies[rel] = local
                walk(local, rel_dir / p.parent)
            elif gmx_dir is not None and (gmx_dir / p).is_file():
                from_lib.append(inc)
            elif active:
                if gmx_dir is None:
                    from_lib.append(inc + L(" (未確認)", " (unverified)"))
                else:
                    problems.append(L(f"{f.name} の #include \"{inc}\" が、{f.parent} にも GROMACS の力場の置き場所 {gmx_dir} にもありません",
                                      f"#include \"{inc}\" in {f.name} is neither in {f.parent} nor in the GROMACS force-field folder {gmx_dir}"))

    walk(top, Path("."))
    return copies, problems, from_lib


def _gro_atom_count(p: Path) -> int | None:
    try:
        if p.suffix.lower() == ".gro":
            with open(p, encoding="utf-8", errors="replace") as f:
                f.readline()
                return int(f.readline().split()[0])
        with open(p, encoding="utf-8", errors="replace") as f:
            return sum(1 for l in f if l.startswith(("ATOM", "HETATM")))
    except (OSError, ValueError, IndexError):
        return None


def _gro_box_nm(p: Path) -> tuple[float, float, float] | None:
    if p.suffix.lower() != ".gro":
        return None
    try:
        with open(p, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 400))
            last = f.read().decode("utf-8", errors="replace").strip().splitlines()[-1].split()
        return (float(last[0]), float(last[1]), float(last[2]))
    except (OSError, ValueError, IndexError):
        return None


def structure_name(m: GromacsMethod) -> str:
    return "conf" + (Path(m.structure_file).suffix.lower() or ".gro")


def _run_conf(spec: CalculationSpec) -> str:
    h = spec.handoff
    return "conf.gro" if h is not None and h.at_run and "conf.gro" in h.files else structure_name(spec.method)


def _continues(spec: CalculationSpec) -> bool:
    h = spec.handoff
    return bool(spec.method.checkpoint_file.strip()) or (h is not None and h.at_run and h.velocities)


class GromacsGenerator(InputGenerator):
    code = "gromacs"
    uses_kpoints = False

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return gmx_top_dir()

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, GromacsMethod):
            return [ValidationError("method.code", L(f"GROMACS の生成器に {m.code!r} の手法が渡されました", f"the GROMACS generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        top = Path(m.topology_file).expanduser() if m.topology_file.strip() else None
        if top is None:
            errs.append(ValidationError("method.topology_file", L("トポロジー (.top) を指定してください (CHARMM-GUI、acpype、pdb2gmx などで作ったもの)",
                                                                  "give the topology (.top) made with CHARMM-GUI, acpype, pdb2gmx, etc.")))
        elif not top.is_file():
            errs.append(ValidationError("method.topology_file", L(f"ファイルがありません: {top}", f"file not found: {top}")))
        else:
            _, problems, _ = scan_includes(top, m.define, gmx_top_dir())
            errs += [ValidationError("method.topology_file", p) for p in problems]
        conf = Path(m.structure_file).expanduser() if m.structure_file.strip() else None
        if conf is None:
            errs.append(ValidationError("method.structure_file", L("構造のファイル (.gro か .pdb) を指定してください", "give the structure file (.gro or .pdb)")))
        elif not conf.is_file():
            errs.append(ValidationError("method.structure_file", L(f"ファイルがありません: {conf}", f"file not found: {conf}")))
        elif conf.suffix.lower() not in (".gro", ".pdb"):
            errs.append(ValidationError("method.structure_file", L("構造のファイルは .gro か .pdb にしてください", "the structure file must be .gro or .pdb")))
        else:
            box = _gro_box_nm(conf)
            rc = max(m.rcoulomb_nm, m.rvdw_nm)
            rl = rc * (1.05 if t.type != "molecular_dynamics" else 1.0)
            if box is not None and rl >= min(box) / 2:
                errs.append(ValidationError("method.rcoulomb_nm", L(
                    f"近傍リストの半径 {rl:.4g} nm (カットオフ {rc:g} nm{' の 1.05 倍' if rl != rc else ''}) が箱の短い辺 ({min(box):g} nm) の半分以上です。"
                    "grompp が止まります (カットオフを短くするか、箱を大きくしてください)",
                    f"the neighbor-list radius {rl:.4g} nm (cut-off {rc:g} nm{' x 1.05' if rl != rc else ''}) is at least half the shortest box edge "
                    f"({min(box):g} nm); grompp stops (shorten the cut-off or enlarge the box)")))
            n = _gro_atom_count(conf)
            if n is not None and n != len(st.atoms.symbols):
                errs.append(ValidationError("structure.atoms", L(
                    f"構造の原子数 ({len(st.atoms.symbols)}) と {conf.name} の原子数 ({n}) が違います。構造にも同じファイルを読み込んでください",
                    f"the structure has {len(st.atoms.symbols)} atoms but {conf.name} has {n}; load the same file as the structure")))
        if m.checkpoint_file.strip() and not Path(m.checkpoint_file).expanduser().is_file():
            errs.append(ValidationError("method.checkpoint_file", L(f"ファイルがありません: {m.checkpoint_file}", f"file not found: {m.checkpoint_file}")))
        if m.rcoulomb_nm <= 0 or m.rvdw_nm <= 0:
            errs.append(ValidationError("method.rcoulomb_nm", L("カットオフ [nm] は 0 より大きい値が必要です", "the cut-offs [nm] must be greater than 0")))
        if st.charge != 0 or st.multiplicity != 1:
            errs.append(ValidationError("structure.charge", L("GROMACS では全電荷と多重度を使いません (電荷はトポロジーが決めます)。0 と 1 のままにしてください",
                                                              "GROMACS does not use the total charge or multiplicity (the topology sets the charges); leave them at 0 and 1")))
        if st.fixed_atoms or st.fixed_axes:
            errs.append(ValidationError("structure.fixed_atoms", L("GROMACS の固定 (freezegrps) は index ファイルの群で指定するので、この生成器では扱いません。位置の拘束 (define = -DPOSRES) か extra_mdp を使ってください",
                                                                   "GROMACS freezing (freezegrps) needs index groups, which this generator does not handle; use position restraints (define = -DPOSRES) or extra_mdp")))
        bad = [k for k in m.extra_mdp if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", k.strip())]
        if bad:
            errs.append(ValidationError("method.extra_mdp", L(f"mdp の項目の名前として読めません: {bad}", f"not valid mdp option names: {bad}")))
        if t.type in ("vibrations", "band_structure"):
            errs.append(ValidationError("task.type", L("GROMACS の生成器は一点計算・エネルギー最小化・分子動力学だけです (基準振動は倍精度の GROMACS と integrator = nm が要ります)",
                                                       "the GROMACS generator supports single point, energy minimization and MD only (normal modes need double-precision GROMACS and integrator = nm)")))
        if t.type == "geometry_optimization":
            if t.relax_cell != "no":
                errs.append(ValidationError("task.relax_cell", L("GROMACS のエネルギー最小化は箱を動かしません (NPT の MD で密度を合わせます)", "GROMACS energy minimization does not change the box (use NPT MD)")))
            if t.max_steps < 1:
                errs.append(ValidationError("task.max_steps", L("1 以上が必要です (mdp の nsteps)", "must be at least 1 (mdp nsteps)")))
        if t.type == "molecular_dynamics" and t.md.ensemble != "NVE":
            if t.md.thermostat == "berendsen":
                errs.append(ValidationError("task.md.thermostat", L("GROMACS の grompp は Berendsen 熱浴に警告を出して止まります (csvr = V-rescale / nose_hoover / langevin = sd から選んでください)",
                                                                    "GROMACS grompp stops with a warning for the Berendsen thermostat (choose csvr = V-rescale / nose_hoover / langevin = sd)")))
            if t.md.thermostat == "andersen":
                errs.append(ValidationError("task.md.thermostat", L("GROMACS の Andersen 熱浴は integrator = md-vv だけで使え、この生成器 (integrator = md) では grompp が止まります",
                                                                    "GROMACS Andersen needs integrator = md-vv; with this generator (integrator = md) grompp stops")))
        if t.type == "molecular_dynamics" and t.md.ensemble == "NPT":
            if m.pcoupl == "Parrinello-Rahman" and not _continues(spec):
                errs.append(ValidationError("method.pcoupl", L(
                    "Parrinello-Rahman で初速を作る (前の段階の .cpt が無い) と、grompp が「平衡化には不安定」と警告して止まります。前の段階 (C-rescale の NPT など) の .cpt を指定してください",
                    "Parrinello-Rahman with generated velocities (no previous .cpt) makes grompp stop with an 'unstable for equilibration' warning; give the .cpt of a previous stage (e.g. C-rescale NPT)")))
            if not m.pcoupl:
                errs.append(ValidationError("method.pcoupl", L("NPT では圧力浴 (C-rescale か Parrinello-Rahman) を選んでください", "choose a barostat (C-rescale or Parrinello-Rahman) for NPT")))
            if m.compressibility_per_bar <= 0:
                errs.append(ValidationError("method.compressibility_per_bar", L("NPT では等方圧縮率 [1/bar] が要ります (GROMACS のマニュアルの水の例は 4.5e-5)",
                                                                                "NPT needs the isothermal compressibility in 1/bar (the GROMACS manual uses 4.5e-5 for water)")))
        return errs

    def version_probe(self, spec):
        return ("adit.log", "GROMACS version:")

    def generate(self, spec: CalculationSpec, gmx_dir) -> dict[str, str]:
        out = {MDP_FILE: self.mdp(spec)}
        text = plumed_text(spec)
        if text is not None:
            out[PLUMED_FILE] = text
        return out

    def files_to_copy(self, spec: CalculationSpec, gmx_dir) -> dict[str, Path]:
        m = spec.method
        top = Path(m.topology_file).expanduser()
        copies, problems, _ = scan_includes(top, m.define, gmx_dir)
        if problems:
            raise GenerationError("; ".join(problems))
        out = {TOP_FILE: top, structure_name(m): Path(m.structure_file).expanduser(), **copies}
        if m.checkpoint_file.strip():
            out[CPT_FILE] = Path(m.checkpoint_file).expanduser()
        return out

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        m, r = spec.method, spec.runtime
        cmd = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=r.mpiprocs, omp_threads=r.omp_threads, binary="")
        exe = cmd.split()[-1]
        conf = _run_conf(spec)
        grompp = [exe, "grompp", "-f", MDP_FILE, "-c", conf, "-p", TOP_FILE, "-o", f"{DEFFNM}.tpr"]
        if m.define.strip():
            grompp += ["-r", conf]
        if _continues(spec):
            grompp += ["-t", CPT_FILE]
        mdrun = [cmd, "mdrun", "-deffnm", DEFFNM]
        if Path(exe).name == "gmx":
            mdrun += ["-ntmpi", str(r.mpiprocs)]
        mdrun += ["-ntomp", str(r.omp_threads)]
        if spec.plumed is not None:
            mdrun += ["-plumed", PLUMED_FILE]
        return " ".join(grompp) + " > grompp.log 2>&1 && " + " ".join(mdrun) + " > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, gmx_dir, copies) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        top = Path(m.topology_file).expanduser()
        _, _, from_lib = scan_includes(top, m.define, gmx_dir) if top.is_file() else ({}, [], [])
        incs = sorted(n for n in copies if n not in (TOP_FILE, structure_name(m), CPT_FILE))
        files = [L("  grompp.mdp    計算の設定 (mdp。gmx grompp が読みます)", "  grompp.mdp    run settings (mdp; read by gmx grompp)"),
                 L(f"  topol.top     トポロジー ({top.name} を写したもの。分子の種類と数、力場)", f"  topol.top     topology (a copy of {top.name}; molecule types and counts, force field)"),
                 L(f"  {structure_name(m)}     構造 ({Path(m.structure_file).name} を写したもの)", f"  {structure_name(m)}     structure (a copy of {Path(m.structure_file).name})")]
        if incs:
            files.append(L(f"  {', '.join(incs)}   トポロジーが #include で読むファイル", f"  {', '.join(incs)}   files read by #include from the topology"))
        if from_lib:
            files.append(L(f"                GROMACS に付いている力場から読むもの: {', '.join(from_lib)}", f"                read from the force fields bundled with GROMACS: {', '.join(from_lib)}"))
        if m.checkpoint_file.strip():
            files.append(L(f"  prev.cpt      前の段階のチェックポイント ({Path(m.checkpoint_file).name} を写したもの。速度と圧力浴・熱浴の状態を引き継ぎます)",
                           f"  prev.cpt      state of the previous stage (a copy of {Path(m.checkpoint_file).name}; carries over velocities and the thermostat and barostat state)"))
        files.append(L("                submit.sh が gmx grompp (設定とトポロジーをまとめて adit.tpr にする) → gmx mdrun (計算) の順に実行します",
                       "                submit.sh runs gmx grompp (combines settings and topology into adit.tpr), then gmx mdrun (the calculation)"))
        out = [L("  grompp.log    gmx grompp の出力 (トポロジーの誤りや注意はここに出ます)", "  grompp.log    output of gmx grompp (topology errors and notes appear here)"),
               L("  output.log / adit.log   gmx mdrun の出力 (エネルギーの表)", "  output.log / adit.log   output of gmx mdrun (energy table)"),
               L("  adit.edr     エネルギーの記録。gmx energy -f adit.edr で温度・圧力・密度などを取り出せます",
                 "  adit.edr     energy file; gmx energy -f adit.edr extracts temperature, pressure, density, ..."),
               L("  adit.gro     最後の構造", "  adit.gro     final structure")]
        if t == "molecular_dynamics":
            out += [L("  adit.xtc     MD の軌跡 (dump の間隔ごと。VMD や MDAnalysis で開けます)", "  adit.xtc     MD trajectory (every dump interval; opens in VMD or MDAnalysis)"),
                    L("  adit.cpt     チェックポイント (続きから実行するための状態)。次の段階 (NVT → NPT → 本計算) を ADIT で作るとき「前の段階の .cpt」に指定します",
                      "  adit.cpt     checkpoint; give it as the previous-stage .cpt when making the next stage (NVT → NPT → production) in ADIT")]
        prepare = []
        if t == "molecular_dynamics":
            dt = spec.task.md.timestep_fs
            prepare = [L(f"  この入力は時間刻み dt = {dt:g} fs、constraints = {m.constraints} を使います。時間刻みと拘束の組み合わせは、使用する力場の文書に照らして利用者が決めます。",
                         f"  This input uses dt = {dt:g} fs and constraints = {m.constraints}. Choose the timestep–constraint combination according to the documentation for your force field.")]
            if not m.checkpoint_file.strip():
                prepare.append(L(f"  初速を作る乱数種 gen-seed = {m.gen_seed} は grompp.mdp と spec.json の両方に記録されます。",
                                 f"  The random seed for initial velocities, gen-seed = {m.gen_seed}, is recorded in both grompp.mdp and spec.json."))
        if spec.plumed is not None:
            files.append(L(f"  {PLUMED_FILE}    PLUMED の入力 (利用者が書いたもの)", f"  {PLUMED_FILE}    the PLUMED input (written by you)"))
            prepare += plumed_prepare(spec)
            out += plumed_outputs(spec)
        return ReadmeNotes(program="gmx", files=files, prepare=prepare, outputs=out)

    # ---- mdp ----
    def mdp(self, spec: CalculationSpec) -> str:
        m, t = spec.method, spec.task
        lines = ["; ADIT が生成した GROMACS の mdp (項目は https://manual.gromacs.org/current/user-guide/mdp-options.html で確かめたもの)"]
        opts: dict[str, str] = {}
        if t.type == "geometry_optimization":
            opts["integrator"] = "l-bfgs" if t.optimizer == "LBFGS" else "steep"
            opts["nsteps"] = str(t.max_steps)
            opts["emtol"] = _f(t.force_tolerance_ev_per_ang * KJ_PER_MOL_NM_PER_EV_ANG)  # eV/Å → kJ/mol/nm
        elif t.type == "single_point":
            opts["integrator"] = "md"
            opts["nsteps"] = "0"
        elif t.type == "molecular_dynamics":
            md = t.md
            opts["integrator"] = "sd" if md.ensemble != "NVE" and md.thermostat == "langevin" else "md"
            opts["nsteps"] = str(md.steps)
            opts["dt"] = _f(md.timestep_fs / 1000.0)  # fs → ps
            opts["nstxout-compressed"] = str(md.dump_interval)
        else:
            raise GenerationError(L(f"GROMACS の生成器にない計算の種類です: {t.type}", f"calculation type not supported by the GROMACS generator: {t.type}"))
        every = str(t.md.dump_interval) if t.type == "molecular_dynamics" else "1"
        opts["nstenergy"] = every
        opts["nstlog"] = every
        opts["cutoff-scheme"] = "Verlet"
        opts["coulombtype"] = m.coulombtype
        opts["rcoulomb"] = _f(m.rcoulomb_nm)
        opts["rvdw"] = _f(m.rvdw_nm)
        opts["constraints"] = m.constraints
        if m.define.strip():
            opts["define"] = m.define.strip()
        if t.type == "molecular_dynamics":
            md = t.md
            if md.ensemble != "NVE":
                if opts["integrator"] != "sd":
                    th = TCOUPL.get(md.thermostat)
                    if th is None:
                        raise GenerationError(L(f"GROMACS にない熱浴です: {md.thermostat}", f"thermostat not available in GROMACS: {md.thermostat}"))
                    opts["tcoupl"] = th
                opts["tc-grps"] = "System"
                opts["tau-t"] = _f(md.coupling_time_fs / 1000.0)
                opts["ref-t"] = _f(md.temperature_k)
            else:
                opts["tcoupl"] = "no"
            if md.ensemble == "NPT":
                opts["pcoupl"] = m.pcoupl
                opts["pcoupltype"] = "isotropic"
                opts["tau-p"] = _f(md.barostat_time_fs / 1000.0)
                opts["ref-p"] = _f(md.pressure_bar)
                opts["compressibility"] = _f(m.compressibility_per_bar)
            else:
                opts["pcoupl"] = "no"
            if _continues(spec):
                opts["continuation"] = "yes"
                opts["gen-vel"] = "no"
            else:
                opts["gen-vel"] = "yes"
                opts["gen-temp"] = _f(md.temperature_k)
                opts["gen-seed"] = str(m.gen_seed)
        for k, v in m.extra_mdp.items():
            opts[k.strip()] = ("yes" if v else "no") if isinstance(v, bool) else str(v)
        width = max(len(k) for k in opts)
        lines += [f"{k:<{width}} = {v}" for k, v in opts.items()]
        return "\n".join(lines) + "\n"


register(GromacsGenerator())


# The papers the GROMACS reference manual asks for ("Citation information"; Bekker 1993 has no DOI and is omitted)
_GMX_CITE_URL = "https://manual.gromacs.org/current/reference-manual/preface.html"
CITATIONS = (
    Citation("gromacs_abraham2015", r"""@article{gromacs_abraham2015,
  author  = {Abraham, Mark James and Murtola, Teemu and Schulz, Roland and P{\'a}ll, Szil{\'a}rd and Smith, Jeremy C. and Hess, Berk and Lindahl, Erik},
  title   = {{GROMACS}: High performance molecular simulations through multi-level parallelism from laptops to supercomputers},
  journal = {SoftwareX},
  volume  = {1-2},
  pages   = {19--25},
  year    = {2015},
  doi     = {10.1016/j.softx.2015.06.001}
}""", doi="10.1016/j.softx.2015.06.001", source=_GMX_CITE_URL),
    Citation("gromacs_pall2015", r"""@incollection{gromacs_pall2015,
  author    = {P{\'a}ll, Szil{\'a}rd and Abraham, Mark James and Kutzner, Carsten and Hess, Berk and Lindahl, Erik},
  title     = {Tackling Exascale Software Challenges in Molecular Dynamics Simulations with {GROMACS}},
  booktitle = {Solving Software Challenges for Exascale},
  editor    = {Markidis, Stefano and Laure, Erwin},
  publisher = {Springer International Publishing},
  pages     = {3--27},
  year      = {2015},
  doi       = {10.1007/978-3-319-15976-8_1}
}""", doi="10.1007/978-3-319-15976-8_1", source=_GMX_CITE_URL),
    Citation("gromacs_pronk2013", r"""@article{gromacs_pronk2013,
  author  = {Pronk, Sander and P{\'a}ll, Szil{\'a}rd and Schulz, Roland and Larsson, Per and Bjelkmar, P{\"a}r and Apostolov, Rossen and Shirts, Michael R. and Smith, Jeremy C. and Kasson, Peter M. and van der Spoel, David and Hess, Berk and Lindahl, Erik},
  title   = {{GROMACS} 4.5: a high-throughput and highly parallel open source molecular simulation toolkit},
  journal = {Bioinformatics},
  volume  = {29},
  number  = {7},
  pages   = {845--854},
  year    = {2013},
  doi     = {10.1093/bioinformatics/btt055}
}""", doi="10.1093/bioinformatics/btt055", source=_GMX_CITE_URL),
    Citation("gromacs_hess2008", r"""@article{gromacs_hess2008,
  author  = {Hess, Berk and Kutzner, Carsten and van der Spoel, David and Lindahl, Erik},
  title   = {{GROMACS} 4: Algorithms for Highly Efficient, Load-Balanced, and Scalable Molecular Simulation},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {4},
  number  = {3},
  pages   = {435--447},
  year    = {2008},
  doi     = {10.1021/ct700301q}
}""", doi="10.1021/ct700301q", source=_GMX_CITE_URL),
    Citation("gromacs_vanderspoel2005", r"""@article{gromacs_vanderspoel2005,
  author  = {van der Spoel, David and Lindahl, Erik and Hess, Berk and Groenhof, Gerrit and Mark, Alan E. and Berendsen, Herman J. C.},
  title   = {{GROMACS}: Fast, flexible, and free},
  journal = {Journal of Computational Chemistry},
  volume  = {26},
  number  = {16},
  pages   = {1701--1718},
  year    = {2005},
  doi     = {10.1002/jcc.20291}
}""", doi="10.1002/jcc.20291", source=_GMX_CITE_URL),
    Citation("gromacs_lindahl2001", r"""@article{gromacs_lindahl2001,
  author  = {Lindahl, Erik and Hess, Berk and van der Spoel, David},
  title   = {{GROMACS} 3.0: a package for molecular simulation and trajectory analysis},
  journal = {Journal of Molecular Modeling},
  volume  = {7},
  number  = {8},
  pages   = {306--317},
  year    = {2001},
  doi     = {10.1007/s008940100045}
}""", doi="10.1007/s008940100045", source=_GMX_CITE_URL),
    Citation("gromacs_berendsen1995", r"""@article{gromacs_berendsen1995,
  author  = {Berendsen, H. J. C. and van der Spoel, D. and van Drunen, R.},
  title   = {{GROMACS}: A message-passing parallel molecular dynamics implementation},
  journal = {Computer Physics Communications},
  volume  = {91},
  number  = {1-3},
  pages   = {43--56},
  year    = {1995},
  doi     = {10.1016/0010-4655(95)00042-E}
}""", doi="10.1016/0010-4655(95)00042-E", source=_GMX_CITE_URL),
)
