"""Scan the output logs of a run for the warnings and failure messages each code is known to print."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from adit.lang import L

SCF = "scf_not_converged"
GEOMETRY = "geometry_not_converged"
WALLTIME = "walltime"
MEMORY = "out_of_memory"
MPI = "mpi_abort"
DIAGONALIZATION = "diagonalization"
LOST_ATOMS = "lost_atoms"
UNSTABLE = "unstable_dynamics"
INPUT = "input_error"
RESOURCE = "resource_limit"
SIGNAL = "killed_by_signal"
ERROR = "error_message"

MAX_STORED_HITS = 30
EXCERPT_LENGTH = 120


@dataclass(frozen=True)
class Pattern:
    text: str
    cause: str | None = None
    regex: bool = False
    warning: bool = False

    def matches(self, line: str) -> bool:
        return re.search(self.text, line) is not None if self.regex else self.text in line


# Message strings were taken from the programs themselves (binaries or sources) or from their
# documentation; see remedy_lines() for the references. Scheduler, MPI and OS messages apply to every code.
COMMON: tuple[Pattern, ...] = (
    Pattern("DUE TO TIME LIMIT", WALLTIME),                                   # slurmstepd
    Pattern("=>> PBS: job killed: walltime", WALLTIME),                       # OpenPBS pbs_mom
    Pattern("=>> PBS: job killed: mem", MEMORY),
    Pattern("=>> PBS: job killed: vmem", MEMORY),
    Pattern("=>> PBS: job killed:", RESOURCE),
    Pattern("oom_kill event", MEMORY),                                        # slurmstepd (cgroup)
    Pattern("OOM Killed", MEMORY),
    Pattern("Out of memory: Killed process", MEMORY),                         # Linux kernel
    Pattern(r"(^|\s)Killed\s{2,}\S", SIGNAL, regex=True),                     # shell report of SIGKILL
    Pattern("Segmentation fault", SIGNAL),
    Pattern("MPI_ABORT was invoked", MPI),                                    # Open MPI
    Pattern("MPI_ERRORS_ARE_FATAL", MPI),
    Pattern("application called MPI_Abort", MPI),                             # MPICH
    Pattern("BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES", MPI),     # Hydra
)

PATTERNS: dict[str, tuple[Pattern, ...]] = {
    "dftbplus": (
        Pattern("SCC is NOT converged, maximal SCC iterations exceeded", SCF),
        Pattern("Geometry did NOT converge", GEOMETRY),
        Pattern("ERROR!", ERROR),
    ),
    "xtb": (
        Pattern("Self consistent charge iterator did not converge", SCF),
        Pattern("convergence criteria cannot be satisfied within", SCF, warning=True),
        Pattern("FAILED TO CONVERGE GEOMETRY OPTIMIZATION", GEOMETRY),
        Pattern("Optimization did not converge, aborting", GEOMETRY),
        Pattern("abnormal termination of xtb", ERROR),
        Pattern("[ERROR]", ERROR),
    ),
    "espresso": (
        Pattern("convergence NOT achieved", SCF),
        Pattern("Maximum CPU time exceeded", WALLTIME),
        Pattern("eigenvalues not converged", DIAGONALIZATION, warning=True),
        Pattern("too many bands are not converged", DIAGONALIZATION),
        Pattern("eigenvectors failed to converge", DIAGONALIZATION),
        Pattern("S matrix not positive definite", DIAGONALIZATION),
        Pattern("problems computing cholesky", DIAGONALIZATION),
        Pattern("zhegvd failed", DIAGONALIZATION),
        Pattern("bfgs failed", GEOMETRY),
        Pattern("history already reset at previous step", GEOMETRY),
        Pattern("Error in routine broyden", DIAGONALIZATION),
        Pattern("Not enough space allocated for radial FFT", ERROR),
        Pattern("Error in routine", ERROR),
    ),
    "vasp": (
        Pattern("ZBRENT: fatal error in bracketing", GEOMETRY),
        Pattern("EDDDAV: Call to ZHEGV failed", DIAGONALIZATION),
        Pattern("Sub-Space-Matrix is not hermitian", DIAGONALIZATION, warning=True),
        Pattern("EDDRMM: call to ZHEGV failed", DIAGONALIZATION, warning=True),
    ),
    "orca": (
        Pattern("OUT OF MEMORY ERROR", MEMORY),
        Pattern("Please increase MaxCore", MEMORY),
        Pattern("INPUT ERROR", INPUT),
    ),
    "cp2k": (
        Pattern("SCF run NOT converged", SCF),
        Pattern("MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED", GEOMETRY),
        Pattern("CPASSERT failed", ERROR),
        Pattern("[ABORT]", ERROR),
    ),
    "lammps": (
        Pattern("Lost atoms:", LOST_ATOMS),
        Pattern("Out of range atoms - cannot compute", UNSTABLE),
        Pattern("Bond atoms missing at step", UNSTABLE),
        Pattern("Non-numeric", UNSTABLE),
        Pattern("ERROR on proc", ERROR),
        Pattern("ERROR:", ERROR),
    ),
    "gromacs": (
        Pattern("Too many LINCS warnings", UNSTABLE),
        Pattern("Too many SETTLE warnings", UNSTABLE),
        Pattern("LINCS WARNING", UNSTABLE, warning=True),
        Pattern("1-4 interaction not within cut-off", UNSTABLE),
        Pattern("is not finite", UNSTABLE),
        Pattern("Fatal error:", ERROR),
        Pattern("Error in user input:", INPUT),
    ),
}

LOG_FILES: dict[str, tuple[str, ...]] = {
    "dftbplus": ("output.log",),
    "xtb": ("output.log",),
    "espresso": ("output.log", "bands/output.log"),
    "vasp": ("output.log", "OUTCAR"),
    "orca": ("output.log",),
    "cp2k": ("output.log",),
    "lammps": ("log.lammps", "output.log"),
    "gromacs": ("output.log", "grompp.log", "md.log", "adit.log"),
}
SCHEDULER_GLOBS = ("*.o[0-9]*", "*.e[0-9]*", "slurm-*.out")
_WARNING = re.compile(r"\bWARNING\b")
_BARE_MARKERS = ("WARNING!", "ERROR!")


def supported_codes() -> list[str]:
    return sorted(PATTERNS)


def cause_label(key: str) -> str:
    return {
        SCF: L("SCF (電子状態) が収束していない", "the SCF (electronic) cycle did not converge"),
        GEOMETRY: L("構造最適化が収束していない", "the geometry optimization did not converge"),
        WALLTIME: L("制限時間 (walltime) で止められた", "stopped at the walltime limit"),
        MEMORY: L("メモリ不足で止められた", "stopped for lack of memory"),
        MPI: L("MPI の異常終了", "MPI aborted"),
        DIAGONALIZATION: L("対角化の失敗", "diagonalization failed"),
        LOST_ATOMS: L("原子が失われた (LAMMPS の Lost atoms)", "atoms were lost (LAMMPS 'Lost atoms')"),
        UNSTABLE: L("力が発散した (系が不安定)", "forces diverged (the system is unstable)"),
        INPUT: L("入力の誤り", "input error"),
        RESOURCE: L("ジョブスケジューラの資源の上限で止められた", "stopped by a scheduler resource limit"),
        SIGNAL: L("シグナルで強制終了された (メモリ不足や上限超過のことが多い)", "killed by a signal (often memory or a limit)"),
        ERROR: L("コードがエラーを出して止まった", "the code stopped with an error message"),
    }.get(key, key)


@dataclass
class Hit:
    file: str
    line: int
    text: str
    cause: str | None = None
    warning: bool = False

    @property
    def where(self) -> str:
        return f"{self.file}:{self.line}"

    def as_dict(self) -> dict:
        return {"file": self.file, "line": self.line, "text": self.text, "cause": self.cause}


@dataclass
class Diagnostics:
    code: str
    supported: bool
    files: list[str] = field(default_factory=list)
    warnings: list[Hit] = field(default_factory=list)
    errors: list[Hit] = field(default_factory=list)
    n_warnings: int = 0

    def causes(self) -> list[str]:
        out: list[str] = []
        for h in self.errors:
            if h.cause and h.cause not in out:
                out.append(h.cause)
        return out

    def as_dict(self) -> dict:
        first = {}
        for h in self.errors:
            first.setdefault(h.cause, h)
        return {"code": self.code, "supported": self.supported, "files": list(self.files),
                "n_warnings": self.n_warnings, "warnings": [h.as_dict() for h in self.warnings],
                "errors": [h.as_dict() for h in self.errors],
                "causes": [{"key": k, "label": cause_label(k), "where": first[k].where, "text": first[k].text} for k in self.causes()]}

    def summary_lines(self) -> list[str]:
        return summary_lines_from_dict(self.as_dict())


def summary_lines_from_dict(d: dict) -> list[str]:
    lines: list[str] = []
    n = int(d.get("n_warnings") or 0)
    if n:
        shown = " / ".join(f"{h['file']}:{h['line']}  {h['text']}" for h in d.get("warnings", [])[:3])
        more = L(f"、ほか {n - 3} 件", f", {n - 3} more") if n > 3 else ""
        lines.append(L(f"警告: {n} 件 ({shown}{more})", f"warnings: {n} ({shown}{more})"))
    causes = d.get("causes") or []
    if causes:
        parts = [f"{cause_label(c['key'])} ({c['where']}  {c['text']})" for c in causes]
        lines.append(L("失敗の原因の候補 (出力の行をそのまま写しています): ", "possible causes of failure (quoted from the output): ") + " / ".join(parts))
    return lines


def _code_of(run_dir: Path) -> str:
    spec = run_dir / "spec.json"
    if spec.is_file():
        import json

        try:
            code = json.loads(spec.read_text(encoding="utf-8")).get("method", {}).get("code", "")
            if code:
                return code
        except ValueError:
            pass
    from adit.analysis.readers import detect_code

    return detect_code(run_dir) or ""


def _excerpt(text: str) -> str:
    text = text.strip()
    return text if len(text) <= EXCERPT_LENGTH else text[:EXCERPT_LENGTH - 1] + "…"


def _scan_file(path: Path, name: str, patterns: tuple[Pattern, ...], diag: Diagnostics) -> None:
    pending: Hit | None = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for number, raw in enumerate(f, 1):
            line = raw.rstrip("\n")
            text = line.strip()
            if not text:
                continue
            if pending is not None:
                # DFTB+ prints "WARNING!" / "ERROR!" alone and the message on the next line.
                pending.text = _excerpt(f"{pending.text} {text}")
                pending = None
            hit = next((p for p in patterns if p.matches(line)), None)
            if hit is not None and not hit.warning:
                if len(diag.errors) < MAX_STORED_HITS:
                    diag.errors.append(Hit(name, number, _excerpt(text), hit.cause))
                    pending = diag.errors[-1] if text in _BARE_MARKERS else None
                continue
            if hit is not None or _WARNING.search(line):
                diag.n_warnings += 1
                if len(diag.warnings) < MAX_STORED_HITS:
                    diag.warnings.append(Hit(name, number, _excerpt(text), hit.cause if hit else None, True))
                    pending = diag.warnings[-1] if text in _BARE_MARKERS else None


def scan_diagnostics(run_dir: Path | str, code: str = "") -> Diagnostics:
    """Scan the logs of a run directory for known warnings and failure messages."""
    d = Path(run_dir).expanduser()
    code = code or _code_of(d)
    diag = Diagnostics(code=code, supported=code in PATTERNS)
    patterns = PATTERNS.get(code, ()) + COMMON
    names = list(LOG_FILES.get(code, ("output.log",)))
    for pattern in SCHEDULER_GLOBS:
        names += sorted(p.name for p in d.glob(pattern) if p.is_file())
    seen: set[str] = set()
    for name in names:
        path = d / name
        if name in seen or not path.is_file():
            continue
        seen.add(name)
        diag.files.append(name)
        _scan_file(path, name, patterns, diag)
    return diag


def diagnostics_lines(run_dir: Path | str, code: str = "") -> list[str]:
    try:
        return scan_diagnostics(run_dir, code).summary_lines()
    except OSError:
        return []


_QE_BASE = "https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py"


def remedy_lines(code: str) -> list[str]:
    """README lines: what the code's own documents say to do when a run fails, with the source of each."""
    items = _REMEDIES.get(code)
    if items is None:
        return []
    head = [L("== 失敗したときに文献が挙げる対処 ==", "== What the documents suggest when a run fails =="),
            L("  ADIT が出力を走査して見つけた印 (解析の要約の「失敗の原因の候補」) ごとに、その計算コードの文書が挙げる対処を写したものです。",
              "  For each message ADIT finds in the output (the 'possible causes of failure' line of the analysis summary), this quotes what the code's own documents suggest."),
            L("  ADIT の推奨ではありません。数値や手順は出典のとおりで、系に合うかはご自身で判断してください。",
              "  These are not ADIT's recommendations; the numbers and steps are as in the sources, and whether they suit your system is your call.")]
    return head + [f"  - {L(ja, en)}" for ja, en in items()] + _common_remedies() + [""]


def _common_remedies() -> list[str]:
    return [L("  - 「DUE TO TIME LIMIT」(Slurm) / 「=>> PBS: job killed: walltime … exceeded limit …」(PBS): 制限時間で止められた印。出典: Slurm slurmstepd "
              "https://github.com/SchedMD/slurm/blob/master/src/slurmd/slurmstepd/req.c、OpenPBS pbs_mom https://github.com/openpbs/openpbs/blob/master/src/resmom/mom_main.c",
              "  - 'DUE TO TIME LIMIT' (Slurm) / '=>> PBS: job killed: walltime … exceeded limit …' (PBS): the walltime limit was hit. Sources: Slurm slurmstepd "
              "https://github.com/SchedMD/slurm/blob/master/src/slurmd/slurmstepd/req.c, OpenPBS pbs_mom https://github.com/openpbs/openpbs/blob/master/src/resmom/mom_main.c"),
            L("  - 「oom_kill event … OOM Killed」(Slurm) / 「=>> PBS: job killed: mem …」(PBS) / 「Out of memory: Killed process」(Linux) / 行頭近くの「Killed」: メモリ不足で強制終了された印。"
              "出典: Slurm https://github.com/SchedMD/slurm/blob/master/src/plugins/task/cgroup/task_cgroup_memory.c、Linux https://github.com/torvalds/linux/blob/master/mm/oom_kill.c",
              "  - 'oom_kill event … OOM Killed' (Slurm) / '=>> PBS: job killed: mem …' (PBS) / 'Out of memory: Killed process' (Linux) / a 'Killed' line: killed for lack of memory. "
              "Sources: Slurm https://github.com/SchedMD/slurm/blob/master/src/plugins/task/cgroup/task_cgroup_memory.c, Linux https://github.com/torvalds/linux/blob/master/mm/oom_kill.c"),
            L("  - 「MPI_ABORT was invoked on rank …」(Open MPI) / 「application called MPI_Abort」「BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES」(MPICH): MPI のプロセスが異常終了した印。"
              "原因はその前の計算コードの行にあります。出典: Open MPI help-mpi-api.txt (配布物に同梱)、MPICH https://github.com/pmodels/mpich/blob/main/src/pm/hydra/mpiexec/pmiserv_cb.c",
              "  - 'MPI_ABORT was invoked on rank …' (Open MPI) / 'application called MPI_Abort', 'BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES' (MPICH): an MPI process aborted; "
              "the reason is in the code's lines before it. Sources: Open MPI help-mpi-api.txt (shipped with Open MPI), MPICH https://github.com/pmodels/mpich/blob/main/src/pm/hydra/mpiexec/pmiserv_cb.c")]


def _dftbplus() -> list[tuple[str, str]]:
    recipes = "https://dftbplus-recipes.readthedocs.io/en/latest/basics/firstcalc.html"
    manual = "https://github.com/dftbplus/dftbplus/blob/main/doc/dftb+/manual/dftbp.tex"
    return [
        (f"「SCC is NOT converged, maximal SCC iterations exceeded」: DFTB+ recipes は、収束の判定が SccTolerance (既定 1e-5 電子)、反復の上限が MaxSccIterations (既定 100) で、"
         f"上限に達すると「そこまでの電荷で全エネルギーを計算して警告を出して止まる」と説明しています。出典: {recipes}",
         f"'SCC is NOT converged, maximal SCC iterations exceeded': the DFTB+ recipes explain that convergence is judged by SccTolerance (default 1e-5 electrons) with MaxSccIterations "
         f"(default 100) as the limit, and that on reaching it the code 'calculates the total energy using the charges obtained so far and stops with an appropriate warning message'. Source: {recipes}"),
        (f"同上: DFTB+ manual は、静的な構造でない計算 (Driver が空でない) では ConvergentSccOnly で「収束しなくても先へ進む」ようにできると書いています (MaxSCCIterations の項)。出典: {manual}",
         f"same: the DFTB+ manual says that for non-static runs (Driver not empty) ConvergentSccOnly can override the stop (MaxSCCIterations entry). Source: {manual}"),
        (f"「!!! Geometry did NOT converge!」: manual の MaxSteps は「収束する前にこの回数で最適化を止める」上限で、-1 で事実上無制限です。出典: {manual}",
         f"'!!! Geometry did NOT converge!': in the manual MaxSteps is the number of steps after which the optimization stops unless converged; -1 means practically unlimited. Source: {manual}"),
        ("「WARNING! -> Current stacksize not set to unlimited」: DFTB+ 自身の警告で、ulimit -s unlimited を勧めています (クラスタ用の submit.sh には入れてあります)。出典: DFTB+ の出力そのもの",
         "'WARNING! -> Current stacksize not set to unlimited': DFTB+'s own warning, advising ulimit -s unlimited (the cluster submit.sh already has it). Source: the DFTB+ output itself"),
    ]


def _vasp() -> list[tuple[str, str]]:
    wiki_opt = "https://www.vasp.at/wiki/index.php/Structure_optimization"
    wiki_el = "https://www.vasp.at/wiki/index.php/Troubleshooting_electronic_convergence"
    wiki_nelm = "https://www.vasp.at/wiki/index.php/NELM"
    f_zbrent = "https://www.vasp.at/forum/viewtopic.php?t=13434"
    f_edddav = "https://www.vasp.at/forum/viewtopic.php?t=10409"
    f_herm = "https://www.vasp.at/forum/viewtopic.php?t=1192"
    f_eddrmm = "https://www.vasp.at/forum/viewtopic.php?t=17822"
    return [
        (f"「ZBRENT: fatal error in bracketing」: VASP の出力自身が「please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue」と続けます。"
         f"VASP wiki は ZBRENT の出力を「補正ステップの誤差が大きすぎて線探索をさらに細かくしている印」と説明しています。出典: {wiki_opt}、公式フォーラム {f_zbrent}",
         f"'ZBRENT: fatal error in bracketing': VASP's own output continues 'please rerun with smaller EDIFF, or copy CONTCAR to POSCAR and continue'. "
         f"The VASP wiki describes ZBRENT output as the sign that 'the error in the corrector step was too large and the line search is further refined'. Sources: {wiki_opt}, official forum {f_zbrent}"),
        (f"「Error EDDDAV: Call to ZHEGV failed」: 公式フォーラムで VASP の開発者は、ALGO を変える、直前のイオンステップの座標から始め直す、原子が近すぎないか見る、"
         f"イオンループの中で電子状態が収束しているか確かめる、を挙げています。出典: {f_edddav}",
         f"'Error EDDDAV: Call to ZHEGV failed': on the official forum a VASP developer lists switching ALGO, restarting from an earlier ionic step, checking for atoms too close together, "
         f"and making sure the electronic steps converge inside the ionic loop. Source: {f_edddav}"),
        (f"「WARNING: Sub-Space-Matrix is not hermitian」: 公式フォーラムの管理者は、入力の構造か緩和の途中の構造が不合理なことが主因で、IBRION を変えるか POTIM を小さくする、"
         f"LAPACK の入れ方を確かめる、を挙げています。出典: {f_herm}",
         f"'WARNING: Sub-Space-Matrix is not hermitian': the forum administrator names an unreasonable input or intermediate geometry as the main cause and suggests changing IBRION or "
         f"reducing POTIM, and checking the LAPACK installation. Source: {f_herm}"),
        (f"「WARNING in EDDRMM: call to ZHEGV failed」: VASP の開発者は「警告であってエラーではなく、機械精度に近づいたときに出る。無視してよい」と答えています。出典: {f_eddrmm}",
         f"'WARNING in EDDRMM: call to ZHEGV failed': a VASP developer answers that it is a warning, not an error, occurring close to machine precision, and can be ignored. Source: {f_eddrmm}"),
        (f"電子状態が NELM 回で収束しない (OSZICAR の反復が NELM に達する): VASP wiki は、INCAR を最小にして k 点・ENCUT・PREC を下げて切り分ける、ISMEAR を確かめる、"
         f"NBANDS を増やす、ALGO を切り替える、IALGO=4X/5X なら TIME を変える、の順を挙げています。NELM の項は「IALGO/ALGO、LSUBROT、混合のパラメータを見直す」。出典: {wiki_el}、{wiki_nelm}",
         f"the electronic loop does not converge within NELM steps (OSZICAR reaches NELM): the VASP wiki lists, in order, a minimal INCAR with reduced k-points/ENCUT/PREC, checking ISMEAR, "
         f"increasing NBANDS, switching ALGO, and changing TIME for IALGO=4X/5X; the NELM page says to reconsider 'IALGO or ALGO, LSUBROT, and the mixing parameters'. Sources: {wiki_el}, {wiki_nelm}"),
    ]


def _espresso() -> list[tuple[str, str]]:
    b = _QE_BASE
    note_ja = "数値は aiida-quantumespresso の実装値"
    note_en = "the numbers are the values implemented in aiida-quantumespresso"
    return [
        (f"「convergence NOT achieved after N iterations」: aiida-quantumespresso の PwBaseWorkChain は mixing_beta を 0.8 倍にして、前の計算から再開します ({note_ja}: delta_factor_mixing_beta = 0.8)。"
         f"relax の途中で起きた場合は、出力の構造に更新して最初から計算し直します。出典: {b}",
         f"'convergence NOT achieved after N iterations': the PwBaseWorkChain of aiida-quantumespresso multiplies mixing_beta by 0.8 and restarts from the previous run ({note_en}: "
         f"delta_factor_mixing_beta = 0.8); when it happens during a relax it also takes the output structure and restarts from scratch. Source: {b}"),
        (f"「Maximum CPU time exceeded」: 出力の構造に更新し、CONTROL の restart_mode = 'restart' で続きを計算します。同ワークチェーンは max_seconds を制限時間の 0.95 倍に設定しています "
         f"({note_ja}: delta_factor_max_seconds = 0.95)。出典: {b}",
         f"'Maximum CPU time exceeded': update to the output structure and continue with restart_mode = 'restart' in CONTROL; the work chain sets max_seconds to 0.95 x the walltime "
         f"({note_en}: delta_factor_max_seconds = 0.95). Source: {b}"),
        (f"対角化の失敗 (「S matrix not positive definite」「problems computing cholesky」「zhegvd failed」「too many bands are not converged」「eigenvectors failed to converge」「Error in routine broyden」): "
         f"ELECTRONS の diagonalization を david / paro / cg の別のものに順に切り替え、全部試したら止めます。出典: {b}",
         f"diagonalization failures ('S matrix not positive definite', 'problems computing cholesky', 'zhegvd failed', 'too many bands are not converged', 'eigenvectors failed to converge', "
         f"'Error in routine broyden'): switch diagonalization in ELECTRONS to another of david / paro / cg in turn, and give up once all are tried. Source: {b}"),
        (f"「history already reset at previous step」(BFGS の履歴の失敗) / 「bfgs failed」: relax では ion_dynamics = 'damp' にする (対称性の低い系では、まず構造を標準偏差 0.01 の乱数で揺らす: "
         f"{note_ja} rattle_stdev = 0.01)。vc-relax では trust_radius_min を 0.1 倍にし ({note_ja} delta_factor_trust_radius_min = 0.1)、それでも駄目なら ion_dynamics = 'damp' と cell_dynamics = 'damp-w'。出典: {b}",
         f"'history already reset at previous step' (BFGS history failure) / 'bfgs failed': for relax set ion_dynamics = 'damp' (for low-symmetry systems first rattle the structure with a random "
         f"displacement of standard deviation 0.01: {note_en} rattle_stdev = 0.01); for vc-relax multiply trust_radius_min by 0.1 ({note_en} delta_factor_trust_radius_min = 0.1), "
         f"then fall back to ion_dynamics = 'damp' with cell_dynamics = 'damp-w'. Source: {b}"),
        (f"構造最適化が nstep 内に収束しない: 出力の構造に更新し、電荷密度を読んで再開します。出典: {b}",
         f"the ionic loop does not converge within nstep: update to the output structure and restart from the charge density. Source: {b}"),
        (f"「Not enough space allocated for radial FFT」(vc-relax で体積が大きく縮んだ): CELL の cell_factor を 2 倍にして (既定 2 → 4) 最初から計算し直します。出典: {b}",
         f"'Not enough space allocated for radial FFT' (large volume contraction in vc-relax): double cell_factor in CELL (default 2 → 4) and restart from scratch. Source: {b}"),
        (f"relax の途中で SCF の精度が上がらなくなった: IONS の upscale を 0.3 倍にして (最小 1) 再開します ({note_ja}: delta_factor_upscale = 0.3)。出典: {b}",
         f"the SCF accuracy gets stuck during a relax: multiply upscale in IONS by 0.3 (minimum 1) and restart ({note_en}: delta_factor_upscale = 0.3). Source: {b}"),
        (f"最高のバンドが占有されている (バンド数の不足): nbnd を「5 % か 4 本の大きいほう」だけ増やします ({note_ja}: delta_factor_nbnd = 0.05、delta_minimum_nbnd = 4)。出典: {b}",
         f"the highest band is occupied (too few bands): increase nbnd by the larger of 5 % and 4 bands ({note_en}: delta_factor_nbnd = 0.05, delta_minimum_nbnd = 4). Source: {b}"),
    ]


def _xtb() -> list[tuple[str, str]]:
    sp = "https://xtb-docs.readthedocs.io/en/latest/sp.html"
    opt = "https://xtb-docs.readthedocs.io/en/latest/optimization.html"
    return [
        (f"「Self consistent charge iterator did not converge」: xtb の文書は「電子温度を上げて収束させ、その結果から通常の温度で再開する」ことを挙げています: "
         f"xtb coord --etemp 1000.0 && xtb coord --restart (--etemp の既定は 300 K)。反復の上限は --iterations (既定 250)。出典: {sp}",
         f"'Self consistent charge iterator did not converge': the xtb docs suggest raising the electronic temperature and restarting from that result at the normal temperature: "
         f"xtb coord --etemp 1000.0 && xtb coord --restart (default --etemp is 300 K); the iteration limit is --iterations (default 250). Source: {sp}"),
        (f"「*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN N ITERATIONS ***」(空のファイル NOT_CONVERGED も書かれる): 反復の上限は --cycles、収束の厳しさは --opt のレベル (crude 〜 extreme、既定 normal)。"
         f"文書は「まず GFN0-xTB で最適化してから GFN2-xTB に切り替える」「気相で問題が出るときは ALPB 溶媒モデルを使う」ことを挙げています。出典: {opt}",
         f"'*** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN N ITERATIONS ***' (an empty file NOT_CONVERGED is also written): the limit is --cycles and the strictness is the --opt level "
         f"(crude to extreme, default normal); the docs suggest optimizing with GFN0-xTB first and then GFN2-xTB, and using the ALPB solvation model when the gas phase gives trouble. Source: {opt}"),
    ]


def _orca() -> list[tuple[str, str]]:
    ts = "https://www.faccts.de/docs/orca/6.1/manual/contents/quickstartguide/troubleshooting.html"
    opt = "https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/optimizations.html"
    return [
        (f"SCF が収束しない: ORCA manual 1.7.13 は !SlowConv を挙げ (「初期推定に近い局所解に収束することがある」と注意)、遷移金属では BP86 + SlowConv で先に計算して "
         f"その軌道を読み込み (ReadMOs、GuessMode CMatrix) B3LYP で計算し直す 2 段の Compound の例を載せています。出典: {ts}",
         f"the SCF does not converge: ORCA manual 1.7.13 suggests !SlowConv (noting it 'may converge to a local minimum solution that is closer to that of the initial guess') and shows a "
         f"two-step Compound example for transition metals: BP86 + SlowConv first, then B3LYP reading those orbitals (ReadMOs, GuessMode CMatrix). Source: {ts}"),
        (f"構造最適化が収束しない: manual 1.7.14.3 は、最後の構造から最適化をやり直す、安価な手法や小さい基底で先に最適化する、緩い収束条件から始めて締める、"
         f"軌跡を見て問題の部分に拘束を掛ける、trust radius や最大ステップを変える、!COPT (直交座標) を挙げています。出典: {ts}",
         f"the geometry optimization does not converge: manual 1.7.14.3 suggests rerunning from the final geometry, optimizing first with a cheaper method or smaller basis, starting with "
         f"looser thresholds before tightening, constraining the troublesome part found in the trajectory, changing the trust radius or maximum step, and !COPT (Cartesian). Source: {ts}"),
        (f"「ORCA TERMINATED NORMALLY」があっても最適化が収束したとは限りません。manual 4.1.12 は「THE OPTIMIZATION HAS CONVERGED」の行で判断するよう書いています。出典: {opt}",
         f"'ORCA TERMINATED NORMALLY' does not mean the optimization converged; manual 4.1.12 says to look for 'THE OPTIMIZATION HAS CONVERGED' instead. Source: {opt}"),
        (f"「Please increase MaxCore」「OUT OF MEMORY ERROR!」: manual 1.7.5 は MaxCore が「コア 1 つあたり MB」で、MaxCore × nprocs ≤ 0.75 × (使えるメモリ) にするよう書いています。出典: {ts}",
         f"'Please increase MaxCore' / 'OUT OF MEMORY ERROR!': manual 1.7.5 says MaxCore is in MB per processing core and should satisfy MaxCore x nprocs <= 0.75 x (available memory). Source: {ts}"),
        (f"突然止まった: manual 1.7.3 は、出力の最後の行と、先頭近くの WARNINGS の節を読む、時間・ディスク・メモリの上限ならジョブスケジューラの出力を見る、と書いています。出典: {ts}",
         f"sudden termination: manual 1.7.3 says to read the last lines of the output and the WARNINGS section near the start, and the scheduler output for time, disk or memory limits. Source: {ts}"),
    ]


def _cp2k() -> list[tuple[str, str]]:
    scf = "https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/SCF.html"
    chol = "https://www.cp2k.org/faq:cholesky_decomp_failed"
    return [
        (f"「SCF run NOT converged」: CP2K manual の SCF 節では、内側の反復の上限が MAX_SCF (既定 50)、IGNORE_CONVERGENCE_FAILURE = T にすると止めずに警告だけにする、"
         f"smearing を使うなら ADDED_MOS が要る、初期推定は SCF_GUESS で選ぶ、とあります。CP2K 自身の停止メッセージも IGNORE_CONVERGENCE_FAILURE を案内します。出典: {scf}",
         f"'SCF run NOT converged': the SCF section of the CP2K manual gives MAX_SCF (default 50) as the inner-loop limit, IGNORE_CONVERGENCE_FAILURE = T to continue with a warning only, "
         f"ADDED_MOS as required for smearing, and SCF_GUESS for the initial guess; CP2K's own abort message also points to IGNORE_CONVERGENCE_FAILURE. Source: {scf}"),
        (f"「CPASSERT failed」が cp_fm_cholesky で出る: CP2K FAQ は、OT の前処理で重なり行列がほぼ特異になったのが原因で、別の PRECONDITIONER を試す、EPS_DEFAULT か EPS_PGF_ORB を小さくする、"
         f"基底が過剰なら小さい基底にする、を挙げています。出典: {chol}",
         f"'CPASSERT failed' in cp_fm_cholesky: the CP2K FAQ attributes it to a (nearly) singular overlap matrix in the OT preconditioner and suggests another PRECONDITIONER, "
         f"a smaller EPS_DEFAULT or EPS_PGF_ORB, and a smaller basis when the basis is over-complete. Source: {chol}"),
        ("「*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***」: 最適化の上限 (MOTION / GEO_OPT の MAX_ITER) に達した印です。出典: CP2K の出力そのもの",
         "'*** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***': the optimization hit its step limit (MAX_ITER in MOTION / GEO_OPT). Source: the CP2K output itself"),
    ]


def _lammps() -> list[tuple[str, str]]:
    det = "https://docs.lammps.org/Errors_details.html"
    return [
        (f"「Lost atoms: original N current M」: LAMMPS の文書 5.2.9 は、多くの場合は速度が極端に大きくなった不安定化が原因で、「Fast moving atoms」(近すぎる原子を消すか先に最小化する、"
         f"初期の緩和に fix nve/limit や fix dt/reset) と「Neighbor list settings」を見るよう書いています。スパッタリングのように失われて当然の系では thermo_modify lost warn/ignore。出典: {det}",
         f"'Lost atoms: original N current M': LAMMPS docs 5.2.9 say it is mostly an instability with extreme velocities and point to 'Fast moving atoms' (delete or minimize close contacts, "
         f"fix nve/limit or fix dt/reset during initial relaxation) and 'Neighbor list settings'; when losing atoms is expected (sputtering) use thermo_modify lost warn/ignore. Source: {det}"),
        (f"「Out of range atoms - cannot compute PPPM」: 5.2.5 は、原子が速すぎる、neighbor の skin が小さすぎる、neighbor list の更新が少なすぎる、が原因で、skin を増やす、"
         f"timestep を小さくする、更新を増やす、初期配置の近接を避ける、を挙げています。出典: {det}",
         f"'Out of range atoms - cannot compute PPPM': 5.2.5 names fast atoms, a too-small neighbor skin and too infrequent neighbor list updates, and suggests a larger skin, "
         f"a smaller timestep, more frequent updates and avoiding close contacts in the setup. Source: {det}"),
        (f"「Non-numeric positions - simulation unstable」「Bond atoms missing」: 5.2.6〜5.2.7 は、力が数として表せなくなった (NaN / Inf) 不安定化で、重なった原子を delete_atoms で消す、"
         f"先に minimize する、ソフトなポテンシャルで初期平衡化する、を挙げています。出典: {det}",
         f"'Non-numeric positions - simulation unstable' / 'Bond atoms missing': 5.2.6-5.2.7 describe forces that can no longer be represented (NaN / Inf) and suggest removing overlapping atoms "
         f"with delete_atoms, minimizing first, and equilibrating with soft potentials. Source: {det}"),
        ("「ERROR: … (file.cpp:line)」: LAMMPS のエラー行の形で、file:line はその判定を出した LAMMPS のソースの場所です。出典: https://docs.lammps.org/Errors_messages.html",
         "'ERROR: … (file.cpp:line)': the form of a LAMMPS error line; file:line is where in the LAMMPS source the check lives. Source: https://docs.lammps.org/Errors_messages.html"),
    ]


def _gromacs() -> list[tuple[str, str]]:
    err = "https://manual.gromacs.org/current/user-guide/run-time-errors.html"
    blow = "https://manual.gromacs.org/current/user-guide/terminology.html#blowing-up"
    feat = "https://manual.gromacs.org/current/user-guide/mdrun-features.html"
    return [
        (f"「LINCS WARNING」「Too many LINCS warnings」(SETTLE / SHAKE も同じ): GROMACS の文書は、系が発散 (blowing up) するときに拘束が最初に壊れるだけで、拘束そのものより系の側に物理的に無理がある印だと書いています。"
         f"blowing up の項が挙げる原因は、最小化の不足、立体的に重なった初期構造、拘束に対して大きすぎる timestep、不適切な圧力・温度制御、実座標から遠い位置拘束など。対処は「力が大きくならないようにするか、timestep を小さくする」。出典: {err}、{blow}",
         f"'LINCS WARNING' / 'Too many LINCS warnings' (also SETTLE / SHAKE): the GROMACS docs say the constraints are just the first thing to fail when the system blows up, a sign of something "
         f"physically unrealistic in the system; the 'blowing up' entry lists insufficient minimization, steric clashes in the start structure, a timestep too large for the constraints, "
         f"inappropriate pressure/temperature coupling and position restraints far from the coordinates, and the remedy 'make sure the forces do not get that large, or use a smaller timestep'. Sources: {err}, {blow}"),
        (f"「1-4 interaction not within cut-off」: もう一度エネルギー最小化する、温度が高すぎないか、timestep が大きすぎないか。出典: {err}",
         f"'1-4 interaction not within cut-off': another round of energy minimization; otherwise the temperature may be too high or the timestep too large. Source: {err}"),
        (f"「Stepsize too small, or no change in energy. Converged to machine precision, but not to the requested Fmax」: エラーとは限らず、別の最小化法か倍精度版で先へ進めるとあります。出典: {err}",
         f"'Stepsize too small, or no change in energy. Converged to machine precision, but not to the requested Fmax': not necessarily an error; a different minimizer or double precision may go further. Source: {err}"),
        (f"「Energy minimization has stopped because the force on at least one atom is not finite」: 2 つの原子が近すぎる印で、ソフトコアのポテンシャルで最小化できることがあるとあります。出典: {err}",
         f"'Energy minimization has stopped because the force on at least one atom is not finite': two atoms are too close; soft-core potentials can sometimes minimize such systems. Source: {err}"),
        (f"制限時間: gmx mdrun -maxh で時間切れの少し前に止めてチェックポイントを書き、続きから再開できます。出典: {feat}",
         f"walltime: gmx mdrun -maxh stops shortly before the limit and writes a checkpoint to continue from. Source: {feat}"),
    ]


_REMEDIES = {"dftbplus": _dftbplus, "vasp": _vasp, "espresso": _espresso, "xtb": _xtb, "orca": _orca,
             "cp2k": _cp2k, "lammps": _lammps, "gromacs": _gromacs}
