"""Write a run directory (inputs, submit.sh, spec.json, README.txt, parameters) and read it back."""

from __future__ import annotations

from adit.errors import AditError
import shlex
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from adit.lang import L, pick
from adit.codes import GENERATORS, GenerationError
from adit.codes.base import ReadmeNotes
from adit.config import Config, config_path
from adit.spec import CalculationSpec
from adit.validate import ValidationError, validate

SPEC_FILE = "spec.json"
SUBMIT_FILE = "submit.sh"
README_FILE = "README.txt"
TRANSFER_FILE = "transfer_and_submit.sh"
ANALYZE_FILE = "analyze.py"


class ProjectError(AditError):

    def __init__(self, message: str, errors: list[ValidationError] | None = None):
        super().__init__(message)
        self.errors = errors or []


class OutputNotEmpty(ProjectError):
    pass


@dataclass
class ProjectFiles:

    texts: dict[str, str] = field(default_factory=dict)
    copies: dict[str, Path] = field(default_factory=dict)


def version_trap(probe: tuple[str, str]) -> str:
    from adit.provenance import VERSION_FILE

    name, pattern = probe
    return (f"adit_dir=$(pwd)\n"
            f"trap 'grep -m1 -E \"{pattern}\" \"$adit_dir/{name}\" > \"$adit_dir/{VERSION_FILE}\" 2>/dev/null' EXIT")


def _transfer_commands(profile, output_dir) -> str:
    """The copy-and-paste commands for sending this directory to the cluster and submitting it there.

    ADIT never runs them: the person does, in a terminal."""
    here = Path(output_dir).expanduser().resolve() if output_dir else None
    name = shlex.quote(here.name) if here else L("<このディレクトリの名前>", "<directory name>")
    host = profile.target or L("<ユーザー名>@<クラスタのホスト名>", "<user>@<cluster host name>")
    work = profile.remote_dir.strip() or L("<クラスタでの作業ディレクトリ>", "<work directory on the cluster>")
    todo = "" if profile.target and profile.remote_dir.strip() else L(
        "# <...> の部分は自分の値に置き換えてください (環境設定のプロファイルに host と remote_dir を書くと、ここが埋まります)。",
        "# Replace the <...> parts with your own values (set host and remote_dir in the profile to fill them in).")
    return "\n".join(x for x in [
        L("# ADIT が書き出したコマンドです。ADIT は実行しません。ターミナルに貼って使います。",
          "# Commands written by ADIT. ADIT does not run them; paste them into a terminal."),
        todo,
        "",
        L("# このファイルのある場所を自分で調べます (フォルダを移動しても、そのまま使えます)",
          "# Work out where this file is, so the commands keep working if the folder is moved"),
        'HERE="$(cd "$(dirname "$0")" && pwd)"',
        "",
        L("# 1. このディレクトリをクラスタへ送る", "# 1. Send this directory to the cluster"),
        f'rsync -av "$HERE" {host}:{work}/',
        "",
        L("# 2. クラスタにログインして投入する", "# 2. Log in to the cluster and submit"),
        f"ssh {host}",
        f"cd {work}/{name} && {profile.submit} {SUBMIT_FILE}",
        f"{profile.status}",
        "",
        L("# 3. 終わったら結果を手元へ戻す", "# 3. Bring the results back when it has finished"),
        f'rsync -av {host}:{work}/{name} "$(dirname "$HERE")"/',
        "",
    ] if x is not None) + ""


def build_project(spec: CalculationSpec, cfg: Config, *, output_dir: Path | str | None = None,
                  pre_command: str = "", extra_readme: list[str] | None = None,
                  extra_texts: dict[str, str] | None = None, drop: tuple[str, ...] = (), run_transform=None) -> ProjectFiles:
    """Assemble the files without writing them. Raises ProjectError when the settings are not valid."""
    errs = validate(spec, cfg, output_dir=output_dir)
    if errs:
        raise ProjectError(L("次の点を直すと生成できます:", "fix the following to generate:") + "\n" + "\n".join(f"  - {e}" for e in errs), errs)
    profile = cfg.profile(spec.runtime.profile)
    gen = GENERATORS[spec.method.code]
    from adit.scripts.render import render_submit
    from adit import provenance

    files = ProjectFiles()
    try:
        res = gen.resolve(spec, cfg)
        files.texts.update(gen.generate(spec, res))
        for name in drop:
            files.texts.pop(name, None)
        if extra_texts:
            files.texts.update(extra_texts)
        files.copies.update(gen.files_to_copy(spec, res))
        h = spec.handoff
        if h is not None and not h.at_run:
            for dest, src in h.files.items():
                files.copies[dest] = Path(h.previous_dir).expanduser() / src
        notes = gen.readme_notes(spec, res, files.copies)
        run_command = gen.run_command(spec, profile)
        if run_transform is not None:
            run_command = run_transform(run_command)
        sources = gen.parameter_sources(spec, res)
        probe = gen.version_probe(spec)
    except GenerationError as ex:
        raise ProjectError(str(ex)) from ex
    if pre_command:
        run_command = f"{pre_command} && {run_command}"
    if probe:
        run_command = version_trap(probe) + "\n" + run_command
    files.texts[SUBMIT_FILE] = render_submit(spec, profile, run_command)
    replaced: set[str] = set()
    if spec.handoff is not None and spec.handoff.at_run:
        from adit.handoff import rewritten_at_run

        replaced = set(spec.handoff.files or {}) | set(rewritten_at_run(spec.method.code, spec.handoff.velocities))
    prov = provenance.collect(files.copies, files.texts, sources, probe, replaced_at_run=replaced)
    files.texts[SPEC_FILE] = provenance.spec_json_with(spec, prov)
    from adit.analysis.report import analyze_script_text
    if gen.supports_analysis:
        files.texts[ANALYZE_FILE] = analyze_script_text()
    files.texts[README_FILE] = _readme(spec, profile, notes, output_dir, prov=prov, extra=extra_readme,
                                     settings_path=cfg.source_path)
    if profile.kind != "direct":
        files.texts[TRANSFER_FILE] = _transfer_commands(profile, output_dir)
    _check_names_stay_inside(list(files.texts) + list(files.copies))
    return files


def _check_names_stay_inside(names: list[str]) -> None:
    # Destination names come from generators and from handoff.files (user data): never leave output_dir.
    from pathlib import PureWindowsPath

    bad = [n for n in names if not n.strip() or Path(n).is_absolute() or PureWindowsPath(n).is_absolute()
           or PureWindowsPath(n).drive or ".." in Path(n).parts or ".." in PureWindowsPath(n).parts]
    if bad:
        raise ProjectError(L(f"出力ディレクトリの外を指すファイル名は書けません: {bad}",
                             f"file names that point outside the output directory cannot be written: {bad}"))


def write_project(spec: CalculationSpec, cfg: Config, output_dir: Path | str, *, overwrite: bool = False,
                  pre_command: str = "", extra_readme: list[str] | None = None,
                  extra_texts: dict[str, str] | None = None, drop: tuple[str, ...] = (), run_transform=None) -> list[Path]:
    out = Path(output_dir).expanduser()
    files = build_project(spec, cfg, output_dir=out, pre_command=pre_command, extra_readme=extra_readme,
                          extra_texts=extra_texts, drop=drop, run_transform=run_transform)
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise OutputNotEmpty(L(f"出力先 {out} には、すでにファイルがあります。黙って上書きしないよう、何も書かずに止めました。",
                               f"The output directory {out} already contains files; nothing was written, to avoid overwriting them silently."))
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for rel, text in files.texts.items():
        p = out / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8", newline="\n")
        written.append(p)
    for name in (SUBMIT_FILE, "make_potcar.sh"):
        if (out / name).exists():
            try:
                (out / name).chmod(0o755)
            except OSError:
                pass
    for rel, src in files.copies.items():
        p = out / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, p)
        written.append(p)
    return written


def load_project(output_dir: Path | str) -> CalculationSpec:
    p = Path(output_dir).expanduser() / SPEC_FILE
    if not p.is_file():
        raise ProjectError(L(f"{p} がありません", f"{p} not found"))
    return CalculationSpec.load(p)


def _local_time(iso: str) -> str:
    from datetime import datetime

    try:
        dt = datetime.fromisoformat(iso)
    except (TypeError, ValueError):
        return iso
    if dt.tzinfo is None:
        return dt.strftime("%Y-%m-%d %H:%M")
    dt = dt.astimezone()
    off = dt.strftime("%z")
    return dt.strftime("%Y-%m-%d %H:%M") + (f" (UTC{off[:3]}:{off[3:]})" if off else "")


def _path_hint(code: str, program: str) -> list[str]:
    if code == "vasp":
        return [L("     何も出なければ、VASP の実行ファイル (vasp_std など。VASP のライセンスを持つ研究室で入手・ビルドしたもの) の",
                  "     If nothing appears, add the directory holding the VASP executables (vasp_std etc., obtained or built by a group"),
                L("     あるディレクトリを PATH に足してから、もう一度確かめます。例: export PATH=<VASP の実行ファイルのあるディレクトリ>:$PATH",
                  "     with a VASP license) to PATH, then check again. e.g. export PATH=<directory with the VASP executables>:$PATH")]
    if code == "dcdftbmd":
        return [L("     DCDFTBMD 本体は開発元から入手してください。配布物には含まれません。実行ファイルの場所はご自身の環境設定に指定します。",
                  "     Obtain DCDFTBMD from its developers; it is not bundled. Set the executable path in your own settings.")]
    if code == "orca":
        return [L("     何も出なければ、入手した ORCA を展開したディレクトリを PATH に足してから、もう一度確かめます。",
                  "     If nothing appears, add the directory where you unpacked ORCA to PATH, then check again."),
                L("     例: export PATH=<ORCA を展開したディレクトリ>:$PATH",
                  "     e.g. export PATH=<directory where ORCA was unpacked>:$PATH")]
    if code == "mlip":
        return [L("     python3 が出ても、ase とモデルのパッケージが入っていなければ走りません。下の「実行する前に用意すること」の pip の行で入れます。",
                  "     Even if python3 is found, it runs only when ase and the model package are installed; see the pip line under 'Before running'.")]
    return [L("     何も出なければ、計算ソフトを入れた conda 環境 (ソフトごとに分けたインストール先) を有効にしてから、",
              "     If nothing appears, activate the conda environment (a separate install location per program) that has it,"),
            L(f"     もう一度確かめます。例: conda activate <環境名>  (入れていなければ: conda install -c conda-forge {_CONDA_PKG.get(code, code)})",
              f"     then check again. e.g. conda activate <env name>  (if it is not installed: conda install -c conda-forge {_CONDA_PKG.get(code, code)})")]


_CONDA_PKG = {"dftbplus": "dftbplus", "xtb": "xtb", "espresso": "qe", "cp2k": "cp2k", "lammps": "lammps", "gromacs": "gromacs",
              "openmm": "openmm", "psi4": "psi4", "abinit": "abinit"}


def _cluster_profile_example(code: str) -> list[str]:
    return [
        "       [profiles.remote]",
        L('       kind = "pbs"                                # Slurm のクラスタなら "slurm"',
          '       kind = "pbs"                                # "slurm" for a Slurm cluster'),
        L('       header_extra = ["#PBS -q <キュー名>"]        # Slurm なら ["#SBATCH --partition=<パーティション名>"]',
          '       header_extra = ["#PBS -q <queue>"]           # Slurm: ["#SBATCH --partition=<partition>"]'),
        "       [profiles.remote.code_modules]",
        L(f'       {code} = ["<module 名>"]                   # 計算ソフトを使えるようにする module の名前 (クラスタで module avail と打つと一覧が出ます)',
          f'       {code} = ["<module name>"]                 # module that makes the program available (module avail on the cluster lists them)'),
    ]


def _continuation_lines(spec: CalculationSpec) -> list[str]:
    c = spec.meta.continued_from
    if not c:
        return []
    lines = [L("== 前の計算から引き継いだもの ==", "== Carried over from the previous run =="),
             L(f"  前の計算: {c.get('dir')} ({c.get('code')}、{c.get('task')})", f"  Previous run: {c.get('dir')} ({c.get('code')}, {c.get('task')})")]
    lines += [f"  + {x}" for x in c.get("carried", [])]
    lines += [f"  - {x}" for x in c.get("not_carried", [])]
    return lines + [""]


CODE_NAMES = {"dftbplus": "DFTB+", "dcdftbmd": "DCDFTBMD", "nwchem": "NWChem", "gaussian": "Gaussian",
              "gamess": "US GAMESS", "qchem": "Q-Chem", "grrm": "GRRM17", "openmx": "OpenMX", "amber": "Amber",
              "namd": "NAMD", "vasp": "VASP", "xtb": "xtb", "espresso": "Quantum ESPRESSO (pw.x)", "orca": "ORCA",
              "cp2k": "CP2K", "lammps": "LAMMPS", "gromacs": "GROMACS", "openmm": "OpenMM", "abinit": "ABINIT",
              "psi4": "Psi4"}


def code_display_name(code: str) -> str:
    if code == "mlip":
        return L("機械学習ポテンシャル (ASE)", "machine-learning potential (ASE)")
    return CODE_NAMES.get(code, code)


def _readme(spec: CalculationSpec, profile, notes: ReadmeNotes, output_dir: Path | str | None = None, *,
            prov: dict | None = None, extra: list[str] | None = None, settings_path: Path | None = None) -> str:
    t, r, m = spec.task, spec.runtime, spec.method
    gen = GENERATORS[m.code]
    task = pick(
        {"geometry_optimization": "構造最適化", "molecular_dynamics": "分子動力学 (MD)", "vibrations": "振動解析", "band_structure": "バンド計算"},
        {"geometry_optimization": "geometry optimization", "molecular_dynamics": "molecular dynamics (MD)", "vibrations": "vibrations", "band_structure": "band structure"},
    ).get(t.type, L("一点計算", "single point"))
    code_name = code_display_name(m.code)
    natoms = len(spec.structure.atoms.symbols)
    if output_dir:
        dir_name = shlex.quote(Path(output_dir).name)
    else:
        dir_name = L("<このディレクトリの名前>", "<directory name>")
    local_dir = L("'<この計算ディレクトリのパス>'", "'<path to this calculation directory>'")
    host = profile.target or L("<ユーザー名>@<クラスタのホスト名>", "<user>@<cluster host name>")
    work = profile.remote_dir.strip() or L("<クラスタでの作業ディレクトリ>", "<work directory on the cluster>")
    back = L("<手元の置き場所>", "<local folder>")
    kind = profile.kind
    where = {"direct": L("この PC で直接実行する", "run directly on this PC"), "pbs": L("PBS のクラスタ", "a PBS cluster"),
             "slurm": L("Slurm のクラスタ", "a Slurm cluster")}[kind]

    lines = [
        L(f"ADIT {spec.meta.app_version} が生成した {code_name} の計算ディレクトリです ({_local_time(spec.meta.created)})",
          f"{code_name} calculation directory generated by ADIT {spec.meta.app_version} ({_local_time(spec.meta.created)})"),
        "",
        L(f"計算の種類: {task} / 元素: {' '.join(spec.elements)} / 原子数: {natoms}",
          f"Task: {task} / elements: {' '.join(spec.elements)} / atoms: {natoms}"),
        L(f"計算コード: {code_name}", f"Code: {code_name}"),
        L(f"実行先: プロファイル {r.profile} ({where})", f"Target: profile {r.profile} ({where})"),
        (L("  (プロファイル = ADIT の環境設定に書いた実行先。CLI では spec.json の runtime.profile で選びます)",
           "  (profile = a named execution target in the ADIT settings; the CLI selects it with runtime.profile in spec.json)")
         if gen.cli_only else L("  (プロファイル = ADIT の環境設定に書いた「どこで・どう実行するか」の組。ADIT の画面の「プロファイル」で選びます)",
                             "  (profile = a named 'where and how to run' entry in the ADIT settings; chosen in the Profile field of ADIT)")),
        "",
        L("== このディレクトリのファイル ==", "== Files in this directory =="),
        L("  submit.sh     ジョブスクリプト (計算を実行する手順を書いたシェルスクリプト。下の「実行する」で使います)",
          "  submit.sh     job script (a shell script with the steps that run the calculation; used below)"),
        (L("  spec.json     設定一式。CLI の adit-gen で再生成するときに使います (GUI では編集できません)",
           "  spec.json     all settings; use it with adit-gen to regenerate this project (not editable in the GUI)")
         if gen.cli_only else L("  spec.json     ADIT で決めた設定一式。ADIT の「ファイル」→「計算設定 (spec.json) を開く…」で読み込むと、同じ設定を画面に戻せます",
                             "  spec.json     all the settings chosen in ADIT; File → Open calculation settings (spec.json)… in ADIT brings them back")),
        *([L("  analyze.py    解析スクリプト (下の「結果の見方」を参照)",
             "  analyze.py    analysis script (see 'Reading the results' below)")] if gen.supports_analysis else []),
        *notes.files,
        "",
    ]
    if t.type == "vibrations":
        lines += [L("== 振動解析に使う構造 ==", "== Structure used for the vibrational analysis =="),
                  L("  この計算は、入力した構造をそのまま基準にして振動を計算します。ADIT は、この構造が最適化済みかどうかを判断しません。",
                    "  This calculation uses the input structure as-is as the reference for the vibrations. ADIT does not decide whether it has been optimized.")]
        if not gen.cli_only:
            lines += [L("  最適化の出力から始める場合は「前の計算の続き」を使うか、CLI で adit-gen --continue-from <最適化のディレクトリ> を使います。",
                        "  To start from an optimization result, use Continue a previous calculation, or adit-gen --continue-from <optimization directory> on the command line.")]
        lines.append("")
    lines += _continuation_lines(spec)
    if extra:
        lines += [*extra, ""]
    if notes.prepare:
        lines += [L("== 実行する前に用意すること ==", "== Before running =="), *notes.prepare, ""]

    lines.append(L("== この PC で実行する ==", "== Run on this PC =="))
    if kind == "direct":
        lines += [
            L("  以下はターミナル (コマンドを 1 行ずつ打ち込んで PC を操作する画面。Windows なら WSL の Ubuntu) で行います。",
              "  Do the following in a terminal (the window where you type commands line by line; on Windows, Ubuntu on WSL)."),
            L("  1. このディレクトリへ移動します (cd = 作業する場所を変えるコマンド)",
              "  1. Move into this directory (cd = the command that changes where you work)"),
            L("     <...> を、いまこのディレクトリを置いている場所のパスに置き換えてください。別の PC へ写した場合は、写し先のパスを使います。",
              "     Replace <...> with this directory's current path. If you copied it to another PC, use the path on that PC."),
            f"       cd {local_dir}",
            L(f"  2. {notes.program} が PATH (コマンドを探す場所の一覧) にあるか確かめます",
              f"  2. Check that {notes.program} is on PATH (the list of places where commands are looked up)"),
            f"       command -v {notes.program}",
            L(f"     場所 (例: /home/.../bin/{notes.program}) が 1 行出れば準備できています。",
              f"     If one line with a location (e.g. /home/.../bin/{notes.program}) appears, you are ready."),
            *_path_hint(m.code, notes.program),
            L("  3. 実行します", "  3. Run it"),
            "       bash submit.sh",
            L("     画面に何も出なくても動いています。記録は output.log に書かれ、終わると次のコマンドを打てる状態に戻ります。",
              "     Nothing appears on screen while it runs; the log goes to output.log, and the prompt comes back when it finishes."),
            L("     途中経過は別のターミナルで  tail -f output.log  (Ctrl+C で表示だけ止まり、計算は続きます)",
              "     To watch progress, run  tail -f output.log  in another terminal (Ctrl+C stops the display only, not the calculation)"),
            "",
        ]
    else:
        lines += [
            L(f"  この submit.sh はクラスタ ({kind.upper() if kind == 'pbs' else 'Slurm'}) 用です。この PC で実行するときは、ADIT でプロファイルを",
              f"  This submit.sh is for a cluster ({kind.upper() if kind == 'pbs' else 'Slurm'}). To run on this PC, switch the profile in ADIT to"),
            L("  local (kind = \"direct\" のもの) に切り替えて生成し直してください。この README.txt もその手順に変わります。",
              "  local (one with kind = \"direct\") and generate again; this README.txt then describes that instead."),
            "",
        ]

    lines += [
        L("== 研究室のクラスタで実行したいとき ==", "== Running on your group's cluster =="),
        L("  クラスタ = 研究室や計算センターが共同で使う計算機の集まり。計算はジョブスケジューラ (PBS や Slurm。計算の順番待ちを",
          "  cluster = a set of shared computers run by your group or a computing center. You hand your calculation to the job scheduler"),
        L("  管理するソフト) に預けて実行します。預けた 1 件の計算を「ジョブ」と呼びます。",
          "  (PBS or Slurm; software that manages the queue), and each calculation handed over is called a 'job'."),
    ]
    if kind == "direct":
        lines += [
            L("  いまの submit.sh はこの PC 用です。クラスタで使うには、先に次の 2 つを行います。",
              "  The current submit.sh is for this PC. To use a cluster, first do these two things."),
            L("  a. ADIT の環境設定ファイルに、kind = \"pbs\" か \"slurm\" のプロファイルを足します。環境設定ファイルの場所:",
              "  a. Add a profile with kind = \"pbs\" or \"slurm\" to the ADIT settings file, which is at"),
            f"       {settings_path or config_path()}",
            L("     いちばん短い書き方は次のとおりです。ファイルの最後に書き足し、<...> を自分の値に置き換えます (# から行末までは説明で、",
              "     The shortest form is below. Append it to the end of the file and replace <...> with your values (from # to the end of a line"),
            L("     消してもかまいません)。キュー名や module 名 (クラスタで計算ソフトを使えるようにするための名前) はクラスタごとに違うので、",
              "     is a comment and may be removed). Queue names and module names (the names that make a program available on the cluster) differ"),
            L("     管理者か研究室の先輩に確かめてください。",
              "     between clusters, so ask the administrator or a senior member of your group."),
            *_cluster_profile_example(m.code),
            L("     ほかの項目 (投入コマンド、環境変数、実行コマンドなど) は ADIT の README.md の「クラスタで実行する場合」にあります。",
              "     Other items (submit command, environment variables, run commands, ...) are in ADIT's README.md, 'Running on a cluster'."),
            *([L("  b. spec.json の runtime.profile を追加したプロファイル名に変えて、adit-gen で生成し直します。",
                 "  b. Set runtime.profile in spec.json to the added profile name, then regenerate with adit-gen."),
               L("     新しい submit.sh は PBS / Slurm 用になり、README.txt にも投入と確認の手順が入ります。",
                 "     The new submit.sh targets PBS / Slurm, and README.txt includes the submit and status steps.")]
              if gen.cli_only else [
                  L("  b. ADIT でプロファイルをそれに切り替えて生成し直します。submit.sh が qsub (PBS) / sbatch (Slurm) で預ける",
                    "  b. Switch the profile in ADIT to it and generate again. submit.sh becomes a job script for qsub (PBS) / sbatch (Slurm),"),
                  L("     ジョブスクリプトになり、この README.txt にも、そのクラスタでの投入と確認のコマンドが入ります。",
                    "     and this README.txt then shows the submit and status commands for that cluster.")]),
            L("  そのあとの流れは次のとおりです (<...> は自分の値に置き換えます)。",
              "  The steps after that are as follows (replace <...> with your own values)."),
        ]
        submit_cmd = L("qsub submit.sh     (PBS のとき)   /   sbatch submit.sh   (Slurm のとき)",
                       "qsub submit.sh     (PBS)   /   sbatch submit.sh   (Slurm)")
        status_cmd = L("qstat -u $USER     (PBS のとき)   /   squeue -u $USER    (Slurm のとき)",
                       "qstat -u $USER     (PBS)   /   squeue -u $USER    (Slurm)")
    else:
        lines += [L("  ADIT はジョブを預けません。次の手順を自分で行います (<...> は自分の値に置き換えます)。",
                    "  ADIT does not submit jobs; do the following yourself (replace <...> with your own values).")]
        submit_cmd, status_cmd = f"{profile.submit} submit.sh", profile.status
    lines += [
        L("  1. このディレクトリごとクラスタへ写します (手元の PC のターミナルで。scp / rsync = ネットワーク越しにファイルを写すコマンド)",
          "  1. Copy this whole directory to the cluster (in a terminal on your PC; scp / rsync = commands that copy files over the network)"),
        L("     <この計算ディレクトリのパス> は、手元でいま置いている場所に置き換えてください。",
          "     Replace <path to this calculation directory> with its current location on your PC."),
        f"       scp -r {local_dir} {host}:{work}/",
        L("     または", "     or"),
        f"       rsync -av {local_dir} {host}:{work}/",
        L("  2. クラスタにログインして、写したディレクトリへ移動します (ssh = 別の計算機にログインするコマンド)",
          "  2. Log in to the cluster and move into the copied directory (ssh = the command that logs in to another computer)"),
        f"       ssh {host}",
        f"       cd {work}/{dir_name}",
        L("  3. ジョブを預けます (投入)。ジョブ番号が表示されれば受け付けられています",
          "  3. Submit the job. If a job number is printed, it was accepted"),
        f"       {submit_cmd}",
        L("  4. 状態を確かめます ($USER は自分のユーザー名に置き換わります)。一覧から消えたら終わっています",
          "  4. Check the status ($USER expands to your user name). When the job disappears from the list, it has finished"),
        f"       {status_cmd}",
    ]
    if kind != "direct":
        joblog = L("<ジョブ名>.o<ジョブ番号>", "<job name>.o<job number>") if kind == "pbs" else "slurm-<job number>.out"
        lines += [
            L(f"     計算の記録は output.log、submit.sh 自体のエラー (module が読めない等) は {joblog} に書かれます",
              f"     the calculation log is output.log; errors of submit.sh itself (e.g. a module that fails to load) go to {joblog}"),
            L(f"     リソース: ノード数 {r.nodes}、ノードあたりのコア数 {r.ncpus}、MPI プロセス数 {r.mpiprocs}、OpenMP スレッド数 {r.omp_threads}、制限時間 {r.walltime}",
              f"     Resources: nodes {r.nodes}, cores/node {r.ncpus}, MPI processes/node {r.mpiprocs}, OpenMP threads {r.omp_threads}, walltime {r.walltime}"),
        ]
    lines += [
        L("  5. 終わったら結果を手元に取り戻します (手元の PC のターミナルで)", "  5. When it has finished, copy the results back (in a terminal on your PC)"),
        f"       scp -r {host}:{work}/{dir_name} {back}/",
        L("     または", "     or"),
        f"       rsync -av {host}:{work}/{dir_name} {back}/",
        "",
    ]

    analyze = (L("python analyze.py --rdf --msd   (動径分布関数と拡散係数も求めます)", "python analyze.py --rdf --msd   (also radial distribution functions and diffusion)")
               if t.type == "molecular_dynamics" else "python analyze.py")
    lines += [L("== 結果の見方 ==", "== Reading the results =="), *notes.outputs]
    if gen.supports_analysis:
        lines += [L("  analyze.py    ADIT の解析タブと同じ処理で、図と要約を analysis/ に書きます。ADIT が入った Python で実行します",
                    "  analyze.py    writes figures and a summary to analysis/, the same as the ADIT analysis tab; run it with the Python that has ADIT installed"),
                  f"       {analyze}"]
    lines.append("")
    if prov:
        from adit.provenance import readme_lines
        lines += readme_lines(prov)
    return "\n".join(lines)
