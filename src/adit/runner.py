
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from adit.lang import L

SUBMIT = "submit.sh"


@dataclass
class RunTarget:
    run_dir: Path
    kind: str = ""
    exe: str = ""
    code: str = ""


def executable_of(run_command: str) -> str:
    tail = run_command.split("&&")[-1].strip()
    words = tail.split()
    if not words:
        return ""
    first = words[0]
    if Path(first).name in ("mpirun", "mpiexec", "srun") and len(words) > 1:
        for word in words[1:]:
            if not word.startswith("-") and not word.isdigit():
                return word
    return first


def target_from_dir(run_dir: Path | str, cfg) -> RunTarget:
    from adit.codes import GENERATORS
    from adit.spec import CalculationSpec

    run_dir = Path(run_dir).expanduser()
    spec = CalculationSpec.load(run_dir / "spec.json")
    name = spec.runtime.profile
    if name not in cfg.profiles:
        return RunTarget(run_dir, "", "", spec.method.code)
    profile = cfg.profiles[name]
    generator = GENERATORS.get(spec.method.code)
    command = generator.run_command(spec, profile) if generator is not None else ""
    return RunTarget(run_dir, profile.kind, executable_of(command), spec.method.code)


def block_reason(cfg, target: RunTarget | None, *, running: bool = False, cfg_path=None) -> str:
    if not cfg.enable_run:
        where = cfg_path or getattr(cfg, "source_path", None) or "~/.config/adit/cluster.toml"
        return L(f"この PC では実行しない設定です ({where} の enable_run = true にすると実行できます)。"
                 "ADIT はクラスタへの投入はしません",
                 f"running here is switched off (set enable_run = true in {where} to enable it); "
                 "ADIT never submits to a cluster")
    if running:
        return L("実行中", "running")
    if target is None:
        return L("先に「生成」を押してください", "Press Generate first")
    if not (target.run_dir / SUBMIT).is_file():
        return L(f"{target.run_dir} に {SUBMIT} がありません", f"{target.run_dir} has no {SUBMIT}")
    if target.kind != "direct":
        return L("クラスタ用 (PBS / Slurm) の入力は、ADIT からは実行しません",
                 "Cluster (PBS / Slurm) inputs are not run by ADIT")
    if os.name == "nt":
        return L("Windows では実行できません (生成した入力を Linux のサーバーに転送して使います)",
                 "Not run on Windows (transfer the files to a Linux server)")
    if not target.exe:
        return L("実行コマンドが分かりません", "cannot determine the executable")
    if "/" in target.exe or os.sep in target.exe:
        if not Path(target.exe).is_file():
            return L(f"実行ファイルがありません: {target.exe}", f"executable not found: {target.exe}")
    elif shutil.which(target.exe) is None:
        return L(f"{target.exe} が見つかりません (インストールされていないか、コマンドを探す場所 PATH に入っていません。"
                 "README.txt の「実行する前に用意すること」を参照)",
                 f"{target.exe} not found (not installed, or not on PATH, the list of places where commands are looked up; "
                 "see \"Before running\" in README.txt)")
    return ""


def start(target: RunTarget, log_name: str = "run.log") -> subprocess.Popen:
    log = (target.run_dir / log_name).open("wb")
    return subprocess.Popen(["bash", SUBMIT], cwd=str(target.run_dir), stdout=log, stderr=subprocess.STDOUT)


def run_and_wait(target: RunTarget, on_line=None, log_name: str = "run.log") -> int:
    with (target.run_dir / log_name).open("w", encoding="utf-8", errors="replace") as log:
        proc = subprocess.Popen(["bash", SUBMIT], cwd=str(target.run_dir), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, errors="replace", bufsize=1)
        for line in proc.stdout:
            log.write(line)
            if on_line is not None:
                on_line(line.rstrip("\n"))
        proc.stdout.close()
        return proc.wait()
