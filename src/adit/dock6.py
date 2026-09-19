"""Package an existing, user-authored DOCK 6 input without choosing chemistry."""
# DOCK 6.13 input/command and rigid/flexible ligand tutorial:
# https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm
# https://dock.compbio.ucsf.edu/DOCK_6/tutorials/ligand_sampling_dock/ligand_sampling_dock.html

from __future__ import annotations

from adit.errors import AditValueError
import hashlib
import json
import shutil
from pathlib import Path

from adit.lang import L


class Dock6Error(AditValueError):
    pass


_INPUT_FILES = ("ligand_atom_file", "receptor_site_file", "vdw_defn_file", "flex_defn_file", "flex_drive_file")


def _inside(source_dir: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or " " in value or not value:
        raise Dock6Error(L(f"DOCK6 の入力ファイルは作業ディレクトリ内の相対パスで指定してください: {value!r}",
                           f"DOCK6 input files must use paths relative to the working directory: {value!r}"))
    candidate = source_dir / path
    if not candidate.is_file() or not candidate.resolve().is_relative_to(source_dir.resolve()):
        raise Dock6Error(L(f"DOCK6 の入力ファイルがありません: {value}", f"DOCK6 input file not found: {value}"))
    return candidate


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_dock6(source: Path | str, output: Path | str, *, assets: list[str] | None = None) -> Path:
    """Stage dock.in and its declared local inputs, leaving the input deck unchanged."""
    # This does not prepare receptors, assign charges, or validate all DOCK6
    # keywords. The author must provide a complete dock.in and its dependencies.
    src, dst = Path(source).expanduser(), Path(output).expanduser()
    if not src.is_file() or src.name != "dock.in":
        raise Dock6Error(L("入力には既存の dock.in を指定してください", "provide an existing dock.in input file"))
    root = src.parent.resolve()
    if dst.resolve() == root or dst.resolve().is_relative_to(root):
        raise Dock6Error(L("出力ディレクトリは入力ディレクトリの外にしてください", "put the output directory outside the input directory"))
    if dst.exists() and (not dst.is_dir() or any(dst.iterdir())):
        raise Dock6Error(L(f"出力先 {dst} にはすでにファイルがあります (上書きしません)", f"output directory {dst} already contains files (not overwritten)"))
    try:
        content = src.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as ex:
        raise Dock6Error(L(f"dock.in を読めません: {ex}", f"cannot read dock.in: {ex}")) from ex
    params = {}
    for line in content.splitlines():
        words = line.strip().split()
        if len(words) >= 2 and not words[0].startswith("#"):
            params[words[0].lower()] = words[1]
    if not params.get("ligand_atom_file"):
        raise Dock6Error(L("dock.in に ligand_atom_file がありません", "dock.in has no ligand_atom_file"))
    if params.get("orient_ligand", "yes").lower() == "yes" and not params.get("receptor_site_file"):
        raise Dock6Error(L("orient_ligand=yes には receptor_site_file が必要です", "orient_ligand=yes requires receptor_site_file"))
    uses_grid = any(params.get(key, "no").lower() == "yes" for key in ("grid_score_primary", "grid_score_secondary"))
    if uses_grid and not params.get("grid_score_grid_prefix"):
        raise Dock6Error(L("grid_score_primary/secondary=yes には grid_score_grid_prefix が必要です", "grid_score_primary/secondary=yes requires grid_score_grid_prefix"))
    paths = {"dock.in": src}
    for key in _INPUT_FILES:
        if value := params.get(key):
            paths[value] = _inside(root, value)
    if uses_grid:
        prefix = params["grid_score_grid_prefix"]
        for ext in (".nrg", ".bmp"):
            value = prefix + ext
            paths[value] = _inside(root, value)
    for value in assets or []:
        paths[value] = _inside(root, value)
    dst.mkdir(parents=True, exist_ok=True)
    for name, path in paths.items():
        target = dst / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    manifest = {name: _sha256(path) for name, path in paths.items()}
    (dst / "input_sha256.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (dst / "run.sh").write_text("#!/bin/sh\nset -e\ndock6 -i dock.in -o dock.out\n", encoding="utf-8")
    (dst / "run.sh").chmod(0o755)
    (dst / "README.txt").write_text(L(
        "DOCK6 の既存入力をまとめたディレクトリです。ADIT は受容体、配位子の原子型・電荷、球、格子、スコアを作成・変更していません。\n"
        "dock.in と各入力ファイルを確認してから、./run.sh を実行するか、クラスタのジョブスクリプトから呼んでください。ADIT は投入しません。\n"
        "dock.out と配位子の .mol2 出力を確認してください。input_sha256.json に複製元のハッシュを記録しています。\n",
        "This directory packages existing DOCK6 input. ADIT did not create or alter receptor preparation, ligand atom types or charges, spheres, grids, or scoring settings.\n"
        "Review dock.in and the input files, then run ./run.sh yourself or call it from your cluster job script. ADIT does not submit jobs.\n"
        "Inspect dock.out and the ligand .mol2 output. input_sha256.json records hashes of the copied inputs.\n"), encoding="utf-8")
    return dst
