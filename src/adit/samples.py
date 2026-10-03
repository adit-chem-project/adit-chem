
from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from adit.lang import L

REPO_URL = "https://github.com/adit-chem-project/adit-chem"


@dataclass(frozen=True)
class Sample:
    name: str
    path: Path
    code: str
    task: str
    formula: str
    natoms: int
    comment: str

    def line(self) -> str:
        return f"  {self.name:<34} {self.code:<10} {self.task:<22} {self.formula:<12} {self.comment}"


def _has_samples(candidate: Path) -> bool:
    return candidate.is_dir() and any(candidate.glob("*/spec.json"))


def examples_dir() -> Path | None:
    import sys

    bundled = getattr(sys, "_MEIPASS", "")
    if bundled and _has_samples(Path(bundled) / "examples"):
        return Path(bundled) / "examples"
    here = Path(__file__).resolve()
    packaged = here.parent / "examples"
    for parent in here.parents:
        candidate = parent / "examples"
        if candidate != packaged and _has_samples(candidate):
            return candidate
    return packaged if _has_samples(packaged) else None


def list_samples() -> list[Sample]:
    root = examples_dir()
    if root is None:
        return []
    from adit.spec import CalculationSpec

    out: list[Sample] = []
    for path in sorted(root.iterdir()):
        spec_file = path / "spec.json"
        if not spec_file.is_file():
            continue
        try:
            spec = CalculationSpec.load(spec_file)
        except Exception:
            continue
        comment = ""
        sources = path / "SOURCES.md"
        if sources.is_file():
            first = sources.read_text(encoding="utf-8", errors="replace").splitlines()[0]
            comment = first.lstrip("# ").strip()
            if "——" in comment:
                comment = comment.split("——", 1)[1].strip()
        out.append(Sample(name=path.name, path=path, code=spec.method.code, task=spec.task.type,
                          formula=spec.atoms.get_chemical_formula(), natoms=len(spec.atoms), comment=comment))
    return out


def missing_message() -> str:
    return L(f"サンプルが見つかりません (パッケージの中にも、ソースの examples/ にもありません)。"
             f"インストールをやり直すか、ソースを取得してください: {REPO_URL}",
             f"No samples were found (neither inside the package nor in examples/ next to the source). "
             f"Reinstall, or get the source: {REPO_URL}")


def copy_sample(name: str, destination: Path | str, *, overwrite: bool = False) -> Path:
    samples = {s.name: s for s in list_samples()}
    if name not in samples:
        known = ", ".join(sorted(samples)) or L("(サンプルがありません)", "(no samples)")
        raise FileNotFoundError(L(f"サンプル {name!r} はありません。あるもの: {known}",
                                  f"no sample named {name!r}; available: {known}"))
    out = Path(destination).expanduser()
    if out.is_dir():
        out = out / f"{name}.json"
    if out.exists() and not overwrite:
        raise FileExistsError(L(f"すでにあります: {out} (上書きしてよければ --overwrite)",
                                f"already exists: {out} (use --overwrite to replace it)"))
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(samples[name].path / "spec.json", out)
    return out
