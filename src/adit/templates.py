"""Lab templates: a spec.json without the structure. List, load and save them."""

from __future__ import annotations

from adit.errors import AditValueError
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from adit.lang import L
from adit.spec import CalculationSpec, Structure

SUFFIX = ".json"
_SKIP = ("structure", "meta", "template", "handoff")


class TemplateError(AditValueError):
    pass


@dataclass(frozen=True)
class TemplateInfo:
    name: str
    path: Path
    code: str
    comment: str
    error: str = ""


def template_dirs(cfg=None) -> list[Path]:
    from adit.config import config_path

    out = []
    if cfg is not None and getattr(cfg, "templates_dir", ""):
        out.append(Path(cfg.templates_dir).expanduser())
    out.append(config_path().parent / "templates")
    return out


def _read(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as ex:
        raise TemplateError(L(f"雛形 {path} を読めません: {ex}", f"cannot read the template {path}: {ex}")) from ex
    if not isinstance(data, dict) or not isinstance(data.get("method"), dict) or "code" not in data["method"]:
        raise TemplateError(L(f"雛形 {path.name} に method (計算コードと手法の設定) がありません", f"the template {path.name} has no method (code and its settings)"))
    return data


def list_templates(cfg=None) -> list[TemplateInfo]:
    out, seen = [], set()
    for d in template_dirs(cfg):
        if not d.is_dir():
            continue
        for p in sorted(d.glob(f"*{SUFFIX}")):
            if p.stem in seen:
                continue
            seen.add(p.stem)
            try:
                data = _read(p)
                out.append(TemplateInfo(p.stem, p, data["method"]["code"], str((data.get("template") or {}).get("comment", ""))))
            except TemplateError as ex:
                out.append(TemplateInfo(p.stem, p, "", "", str(ex)))
    return out


def find_template(name_or_path: str | Path, cfg=None) -> Path:
    p = Path(name_or_path).expanduser()
    if p.suffix == SUFFIX and p.is_file():
        return p
    for info in list_templates(cfg):
        if info.name == str(name_or_path):
            return info.path
    names = ", ".join(i.name for i in list_templates(cfg)) or L("(なし)", "(none)")
    raise TemplateError(L(f"雛形 {name_or_path!r} が見つかりません (置き場所: {', '.join(str(d) for d in template_dirs(cfg))}。ある雛形: {names})",
                          f"template {name_or_path!r} not found (folders: {', '.join(str(d) for d in template_dirs(cfg))}; available: {names})"))


def flatten(d, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict) and v:
            out.update(flatten(v, key))
        else:
            out[key] = v
    return out


def load_template(name_or_path: str | Path, structure: Structure | CalculationSpec, cfg=None) -> CalculationSpec:
    from pydantic import ValidationError as PydanticError

    from adit.validate_types import friendly_pydantic

    path = find_template(name_or_path, cfg)
    raw = path.read_bytes()
    data = _read(path)
    body = {k: v for k, v in data.items() if k not in _SKIP}
    body = CalculationSpec.migrate(body) if "version" in body else dict(body, version=2)
    values = {k: v for k, v in flatten({k: v for k, v in body.items() if k != "version"}).items()}
    st = structure.structure if isinstance(structure, CalculationSpec) else structure
    meta = {"comment": str((data.get("template") or {}).get("comment", "")),
            "template": {"name": path.stem, "file": str(path), "sha256": hashlib.sha256(raw).hexdigest(), "values": values}}
    try:
        return CalculationSpec.model_validate(dict(body, structure=st.model_dump(mode="json"), meta=meta))
    except PydanticError as ex:
        raise TemplateError(L(f"雛形 {path.name} の値を読めません: ", f"cannot read the values of the template {path.name}: ") + friendly_pydantic(ex)) from ex


def template_origin(spec: CalculationSpec) -> dict[str, str]:
    t = spec.meta.template
    if not t:
        return {}
    now = flatten(spec.model_dump(mode="json"))
    return {k: t["name"] for k, v in t.get("values", {}).items() if k in now and now[k] == v}


def save_template(spec: CalculationSpec, name: str, *, directory: Path | str | None = None, cfg=None, comment: str = "",
                  overwrite: bool = False) -> Path:
    if not re.fullmatch(r"[0-9A-Za-z][0-9A-Za-z._-]*", name or ""):
        raise TemplateError(L(f"雛形の名前は英数字と . _ - だけにしてください (英数字で始める): {name!r}", f"template names may use letters, digits, . _ - (starting with a letter or digit): {name!r}"))
    d = Path(directory).expanduser() if directory else template_dirs(cfg)[0]
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{name}{SUFFIX}"
    if p.exists() and not overwrite:
        raise TemplateError(L(f"雛形 {p} はすでにあります (上書きしない)", f"the template {p} already exists (not overwritten)"))
    from adit import __version__

    data = spec.model_dump(mode="json")
    out = {"template": {"comment": comment, "saved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "app_version": __version__}}
    out.update({k: v for k, v in data.items() if k not in _SKIP})
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return p


__all__ = ["TemplateError", "TemplateInfo", "template_dirs", "list_templates", "find_template", "load_template",
           "save_template", "template_origin", "flatten"]
