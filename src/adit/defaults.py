"""Which fields differ from the code's own defaults (the pydantic defaults in spec.py), and how to put one back."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from pydantic import BaseModel
from pydantic_core import PydanticUndefined

from adit.lang import L
from adit.spec import CalculationSpec

# Not marked: the code is a choice, not a value with a default; runtime values are resource requests the GUI seeds itself.
UNMARKED_PREFIXES = ("method.code", "runtime.", "structure.", "meta.", "handoff", "plumed", "version")


@dataclass(frozen=True)
class Change:
    label: str
    paths: tuple[str, ...]
    value: Any
    default: Any


def _field_default(model: BaseModel | type[BaseModel], name: str) -> tuple[bool, Any]:
    fields = (model if isinstance(model, type) else type(model)).model_fields
    f = fields.get(name)
    if f is None:
        return False, None
    if f.default is not PydanticUndefined:
        return True, f.default
    if f.default_factory is not None:
        return True, f.default_factory()
    return False, None


def default_at(spec: CalculationSpec, path: str) -> tuple[bool, Any]:
    """(known, default) for a dotted path such as task.md.steps; the method's own class decides method.* defaults."""
    parts = path.split(".")
    model: Any = spec
    for i, name in enumerate(parts):
        if isinstance(model, BaseModel):
            known, default = _field_default(model, name)
            if i == len(parts) - 1:
                return known, default
            model = getattr(model, name, None)
            if model is None:
                return False, None
        else:
            return False, None
    return False, None


MISSING = object()


def value_at(spec: CalculationSpec, path: str) -> Any:
    """The value at a dotted path, or MISSING when a parent along the way is absent (kpoints of a molecule)."""
    v: Any = spec
    for name in path.split("."):
        if isinstance(v, BaseModel):
            v = getattr(v, name, MISSING)
        elif isinstance(v, dict):
            v = v.get(name, MISSING)
        else:
            return MISSING
    return v


def same(a: Any, b: Any) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        try:
            return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=0.0)
        except (TypeError, ValueError):
            return False
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, BaseModel) and isinstance(b, BaseModel):
        return a.model_dump() == b.model_dump()
    return a == b


def _paths_by_label(spec: CalculationSpec) -> dict[str, list[str]]:
    from adit.web.codefields import COMMON_LABELS, METHOD_LABELS, label_for_path

    code = spec.method.code
    out: dict[str, list[str]] = {}
    for path in [*COMMON_LABELS, *(f"method.{name}" for name in METHOD_LABELS.get(code, {}))]:
        if path.startswith(UNMARKED_PREFIXES):
            continue
        lab = label_for_path(code, path)
        if lab:
            out.setdefault(lab, []).append(path)
    return out


def changed_fields(spec: CalculationSpec) -> dict[str, Change]:
    """Labels whose value differs from the spec default, in form order of the label tables."""
    out: dict[str, Change] = {}
    for lab, paths in _paths_by_label(spec).items():
        for path in paths:
            known, default = default_at(spec, path)
            if not known:
                continue
            value = value_at(spec, path)
            if value is not MISSING and not same(value, default):
                out[lab] = Change(lab, tuple(paths), value, default)
                break
    return out


def format_value(v: Any) -> str:
    if v is None:
        return L("なし", "none")
    if isinstance(v, bool):
        return L("オン", "on") if v else L("オフ", "off")
    if isinstance(v, float):
        return f"{v:g}"
    if isinstance(v, str):
        return v if v else L("空欄", "empty")
    if isinstance(v, (list, tuple)):
        return " ".join(format_value(x) for x in v) if v else L("なし", "none")
    if isinstance(v, dict):
        return ", ".join(f"{k} = {format_value(x)}" for k, x in v.items()) if v else L("なし", "none")
    if isinstance(v, BaseModel):
        return format_value(v.model_dump())
    return str(v)


def with_default(model: BaseModel, parts: list[str]) -> BaseModel:
    """A copy of the model with the field at parts (relative to it) put back to its default."""
    name = parts[0]
    if len(parts) == 1:
        known, default = _field_default(model, name)
        return model.model_copy(update={name: default}) if known else model
    sub = getattr(model, name)
    return model.model_copy(update={name: with_default(sub, parts[1:])})


__all__ = ["Change", "MISSING", "changed_fields", "default_at", "value_at", "same", "format_value", "with_default", "UNMARKED_PREFIXES"]
