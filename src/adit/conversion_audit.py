"""Mechanical, field-level audit of a code-to-code CalculationSpec retarget."""
# The audit describes where settings came from. It does not assert that two
# engines use the same equations of motion or electronic-structure model.

from __future__ import annotations

from typing import Any


_UNITS = {
    "structure.atoms.positions": "Å",
    "structure.atoms.cell": "Å",
    "structure.velocities": "Å/fs",
    "task.force_tolerance_ev_per_ang": "eV/Å",
    "task.md.timestep_fs": "fs",
    "task.md.temperature_k": "K",
    "task.md.coupling_time_fs": "fs",
    "task.md.pressure_bar": "bar",
    "task.md.barostat_time_fs": "fs",
}


def _value(obj: Any, path: str) -> Any:
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj


def _json_value(value: Any) -> Any:
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if isinstance(value, tuple):
        return list(value)
    return value


def field_audit(source, target_template, converted, report: dict) -> dict:
    """Describe provenance and input applicability without chemical claims."""
    # Atom arrays are deliberately summarized to keep conversion.json bounded;
    # their exact values remain in spec.json and the generated structure files.
    from adit.lang import L

    entries: list[dict] = []
    explicit = set(report.get("changes", {}))
    omitted = set(report.get("omitted", []))
    unapplied = report.get("not_applied_by_target", {})

    def add(path: str, status: str, origin: str, *, source_value: Any = None,
            target_value: Any = None, note: str = "") -> None:
        row = {"path": path, "status": status, "origin": origin}
        if path in _UNITS:
            row["unit"] = _UNITS[path]
        if source_value is not None:
            row["source_value"] = _json_value(source_value)
        if target_value is not None:
            row["target_value"] = _json_value(target_value)
        if note:
            row["note"] = note
        entries.append(row)

    atoms = source.structure.atoms
    for part in ("symbols", "positions", "cell", "pbc"):
        path = f"structure.atoms.{part}"
        src = getattr(atoms, part)
        dst = getattr(converted.structure.atoms, part)
        if src != dst:
            raise ValueError(f"retarget changed {path} without an audit entry")
        add(path, "preserved", "source_spec", note=L(
            f"{len(src)} 件の値を保持。実値は spec.json を参照してください。" if part in {"symbols", "positions"} else "元の値を保持。実値は spec.json を参照してください。",
            f"Preserved {len(src)} values; see spec.json for the exact values." if part in {"symbols", "positions"} else "Preserved; see spec.json for the exact values."))

    for part in ("source", "source_ref", "charge", "multiplicity", "fixed_atoms", "fixed_axes", "velocities"):
        path = f"structure.{part}"
        src = _value(source, path)
        dst = _value(converted, path)
        if path in omitted:
            status, origin = "omitted_by_request", "explicit_option"
        elif path in unapplied:
            status, origin = "retained_in_spec_not_input", "source_spec"
        elif src == dst:
            status, origin = "preserved", "source_spec"
        else:
            raise ValueError(f"retarget changed {path} without an audit entry")
        add(path, status, origin, source_value=src, target_value=dst,
            note=unapplied.get(path, {}).get("reason", ""))

    task = source.task
    paths = [f"task.{p}" for p in type(task).model_fields if p != "md"]
    active_task = {
        "single_point": {"task.type"},
        "geometry_optimization": {"task.type", "task.optimizer", "task.max_steps", "task.force_tolerance_ev_per_ang", "task.relax_cell"},
        "molecular_dynamics": {"task.type"},
        "vibrations": {"task.type"},
        "band_structure": {"task.type", "task.bands"},
    }.get(task.type, {"task.type"})
    if task.type == "molecular_dynamics":
        paths += [f"task.md.{p}" for p in type(task.md).model_fields]
    for path in paths:
        src = _value(source, path)
        dst = _value(converted, path)
        if path not in active_task and not path.startswith("task.md."):
            status, origin = "not_used_by_task", "source_spec"
        elif path in explicit:
            status, origin = "changed_by_request", "target_template"
        elif path in unapplied:
            status, origin = "retained_in_spec_not_input", "source_spec"
        elif src == dst:
            status, origin = "preserved", "source_spec"
        else:
            raise ValueError(f"retarget changed {path} without an audit entry")
        add(path, status, origin, source_value=src, target_value=dst,
            note=(L("この計算種類では使いません。", "Not used by this task type.") if status == "not_used_by_task"
                  else unapplied.get(path, {}).get("reason", "")))

    for part in type(source.runtime).model_fields:
        path = f"runtime.{part}"
        src, dst = _value(source, path), _value(converted, path)
        if src != dst:
            raise ValueError(f"retarget changed {path} without an audit entry")
        add(path, "preserved", "source_spec", source_value=src, target_value=dst)

    add("method", "not_translated", "target_template", source_value=source.method.code,
        target_value=converted.method.code, note=L(
            "力場・電子状態の手法・擬ポテンシャル等は対応付けていません。変換先の雛形の全値を確認してください。",
            "Force fields, electronic-structure methods, and pseudopotentials are not mapped. Review every target-template value."))
    add("kpoints", "from_target_template", "target_template",
        source_value=source.kpoints, target_value=converted.kpoints, note=L(
            "k 点は変換せず、変換先の雛形から取っています。",
            "K-points were not translated; they came from the target template."))

    return {
        "schema_version": 1,
        "equivalence_claim": "none",
        "equivalence_note": L(
            "入力値の出典を記録します。異なる計算コードが同じ物理モデルや数値解を実現するとは主張しません。",
            "This records the origin of input settings. It does not claim that different engines implement the same physical model or produce the same numerical result."),
        "application_scope": L(
            "preserved は Spec 内の値を保持した意味です。生成入力への適用を項目ごとに検証した意味ではありません。既知の未適用項目だけを別記しています。",
            "Preserved means the value remains in the spec, not that its application to the generated input was verified field by field. Known unapplied settings are listed separately."),
        "entries": entries,
        "target_template_code": target_template.method.code,
    }


__all__ = ["field_audit"]
