"""Shared, conservative workflow for inspecting native input bundles."""
# The CLI, desktop GUI, and web UI call this module; it never executes input or
# submits a calculation. Parsing is delegated to :mod:`adit.native_import`.

from __future__ import annotations

from adit.errors import AditValueError
import json
from pathlib import Path

from adit.lang import L
from adit.native_import import NativeImportError, NativeImportResult, import_native
from adit.native_verify import NativeVerificationResult


class NativeWorkflowError(AditValueError):
    """An audit destination is unsafe or cannot be created."""


def write_import_review(source: Path | str, output: Path | str, *,
                        code: str | None = None) -> tuple[NativeImportResult, Path]:
    """Write a new report and, only for fully parsed input, a review draft."""
    # A parser failure still produces a report. The original input is never
    # changed, and no review artifact is written inside its bundle.
    if not str(source).strip() or not str(output).strip():
        raise NativeWorkflowError(L("入力と点検結果の保存先を指定してください",
                                    "Specify both the input and a review-output directory."))
    src = Path(source).expanduser().resolve()
    dst = Path(output).expanduser().resolve()
    bundle = src if src.is_dir() else src.parent if src.is_file() else None
    if dst.exists() or dst == src or dst.is_relative_to(src) or (bundle is not None and dst.is_relative_to(bundle)):
        raise NativeWorkflowError(L(
            f"点検結果の保存先は、既存入力の外にある新しいディレクトリにしてください: {dst}",
            f"Choose a new review-output directory outside the existing input bundle: {dst}"))
    try:
        result = import_native(src, code=code)
    except NativeImportError as exc:
        result = NativeImportResult(code or "undetermined")
        result.issue(src, 1, src.name, str(exc))
    dst.mkdir(parents=True, exist_ok=False)
    (dst / "import_report.json").write_text(
        json.dumps(result.report(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if result.spec is not None:
        result.spec.save(dst / "draft_spec.json")
    return result, dst


def verification_passed(result: NativeVerificationResult) -> bool:
    """Whether mapped, re-importable fields passed; not global equivalence."""
    return bool(result.preserved) and not result.mismatched and bool(result.native_report.get("complete"))


__all__ = ["NativeWorkflowError", "write_import_review", "verification_passed"]
