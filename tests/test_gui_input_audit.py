"""Desktop entry points for the shared native-input and cross-run checks."""

import os
from types import SimpleNamespace

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from adit.config import default_config  # noqa: E402
from adit.gui.input_audit_dialog import InputAuditDialog  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.native_verify import NativeVerificationResult  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_file_ribbon_exposes_inspection(app, tmp_path):
    window = MainWindow(default_config(sk_root=str(tmp_path)), tmp_path / "cluster.toml")
    assert window.act_input_audit in window.ribbon.actions_on_page(0)
    assert "点検" in window.act_input_audit.iconText()


def test_import_uses_shared_review_writer_and_shows_draft_status(app, monkeypatch, tmp_path):
    from adit import native_workflow

    calls = []
    result = SimpleNamespace(spec=object(), report=lambda: {"complete": True, "ready_to_regenerate": False})
    output = tmp_path / "review"

    def write(source, destination, *, code=None):
        calls.append((source, destination, code))
        return result, output

    monkeypatch.setattr(native_workflow, "write_import_review", write)
    dialog = InputAuditDialog()
    dialog.import_source.setText(str(tmp_path / "native"))
    dialog.import_review.setText(str(output))
    dialog.import_code.setCurrentIndex(dialog.import_code.findData("vasp"))
    dialog.inspect_button.click()
    assert calls == [(str(tmp_path / "native"), str(output), "vasp")]
    assert "import_report.json" in dialog.summary.text()
    assert "draft_spec.json" not in dialog.details.toPlainText()
    assert '"ready_to_regenerate": false' in dialog.details.toPlainText()


def test_verify_displays_mismatch_and_unverified_without_writing(app, monkeypatch, tmp_path):
    from adit import native_verify

    result = NativeVerificationResult(
        preserved=[{"field": "method.encut"}], mismatched=[{"field": "kpoints.mesh"}],
        unverifiable=[{"field": "runtime.profile"}], native_report={"complete": True})
    monkeypatch.setattr(native_verify, "verify_generated_bundle", lambda path: result)
    dialog = InputAuditDialog()
    dialog.kind.setCurrentIndex(dialog.kind.findData("verify"))
    dialog.verify_dir.setText(str(tmp_path / "generated"))
    dialog.inspect()
    assert "問題あり" in dialog.summary.text()
    assert "kpoints.mesh" in dialog.details.toPlainText()
    assert "runtime.profile" in dialog.details.toPlainText()
    assert list(tmp_path.iterdir()) == []


def test_cross_run_audit_uses_selected_mode_and_shows_scope(app, monkeypatch, tmp_path):
    from adit.analysis import comparability

    calls = []

    def audit(paths, *, require_md):
        calls.append((paths, require_md))
        return SimpleNamespace(summary_text=lambda: "energy audit: failed\nnote: chemical comparability is not judged")

    monkeypatch.setattr(comparability, "audit_run_dirs", audit)
    dialog = InputAuditDialog()
    dialog.kind.setCurrentIndex(dialog.kind.findData("audit"))
    dialog.audit_first.setText(str(tmp_path / "run_a"))
    dialog.audit_second.setText(str(tmp_path / "run_b"))
    dialog.audit_kind.setCurrentIndex(dialog.audit_kind.findData("energy"))
    dialog.inspect()
    assert calls == [([str(tmp_path / "run_a"), str(tmp_path / "run_b")], False)]
    assert "energy audit: failed" in dialog.details.toPlainText()
    assert "化学的に比べてよいか" in dialog.summary.text()


def test_switching_kind_hides_stale_results(app):
    dialog = InputAuditDialog()
    dialog.summary.setText("stale")
    dialog.details.setPlainText("stale")
    dialog.kind.setCurrentIndex(dialog.kind.findData("verify"))
    assert not dialog.summary.text()
    assert not dialog.details.toPlainText()


def test_dialog_labels_are_available_in_english(app):
    from adit.gui.i18n import set_language

    set_language("en")
    try:
        dialog = InputAuditDialog()
        assert dialog.windowTitle() == "Inspect inputs and results"
        assert dialog.kind.itemText(0) == "Import existing input"
        assert dialog.kind.itemText(2) == "Compare prerequisites across runs"
        assert dialog.audit_kind.itemText(0) == "MSD prerequisites"
    finally:
        set_language("ja")
