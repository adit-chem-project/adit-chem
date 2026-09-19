"""Round trips, restore paths and swallowed errors in the desktop GUI (2026-09-19)."""

import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.panels.method_panel import CODES  # noqa: E402
from adit.spec import AtomsData, HubbardU, KPoints, Structure  # noqa: E402
from adit.textparse import short_number  # noqa: E402
from tests.test_gui import make_window  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def quiet(monkeypatch):
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))


def _with_code(win: MainWindow, code: str):
    win.method.code.setCurrentIndex(list(CODES).index(code))
    if code in ("vasp", "espresso"):
        win.structure.set_source("bulk"); win.structure._rebuild()
    win.refresh_preview()
    return win.current_spec()


def _spin(app, ms: int) -> None:
    from PySide6.QtCore import QEventLoop, QTimer
    loop = QEventLoop(); QTimer.singleShot(ms, loop.quit); loop.exec()


# ---- numbers -----------------------------------------------------------------

def test_short_number_keeps_every_digit():
    assert [short_number(v) for v in (0.1, 10.0, 1e-5, 1.2345678901234567, 3, True)] == \
        ["0.1", "10", "1e-05", "1.2345678901234567", "3", "True"]
    assert float(short_number(0.30000000000000004)) == 0.30000000000000004


def test_vasp_magmom_and_hubbard_values_survive_the_round_trip(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    spec = _with_code(win, "vasp")
    method = spec.method.model_copy(update={"magmom": [1.2345678901, -0.5], "hubbard": {"Cu": HubbardU(orbital="3d", u_ev=4.1234567891, j_ev=0.5)}})
    win.apply_spec(spec.model_copy(update={"method": method}))
    got = win.current_spec().method
    assert got.magmom == [1.2345678901, -0.5] and got.hubbard["Cu"].u_ev == 4.1234567891


# ---- values that the panels have no widget for --------------------------------

def test_vasp_values_outside_the_lists_are_added_not_replaced(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    spec = _with_code(win, "vasp")
    method = spec.method.model_copy(update={"ibrion": 5, "ivdw": 13, "prec": "Fast", "lreal": ".TRUE.", "algo": "Exact"})
    win.apply_spec(spec.model_copy(update={"method": method}))
    got = win.current_spec().method
    assert (got.ibrion, got.ivdw, got.prec, got.lreal, got.algo) == (5, 13, "Fast", ".TRUE.", "Exact")


def test_espresso_dipole_fields_survive_the_round_trip(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    spec = _with_code(win, "espresso")
    method = spec.method.model_copy(update={"dipole_correction": True, "dipole_direction": 3, "dipole_maxpos": 0.7,
                                            "dipole_decrease": 0.1, "dipole_amplitude": 0.01,
                                            "hubbard": {"Cu": HubbardU(orbital="3d", u_ev=4.0, j_ev=0.8)}})
    win.apply_spec(spec.model_copy(update={"method": method}))
    got = win.current_spec().method
    assert got.dipole_correction and got.dipole_direction == 3 and got.dipole_maxpos == 0.7
    assert got.dipole_decrease == 0.1 and got.dipole_amplitude == 0.01
    assert got.hubbard["Cu"].j_ev == 0.8                          # J is not shown for QE but must not turn into 0


def test_lammps_pressure_tensor_and_dftb_seed_survive_the_round_trip(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    spec = win.current_spec()
    win.apply_spec(spec.model_copy(update={"method": spec.method.model_copy(update={"seed": 777})}))
    assert win.current_spec().method.seed == 777
    spec = _with_code(win, "lammps")
    win.apply_spec(spec.model_copy(update={"method": spec.method.model_copy(update={"thermo_pressure_tensor": True})}))
    assert win.current_spec().method.thermo_pressure_tensor is True


def test_kpoint_shift_keeps_three_components_until_the_field_is_changed(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("bulk"); win.structure._rebuild(); win.refresh_preview()
    win.kpoints.set_kpoints(KPoints(mode="mesh", mesh=(2, 2, 2), shift=(0.5, 0.0, 0.0)))
    assert win.kpoints.kpoints().shift == (0.5, 0.0, 0.0)
    win.kpoints.shift.setCurrentText("0")
    assert win.kpoints.kpoints().shift == (0.0, 0.0, 0.0)
    spec = win.current_spec()
    win.apply_spec(spec.model_copy(update={"kpoints": None}))     # a spec without k-points resets the panel
    assert win.kpoints.kpoints().mode == "gamma"


def test_unknown_profile_and_sk_set_are_reported(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    spec = win.current_spec()
    bad = spec.model_copy(update={"runtime": spec.runtime.model_copy(update={"profile": "nowhere"}),
                                  "method": spec.method.model_copy(update={"sk_set": "not-a-set"})})
    win.apply_spec(bad)
    msg = win.statusBar().currentMessage()
    assert "nowhere" in msg and "not-a-set" in msg and "変わりました" in msg


def test_reload_config_reaches_the_code_panels(app, quiet, sk_root, tmp_path):
    from adit.config import save_config
    win = make_window(sk_root, tmp_path)
    save_config(win.cfg, tmp_path / "cluster.toml")
    win.reload_config()
    assert win.method.vasp.cfg is win.cfg and win.method.espresso.cfg is win.cfg


# ---- restoring a structure ------------------------------------------------------

def test_a_spacegroup_spec_is_refused_with_a_message(app, quiet, sk_root, tmp_path, monkeypatch):
    shown = []
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: shown.append(a[2])))
    win = make_window(sk_root, tmp_path)
    spec = win.current_spec()
    win.apply_spec(spec.model_copy(update={"structure": spec.structure.model_copy(update={"source": "spacegroup"})}))
    assert shown and "spacegroup" in shown[0]


def test_smiles_coordinates_from_the_spec_are_not_rebuilt(app, quiet, sk_root, tmp_path):
    pytest.importorskip("rdkit")
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("smiles"); win.structure.smiles.setText("O"); win.structure._rebuild()
    spec = win.current_spec()
    moved = np.asarray(spec.structure.atoms.positions) + 0.37
    st = spec.structure.model_copy(update={"atoms": spec.structure.atoms.model_copy(update={"positions": [tuple(p) for p in moved]})})
    win.apply_spec(spec.model_copy(update={"structure": st}))
    _spin(app, 800)                                               # the 500 ms SMILES timer must not fire
    assert np.allclose(win.structure.structure().atoms.positions, moved)


def test_a_molecule_in_a_box_comes_back_in_its_box(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.box_size.setValue(12.0); win.structure.box.setChecked(True)
    spec = win.current_spec()
    assert spec.structure.periodic
    win2 = make_window(sk_root, tmp_path / "b")
    win2.apply_spec(spec)
    assert win2.structure.box.isChecked() and win2.structure.box_size.value() == 12.0
    win2.structure.charge.setValue(0); win2.structure._rebuild()  # touching a field keeps it periodic
    assert win2.current_spec().structure.periodic


def test_a_leftover_preset_filter_does_not_change_the_molecule(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.preset.setCurrentIndex(win.structure.preset.findData("CH4")); win.structure._rebuild()
    spec = win.current_spec()
    win2 = make_window(sk_root, tmp_path / "b")
    win2.structure.preset_search.setText("H2O")
    win2.apply_spec(spec)
    assert win2.current_spec().structure.source_ref == "CH4"
    assert win2.structure.preset_search.text() == ""


# ---- other swallowed problems ---------------------------------------------------

def test_running_the_analysis_keeps_the_rdf_and_msd_choices(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    (tmp_path / "run").mkdir()
    win.analysis.run_dir.setText(str(tmp_path / "run"))
    win.analysis.cb_rdf.setChecked(True); win.analysis.cb_msd.setChecked(True)
    win.analysis.run()
    assert win.analysis.cb_rdf.isChecked() and win.analysis.cb_msd.isChecked()


def test_a_failed_image_save_is_reported(app, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog, QLabel

    from adit.gui.copy_save import save_dialog

    warned = []
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: warned.append(a[2])))
    target = tmp_path / "missing_dir" / "shot.png"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(target), "")))
    assert save_dialog(None, "shot.png", widget=QLabel("x")) is None
    assert warned and str(target) in warned[0]


def test_the_3d_viewer_draws_heavy_elements(app):
    from ase import Atoms

    from adit.gui.copy_save import widget_image
    from adit.gui.viewer3d import Viewer3D

    view = Viewer3D()
    view.resize(200, 200)
    view.set_atoms(Atoms(numbers=[110, 118, 1], positions=[(0, 0, 0), (2, 0, 0), (0, 2, 0)]))
    assert not widget_image(view).isNull()                             # no IndexError past the color table


def test_settings_keep_the_trailing_comment_when_a_value_changes():
    from adit.gui.settings_dialog import set_top_level

    text = 'theme = "auto"  # 好み\nsk_root = "/a"\n'
    out = set_top_level(text, {"theme": "dark"})
    assert 'theme = "dark"  # 好み' in out and 'sk_root = "/a"' in out


def test_the_report_refuses_to_overwrite_files(tmp_path):
    from adit.gui import report_fields as R
    from adit.report import ReportError

    existing = tmp_path / "methods.md"; existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ReportError, match="methods.md"):
        R.build(R.ReportRequest(run_dirs=[], language="", methods_path=existing, conditions_csv=None, results_csv=None, bundle=None, check=False))
    assert existing.read_text(encoding="utf-8") == "keep"


# ---- English mode -------------------------------------------------------------------

def test_help_is_found_in_english_for_rows_that_used_translated_keys(app):
    from adit.gui.i18n import set_language, translate_widgets
    from adit.gui.panels.cp2k_panel import Cp2kMethodPanel
    from adit.gui.panels.method_panel import DftbMethodPanel
    from adit.gui.panels.task_panel import TaskPanel
    from adit.config import default_config
    from PySide6.QtWidgets import QLabel

    set_language("en")
    try:
        panels = [TaskPanel(), DftbMethodPanel(""), Cp2kMethodPanel(default_config())]
        for p in panels:
            translate_widgets(p)
        labels = {w.property("adit_key"): w for p in panels for w in p.findChildren(QLabel) if w.property("adit_key")}
        for key, shown in (("最大ステップ数 (MaxSteps)", "Maximum steps (MaxSteps)"),
                           ("SCC の収束判定 (SccTolerance)", "SCC tolerance (SccTolerance)"),
                           ("SCC の反復上限 (MaxSccIterations)", "Maximum SCC iterations (MaxSccIterations)"),
                           ("軌道変換法 (OT)", "Orbital transformation (OT)"), ("OT の追加行", "Additional OT lines")):
            assert labels[key].text() == shown and labels[key].toolTip(), key
    finally:
        set_language("ja")


def test_batch_and_report_dialogs_are_translated(app, sk_root, tmp_path):
    from PySide6.QtWidgets import QLabel

    from adit.gui.i18n import set_language
    from adit.gui.prep23_dialogs import DIALOGS
    from adit.gui.report_dialog import ReportDialog
    from adit.gui import report_fields as R
    from adit.web import prep23 as P

    set_language("en")
    try:
        dlg = DIALOGS["compare"](lambda: None, None, "")
        texts = [w.text() for w in dlg.findChildren(QLabel)]
        assert P.LABELS["out"][1] in texts and P.LABELS["out"][0] not in texts
        rep = ReportDialog()
        texts = [w.text() for w in rep.findChildren(QLabel)]
        assert R.LABELS["rep_dirs"][1] in texts and R.LABELS["rep_dirs"][0] not in texts
        assert rep.check.toolTip()
    finally:
        set_language("ja")


def test_every_msd_and_vanhove_row_has_help():
    from adit.gui import analysis_fields as AF
    from adit.gui.help import help_for

    for key in ("msd_axes", "msd_keep_drift", "msd_blocks", "msd_per_atom", "vanhove", "vanhove_taus", "vanhove_displacement"):
        h = help_for(AF.LABELS[key][0])
        assert h is not None and h.ja and h.en, key
