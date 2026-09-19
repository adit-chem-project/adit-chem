
import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from ase.build import bulk  # noqa: E402
from ase.io import write  # noqa: E402
from PySide6.QtWidgets import QApplication, QFormLayout, QMessageBox  # noqa: E402

from adit.config import default_config, load_config, save_config  # noqa: E402
from adit.gui.help import help_for  # noqa: E402
from adit.gui.i18n import _EN  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.panels.method_panel import CODES  # noqa: E402
from tests.test_cp2k import BASIS, POTENTIALS  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
SPCE = REPO / "examples" / "gromacs_spce"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def quiet(monkeypatch):
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))


@pytest.fixture
def cp2k_data(tmp_path) -> Path:
    d = tmp_path / "cp2k_data"; d.mkdir()
    (d / "BASIS_MOLOPT").write_text(BASIS, encoding="utf-8"); (d / "GTH_POTENTIALS").write_text(POTENTIALS, encoding="utf-8")
    return d


def window(sk_root, tmp_path, cp2k_data="") -> MainWindow:
    cfg = default_config(sk_root=str(sk_root)); cfg.cp2k_data = str(cp2k_data)
    win = MainWindow(cfg, tmp_path / "cluster.toml")
    win.runtime.outdir.setText(str(tmp_path / "out"))
    return win


def use_code(win: MainWindow, code: str) -> None:
    win.method.code.setCurrentIndex(list(CODES).index(code))


def use_file(win: MainWindow, path: Path) -> None:
    win._use_structure_file(str(path))


def cu_file(tmp_path) -> Path:
    p = tmp_path / "cu.extxyz"
    write(p, bulk("Cu", "fcc", a=3.615, cubic=True))
    return p


def message(win: MainWindow) -> str:
    return win.preview.status.text() + "\n" + win.preview.detail.text()


def roundtrip(win: MainWindow, sk_root, tmp_path, cp2k_data=""):
    win.refresh_preview()
    spec = win.current_spec()
    win2 = window(sk_root, tmp_path / "b", cp2k_data)
    win2.apply_spec(spec)
    assert win2.current_spec().model_dump(exclude={"meta"}) == spec.model_dump(exclude={"meta"})
    return spec, win2


def test_codes_are_in_the_dropdown(app, quiet, sk_root, tmp_path):
    win = window(sk_root, tmp_path)
    items = [win.method.code.itemData(i) for i in range(win.method.code.count())]
    assert items[-4:] == ["cp2k", "lammps", "gromacs", "mlip"]
    assert all(not win.method.code.itemIcon(i).isNull() for i in range(win.method.code.count()))
    for code, panel in (("cp2k", win.method.cp2k), ("lammps", win.method.lammps), ("gromacs", win.method.gromacs)):
        use_code(win, code)
        assert win.method.stack.currentWidget() is panel


def test_cp2k_roundtrip_and_candidates(app, quiet, sk_root, tmp_path, cp2k_data):
    win = window(sk_root, tmp_path, cp2k_data)
    use_code(win, "cp2k")
    p = win.method.cp2k
    assert set(p.basis_widgets) == {"O", "H"}
    o_basis = [p.basis_widgets["O"].itemText(i) for i in range(p.basis_widgets["O"].count())]
    assert o_basis == ["", "DZVP-MOLOPT-SR-GTH", "SZV-MOLOPT-SR-GTH"]
    assert "GTH-PBE-q6" in p.potential_widgets["O"].lineEdit().placeholderText()
    p.xc.setCurrentText("PBE"); p.cutoff.setValue(400); p.rel_cutoff.setValue(60)
    p.poisson.setCurrentIndex(p.poisson.findData("MT")); p.box.setValue(10)
    p.basis_widgets["O"].setCurrentText("DZVP-MOLOPT-SR-GTH"); p.basis_widgets["H"].setCurrentText("DZVP-MOLOPT-SR-GTH")
    p.dispersion.setCurrentIndex(p.dispersion.findData("none"))
    p.extra.setPlainText("[FORCE_EVAL/DFT/SCF]\nSCF_GUESS ATOMIC")
    win.refresh_preview()
    assert win.btn_generate.isEnabled(), message(win)
    inp = win.preview.editors["cp2k.inp"].toPlainText()
    assert "BASIS_SET DZVP-MOLOPT-SR-GTH" in inp and "POTENTIAL GTH-PBE-q6" in inp and "SCF_GUESS ATOMIC" in inp
    spec, win2 = roundtrip(win, sk_root, tmp_path, cp2k_data)
    assert spec.method.code == "cp2k" and spec.method.basis == {"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"}
    assert spec.method.potential == {}
    assert win2.method.stack.currentWidget() is win2.method.cp2k


def test_lammps_hides_kpoints_and_roundtrips(app, quiet, sk_root, tmp_path):
    win = window(sk_root, tmp_path)
    use_file(win, cu_file(tmp_path))
    win.refresh_preview()
    win.set_mode(win.MODE_SETTINGS)
    assert win.kpoints.isVisibleTo(win)
    use_code(win, "lammps")
    assert not win.kpoints.isVisibleTo(win)
    eam = tmp_path / "Cu_u3.eam"; eam.write_text("fake eam\n", encoding="utf-8")
    p = win.method.lammps
    p.units.setCurrentIndex(p.units.findData("metal")); p.pair_style.setText("eam"); p.pair_coeff.setPlainText("* * Cu_u3.eam")
    p.potential_files.setPlainText(str(eam)); p.seed.setValue(4242)
    win.task.type.setCurrentIndex(2)
    win.refresh_preview()
    assert win.btn_generate.isEnabled(), message(win)
    assert "pair_style eam" in win.preview.editors["in.lammps"].toPlainText()
    spec, _ = roundtrip(win, sk_root, tmp_path)
    assert spec.kpoints is None and spec.method.potential_files == [str(eam)] and spec.method.seed == 4242
    use_code(win, "cp2k")
    assert win.kpoints.isVisibleTo(win)


def test_gromacs_structure_file_and_roundtrip(app, quiet, sk_root, tmp_path):
    win = window(sk_root, tmp_path)
    use_code(win, "gromacs")
    p = win.method.gromacs
    p.topology.setText(str(SPCE / "topol.top")); p.structure_file.setText(str(SPCE / "conf.gro"))
    p.load_structure.click()
    st = win.structure.structure()
    assert st is not None and len(st.atoms.symbols) == 648
    assert not win.kpoints.isVisibleTo(win)
    p.extra.setPlainText("nstcalcenergy = 100")
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and "grompp" in message(win)
    p.rcoulomb.setValue(0.85); p.rvdw.setValue(0.85)
    win.refresh_preview()
    assert win.btn_generate.isEnabled(), message(win)
    assert "nstcalcenergy" in win.preview.editors["grompp.mdp"].toPlainText()
    spec, _ = roundtrip(win, sk_root, tmp_path)
    assert spec.method.extra_mdp == {"nstcalcenergy": 100} and spec.structure.source == "file"


def test_required_fields_and_unsupported_tasks_show_the_reason(app, quiet, sk_root, tmp_path, cp2k_data):
    win = window(sk_root, tmp_path, cp2k_data)
    use_file(win, cu_file(tmp_path)); use_code(win, "lammps")
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and win.error_badge_action.isVisible()
    status = message(win)
    assert "pair_style" in status and "units" in status
    assert win.error_badge.text() == f"{len(win._errors)} 件の不足" and len(win._errors) >= 2
    win.task.type.setCurrentIndex(3)
    win.refresh_preview()
    assert "一点計算・構造最適化 (最小化)・分子動力学だけ" in message(win)
    use_code(win, "cp2k"); win.task.type.setCurrentIndex(4)
    win.refresh_preview()
    status = message(win)
    assert "バンド計算はありません" in status and "汎関数" in status


def test_labels_have_help_and_english(app, quiet, sk_root, tmp_path):
    win = window(sk_root, tmp_path)
    missing = []
    for panel in (win.method.cp2k, win.method.lammps, win.method.gromacs):
        form = panel.layout(); assert isinstance(form, QFormLayout)
        for r in range(form.rowCount()):
            item = form.itemAt(r, QFormLayout.ItemRole.LabelRole)
            text = item.widget().text() if item and item.widget() else ""
            if text and (help_for(text) is None or text not in _EN):
                missing.append(text)
    assert missing == []


def test_english_ui(app, quiet, sk_root, tmp_path):
    from PySide6.QtWidgets import QLabel
    from adit.gui.i18n import set_language, translate_widgets
    set_language("en")
    try:
        win = window(sk_root, tmp_path)
        translate_widgets(win)
        texts = {w.text() for w in win.method.findChildren(QLabel)}
        assert {"Functional", "Cutoff [Ry]", "Units", "Files to copy", "Topology (.top)", "Compressibility [1/bar]"} <= texts
        assert win.method.gromacs.load_structure.text() == "Also load into the structure group"
        assert win.method.lammps.units.itemText(0) == "(choose one)"
        use_file(win, cu_file(tmp_path)); use_code(win, "lammps"); win.refresh_preview()
        assert "choose the unit system" in message(win)
    finally:
        set_language("ja")


def test_settings_dialog_edits_cp2k_data(app, quiet, tmp_path, cp2k_data):
    from adit.gui.settings_dialog import SettingsDialog
    path = tmp_path / "cluster.toml"; save_config(default_config(), path)
    dlg = SettingsDialog(path)
    assert dlg.cp2k_data.text() == ""
    dlg.cp2k_data.setText(str(cp2k_data)); dlg._save()
    assert load_config(path).cp2k_data == str(cp2k_data)
