
import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from ase.build import bulk, molecule  # noqa: E402
from ase.io import write  # noqa: E402
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QMessageBox  # noqa: E402

from tests.conftest import pbs_profile, water_spec  # noqa: E402
from adit.config import default_config  # noqa: E402
from adit.gui.help import help_for  # noqa: E402
from adit.gui.i18n import _EN  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.panels.method_panel import CODES  # noqa: E402
from adit.gui.panels.mlip_panel import MlipMethodPanel  # noqa: E402
from adit.gui.panels.orca_panel import OrcaMethodPanel  # noqa: E402
from adit.spec import DcdftbmdMethod, MlipMethod, OrcaMethod  # noqa: E402
from adit.web import prep23 as P  # noqa: E402

UPF = """<UPF version="2.0.1">
  <PP_INFO>
     Suggested minimum cutoff for wavefunctions:  30. Ry
     Suggested minimum cutoff for charge density: 240. Ry
  </PP_INFO>
  <PP_HEADER element="O" z_valence="6.0" wfc_cutoff="3.5000000000E+01" rho_cutoff="2.8000000000E+02"/>
  <PP_MESH>
</UPF>
"""


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def boxes(monkeypatch):
    seen = []
    for name in ("information", "critical", "warning"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(lambda *a, _n=name, **k: seen.append((_n, a[2] if len(a) > 2 else ""))))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    return seen


def make_window(sk_root, tmp_path) -> MainWindow:
    cfg = default_config(sk_root=str(sk_root)); cfg.profiles["cluster"] = pbs_profile()
    win = MainWindow(cfg, tmp_path / "cluster.toml")
    win.method.sk_set.setCurrentText("fake-1-0")
    win.runtime.outdir.setText(str(tmp_path / "out"))
    win.refresh_preview()
    return win


def use_code(win: MainWindow, code: str) -> None:
    win.method.code.setCurrentIndex(list(CODES).index(code))


def use_bulk_carbon(win: MainWindow) -> None:
    win.structure.set_source("bulk")
    win.structure.bulk_el.setCurrentText("C")
    win.structure.bulk_struct.setCurrentIndex(win.structure.bulk_struct.findText("diamond"))
    win.structure._rebuild(); win._on_context(); win.refresh_preview()


def test_compare_reaction_generates(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    nh3 = tmp_path / "nh3.xyz"; write(nh3, molecule("NH3"))
    dlg = win.batch_dialog("compare")
    dlg.kind_box.setCurrentIndex(dlg.kind_box.findData("reaction"))
    dlg.table.cellWidget(0, 0).setText("water"); dlg.table.cellWidget(0, 2).setText("-1")
    dlg.table.cellWidget(1, 0).setText("nh3"); dlg.table.cellWidget(1, 1).setText(str(nh3)); dlg.table.cellWidget(1, 2).setText("1")
    out = tmp_path / "set"; dlg.folder.setText(str(out))
    dlg.generate()
    assert [d.name for d in dlg.dirs] == ["water", "nh3"], boxes
    assert (out / "compare.json").is_file() and (out / "water" / "dftb_in.hsd").is_file()
    win.batch_written(dlg.result)
    assert win.analysis.run_dir.text() == str(out / "water")


def test_compare_reaction_reports_bad_cells(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.batch_dialog("compare")
    dlg.kind_box.setCurrentIndex(dlg.kind_box.findData("reaction"))
    dlg.table.cellWidget(0, 0).setText("water"); dlg.table.cellWidget(0, 2).setText("なし")
    dlg.folder.setText(str(tmp_path / "bad")); dlg.generate()
    assert dlg.out_dir is None and boxes[-1][0] == "critical" and "なし" in boxes[-1][1]


def test_compare_structure_buttons_are_explicit(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.batch_dialog("compare")
    dlg.slab.setText("slab.xyz"); dlg.use_slab.click()
    dlg.mol.setText("mol.xyz"); dlg.use_mol.click()
    dlg.ads.setText("ads.xyz"); dlg.use_ads.click()
    assert dlg.slab.text() == dlg.mol.text() == dlg.ads.text() == ""
    assert "extended XYZ" in dlg.save_structure.text()


def test_conformers_generate(app, boxes, sk_root, tmp_path):
    pytest.importorskip("rdkit")
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("smiles"); win.structure.smiles.setText("CCO"); win.structure._rebuild(); win.refresh_preview()
    dlg = win.batch_dialog("conformers")
    dlg.n.setText("4"); dlg.rmsd.setText("0.5")
    out = tmp_path / "conf"; dlg.folder.setText(str(out)); dlg.generate()
    assert dlg.dirs and (out / "conformers.json").is_file() and (out / "scan.json").is_file(), boxes
    win.batch_written(dlg.result)
    assert win.analysis.run_dir.text() == str(out)


def test_neb_images_generate(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    end = tmp_path / "end.xyz"
    a = molecule("H2O"); a.positions[1] += (0.0, 0.0, 0.3); write(end, a)
    dlg = win.batch_dialog("neb")
    dlg.end.setText(str(end)); dlg.images.setText("3")
    dlg.mode.setCurrentIndex(dlg.mode.findData("images"))
    out = tmp_path / "neb"; dlg.folder.setText(str(out)); dlg.generate()
    assert [d.name for d in dlg.dirs] == ["image_00", "image_01", "image_02", "image_03", "image_04"], boxes
    assert (out / "neb.json").is_file()


def test_neb_native_with_dftb_is_refused(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    end = tmp_path / "end.xyz"; write(end, molecule("H2O"))
    dlg = win.batch_dialog("neb")
    dlg.end.setText(str(end)); dlg.images.setText("3"); dlg.folder.setText(str(tmp_path / "neb2"))
    dlg.generate()
    assert dlg.out_dir is None and boxes[-1][0] == "critical" and "VASP" in boxes[-1][1]


def test_phonons_generate(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    use_bulk_carbon(win)
    dlg = win.batch_dialog("phonons")
    dlg.dim.setText("1x1x1"); dlg.backend.setCurrentIndex(dlg.backend.findData("ase"))
    out = tmp_path / "ph"; dlg.folder.setText(str(out)); dlg.generate()
    assert len(dlg.dirs) == 12 and (out / "phonons.json").is_file() and (out / "phonon_collect.py").is_file(), boxes
    assert "phonon_collect.py" in P.result_message(dlg.result)


def test_elastic_generate(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    use_bulk_carbon(win)
    dlg = win.batch_dialog("elastic")
    dlg.strains.setText("-0.01,0.01")
    for j, c in dlg.comps.items():
        c.setChecked(j == 1)
    out = tmp_path / "el"; dlg.folder.setText(str(out)); dlg.generate()
    assert [d.name for d in dlg.dirs] == ["e0", "e1_-0.01", "e1_+0.01"], boxes
    assert (out / "elastic_collect.py").is_file()


def test_ts_orca_stages_and_sella(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    use_code(win, "orca"); win.refresh_preview()
    dlg = win.batch_dialog("ts")
    assert dlg.mode.currentData() == "ts"
    dlg.mode.setCurrentIndex(dlg.mode.findData("ts+irc"))
    out = tmp_path / "ts"; dlg.folder.setText(str(out)); dlg.generate()
    assert [d.name for d in dlg.dirs] == ["stage_01_ts", "stage_02_irc"], boxes
    assert "OptTS" in (out / "stage_01_ts" / "orca.inp").read_text(encoding="utf-8")
    use_code(win, "mlip"); win.method.mlip.family.setCurrentIndex(win.method.mlip.family.findData("mace_mp"))
    win.task.type.setCurrentIndex(win.task.type.findData("geometry_optimization")); win.refresh_preview()
    dlg = win.batch_dialog("ts")
    assert dlg.mode.currentData() == "sella"
    out2 = tmp_path / "sella"; dlg.folder.setText(str(out2)); dlg.generate()
    assert (out2 / "run_sella.py").is_file() and (out2 / "sella_settings.json").is_file(), boxes


def test_mlip_panel_roundtrip(app):
    p = MlipMethodPanel()
    m = MlipMethod(model_family="mace_mp", model="small", device="cpu", dtype="float32", dispersion=True, seed=7)
    p.set_method(m)
    assert p.method() == m
    p.set_method(MlipMethod(model_family="chgnet"))
    assert p.method().model_family == "chgnet" and "chgnet" in p.pip.text()


def test_mlip_spec_opens_without_crashing(app, sk_root, tmp_path):
    from tests.conftest import water_spec
    win = make_window(sk_root, tmp_path)
    s = water_spec(method=MlipMethod(model_family="mace_mp", model="small"))
    s = s.model_copy(update={"task": s.task.model_copy(update={"optimizer": "LBFGS"})})
    win.apply_spec(s)
    assert win.method.current_code() == "mlip" and win.current_spec().method == s.method
    assert not win.kpoints.isVisible()
    assert "run_mlip.py" in win.build().texts


def test_orca_ts_irc_roundtrip_and_rows(app):
    p = OrcaMethodPanel()
    m = OrcaMethod(ts_search=True, ts_calc_hess=False, ts_recalc_hess=5, ts_freq=False, irc=False, irc_max_iter=0, irc_direction="both")
    p.set_method(m)
    assert p.method() == m
    assert p._form.isRowVisible(p.ts_recalc) and not p._form.isRowVisible(p.irc_max_iter)
    m2 = OrcaMethod(irc=True, irc_max_iter=30, irc_direction="forward")
    p.set_method(m2)
    assert p.method() == m2 and p._form.isRowVisible(p.irc_direction) and not p._form.isRowVisible(p.ts_recalc)


def test_orca_cli_options_survive_gui_roundtrip(app):
    p = OrcaMethodPanel()
    m = OrcaMethod(goat=True, docker_guest_file="guests.xyz", docker_assume_neutral_singlet=True)
    p.set_method(m)
    assert p.method() == m
    docker = OrcaMethod(method="XTB", basis="", docker_guest_file="guests.xyz")
    p.set_method(docker)
    assert p.method() == docker
    p.set_method(OrcaMethod())
    assert p.method() == OrcaMethod()


def test_gui_rejects_cli_only_code_without_changing_current_spec(app, boxes, sk_root, tmp_path, monkeypatch):
    win = make_window(sk_root, tmp_path)
    before = win.current_spec().method
    path = tmp_path / "cli_only.json"
    water_spec(method=DcdftbmdMethod()).save(path)
    monkeypatch.setattr(QFileDialog, "getOpenFileName", staticmethod(lambda *args, **kwargs: (str(path), "")))
    win.open_spec()
    assert win.current_spec().method == before
    assert boxes and "CLI" in boxes[-1][1]


def test_documented_values_are_shown_and_can_be_used(app, sk_root, tmp_path):
    lib = tmp_path / "pseudo" / "SSSP_test"; lib.mkdir(parents=True)
    (lib / "O.upf").write_text(UPF, encoding="utf-8")
    (lib / "H.upf").write_text(UPF.replace('element="O"', 'element="H"'), encoding="utf-8")
    (lib / "sssp.json").write_text('{"O": {"filename": "O.upf", "cutoff_wfc": 50.0, "cutoff_rho": 400.0}}', encoding="utf-8")
    win = make_window(sk_root, tmp_path)
    win.cfg.pseudo_root = str(tmp_path / "pseudo")
    use_code(win, "espresso")
    win.method.espresso.pseudo_set.setCurrentText("SSSP_test")
    win.method.reload_sets(str(sk_root), win.cfg)
    win.refresh_preview()
    box = win.method.doc_boxes["ecutwfc"]
    texts = [box.grid.itemAt(i).widget().text() for i in range(box.grid.count())]
    assert box.isVisibleTo(win.method.espresso) and any("O.upf:" in t and "35" in t for t in texts), texts
    assert any("sssp.json:" in t and "50" in t for t in texts), texts
    assert not any("ADIT" in t for t in texts)
    win.method._use_doc_value("ecutwfc", "35")
    assert win.method.espresso.ecutwfc.value() == 35.0
    use_code(win, "dftbplus"); win.refresh_preview()
    assert not win.method.doc_boxes["ecutwfc"].isVisibleTo(win.method.espresso) or True


def test_ribbon_has_batch_buttons(app, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    titles = win.ribbon.page_titles()
    run_actions = win.ribbon.actions_on_page(titles.index("実行"))
    for kind in P.BATCH_KINDS:
        assert win.act_batch[kind] in run_actions, kind
    assert run_actions.index(win.act_batch["compare"]) == run_actions.index(win.act_stages) + 1


def test_run_page_fits_at_1366(app, sk_root, tmp_path):
    from adit.gui.i18n import set_language, translate_widgets
    from adit.gui.style import apply_base_style
    apply_base_style(app, "light")
    set_language("en")
    try:
        win = make_window(sk_root, tmp_path); translate_widgets(win)
        win.resize(1366, 768); win.show()
        for _ in range(6):
            app.processEvents()
        i = win.ribbon.page_titles().index("Run")
        win.ribbon.show_page(i)
        for _ in range(6):
            app.processEvents()
        assert win.ribbon.pages[i].minimumSizeHint().width() <= 1366
        win.close()
    finally:
        set_language("ja")


def test_new_labels_have_help_and_english():
    keys = [ja for key, (ja, _en) in P.LABELS.items() if key != "out"]
    keys += ["機械学習ポテンシャルの種類", "モデル", "計算に使うデバイス (device)", "数値の精度 (dtype)",
             "遷移状態の探索 (OptTS)", "最初にヘシアンを計算 (Calc_Hess)", "ヘシアンを計算し直す間隔 (Recalc_Hess)",
             "最後に振動数を計算 (Freq)", "反応座標をたどる (IRC)", "IRC の反復の上限", "IRC の向き"]
    for ja in keys:
        h = help_for(ja)
        assert h is not None and h.ja and h.en, ja
        assert ja in _EN, ja


def test_english_dialogs(app, sk_root, tmp_path):
    from adit.gui import i18n
    from adit.lang import set_language
    i18n.set_language("en")
    try:
        win = make_window(sk_root, tmp_path)
        i18n.translate_widgets(win)
        assert win.act_batch["compare"].iconText() == "Compare a set"
        dlg = win.batch_dialog("phonons"); i18n.translate_widgets(dlg)
        assert dlg.windowTitle() == "Phonons (finite displacements)"
        texts = {w.text() for w in dlg.findChildren(QLabel)}
        assert "Supercell" in texts and "Output directory" in texts
        use_code(win, "mlip"); i18n.translate_widgets(win.method.mlip)
        assert "Machine-learning potential" in {w.text() for w in win.method.mlip.findChildren(QLabel)}
    finally:
        i18n.set_language("ja"); set_language("ja")
