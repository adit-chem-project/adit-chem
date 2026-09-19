
import os
import shutil
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEventLoop, QTimer  # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from adit.config import default_config  # noqa: E402
from tests.conftest import pbs_profile  # noqa: E402
from adit.project import load_project  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from tests.conftest import REAL_SK_ROOT  # noqa: E402

DFTB_EXE = os.environ.get("ADIT_DFTB_EXE") or shutil.which("dftb+") or ""
HAVE_REAL = (REAL_SK_ROOT / "mio-1-1").is_dir() and Path(DFTB_EXE).is_file()



@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _delete_windows_after_test():
    yield
    from PySide6.QtCore import QCoreApplication, QEvent
    app = QApplication.instance()
    if app is None:
        return
    for w in app.topLevelWidgets():
        if w.parentWidget() is None:
            w.close(); w.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


@pytest.fixture
def quiet(monkeypatch):
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))


def make_window(sk_root, tmp_path) -> MainWindow:
    cfg = default_config(sk_root=str(sk_root))
    cfg.profiles["cluster"] = pbs_profile()
    win = MainWindow(cfg, tmp_path / "cluster.toml")
    win.method.sk_set.setCurrentText(sorted(cfg.profiles and win.method.sets)[0])
    win.runtime.outdir.setText(str(tmp_path / "out"))
    win.refresh_preview()
    return win


def test_defaults_generate_without_typing(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0")
    win.refresh_preview()
    assert win.btn_generate.isEnabled(), win.preview.status.text()
    assert win.structure.preset.currentData() == "H2O" and win.structure.preset.currentText() == "H₂O"
    assert 'O = "p"' in win.preview.editors["dftb_in.hsd"].toPlainText()
    win.generate()
    out = tmp_path / "out"
    assert (out / "dftb_in.hsd").is_file() and (out / "skf" / "LICENSE").is_file()
    spec = load_project(out)
    assert spec.structure.source == "preset" and spec.method.sk_set == "fake-1-0"


def test_errors_disable_generate(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0")
    win.structure.set_source("file"); win.structure.file.setText(str(tmp_path / "absent.xyz")); win.structure._rebuild()
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and "absent.xyz" in win.preview.status.text()


def test_spec_roundtrip_through_panels(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0")
    win.task.type.setCurrentIndex(1); win.task.max_steps.setValue(37)
    win.runtime.profile.setCurrentText("cluster"); win.runtime.ncpus.setValue(16); win.runtime.omp.setValue(16)
    win.refresh_preview()
    spec = win.current_spec()
    assert spec.task.type == "geometry_optimization" and spec.task.max_steps == 37
    assert "ncpus=16" in win.preview.editors["submit.sh"].toPlainText()
    win2 = make_window(sk_root, tmp_path / "b")
    win2.apply_spec(spec)
    assert win2.current_spec().model_dump(exclude={"meta"}) == spec.model_dump(exclude={"meta"})






def test_theme_falls_back_without_qdarktheme(app, monkeypatch):
    import builtins
    from adit.gui.style import LIGHT, apply_theme, qss
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "qdarktheme":
            raise ImportError("simulated: fork is gone")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    name = apply_theme(app, "light")
    assert name.startswith("base") and "simulated" in name
    assert app.styleSheet() == qss(LIGHT)


def test_theme_applies_when_available(app):
    pytest.importorskip("qdarktheme")
    from adit.gui.style import apply_theme
    assert apply_theme(app, "light") == "qdarktheme:light"
    assert "QGroupBox" in app.styleSheet()


def test_gui_vasp_bulk_flow(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("bulk"); win.structure.bulk_el.setCurrentText("Si"); win.structure._rebuild()
    assert "周期系" in win.structure.info.text()
    assert win.kpoints.isEnabled()
    win.method.code.setCurrentIndex(1)  # VASP
    win.kpoints.mode.setCurrentIndex(1)
    for w in win.kpoints.mesh:
        w.setValue(4)
    win.method.vasp.encut.setValue(240)
    win.method.vasp.extra.setPlainText("NCORE = 4\nLASPH = .TRUE.")
    win.runtime.profile.setCurrentText("cluster")
    win.refresh_preview()
    assert win.btn_generate.isEnabled(), win.preview.status.text()
    assert "INCAR" in win.preview.editors and "dftb_in.hsd" not in win.preview.editors
    inc = win.preview.editors["INCAR"].toPlainText()
    assert "ENCUT = 240" in inc and "NCORE = 4" in inc and "LASPH = .TRUE." in inc
    assert "Monkhorst" in win.preview.editors["KPOINTS"].toPlainText()
    assert "Si  Si  (未確認)" in win.preview.editors["potcar.spec"].toPlainText()
    win.generate()
    out = tmp_path / "out"
    assert (out / "INCAR").is_file() and (out / "make_potcar.sh").is_file() and not (out / "POTCAR").exists()
    assert "transfer_and_submit.sh" in win.run_hint.text()   # the hint for a cluster profile
    spec = load_project(out)
    assert spec.method.code == "vasp" and spec.kpoints.mesh == (4, 4, 4) and spec.structure.source == "bulk"
    win2 = make_window(sk_root, tmp_path / "b")
    win2.apply_spec(spec)
    assert win2.current_spec().method == spec.method and win2.current_spec().kpoints == spec.kpoints


def test_gui_molecule_in_box_and_fixed_atoms(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.code.setCurrentIndex(1)
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and "周期セル" in win.preview.status.text() + win.preview.detail.text()
    win.structure.box.setChecked(True); win.structure.fixed.setText("1"); win.structure._rebuild()
    win.refresh_preview()
    assert win.btn_generate.isEnabled(), win.preview.status.text()
    pos = win.preview.editors["POSCAR"].toPlainText()
    assert "Selective dynamics" in pos and "F   F   F" in pos
    win.structure.fixed.setText("1-9"); win.structure._rebuild()
    assert "範囲外" in win.structure.error()


def test_gui_extra_incar_error_is_shown(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("bulk"); win.structure._rebuild()
    win.method.code.setCurrentIndex(1)
    win.method.vasp.extra.setPlainText("this is not a key value line")
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and "KEY = value" in win.preview.status.text()




def test_english_ui(app, quiet, sk_root, tmp_path):
    from adit.gui.i18n import set_language, translate_widgets, tr
    set_language("en")
    try:
        win = make_window(sk_root, tmp_path)
        translate_widgets(win)
        assert win.structure.title() == "Structure" and win.btn_generate.text() == "Generate"
        assert win.structure.source.itemText(win.structure.source.findData("bulk")) == "Bulk" and win.kpoints.mode.itemText(0) == "Γ only"
        assert tr("生成できます") == "Ready to generate" and tr("未知の文") == "未知の文"
        assert win.ribbon.page_titles() == ["File", "Insert", "View", "Run", "Settings", "Help"]
        assert win.act_open.iconText() == "Open" and win.act_open.text() == "Open calculation settings (spec.json)…"
        assert win.act_back.toolTip() == "Undo (previous settings)"
        assert [b.text() for b in win.ribbon.pages[0].groups[0].buttons] == ["Open", "Save"]
        assert win.ribbon.pages[0].groups[0].caption.text() == "Calculation settings"
        win.refresh_preview()
        assert win.preview.status.text() == "Ready to generate"
    finally:
        set_language("ja")


def test_analysis_tab(app, quiet, sk_root, tmp_path):
    import shutil
    from pathlib import Path
    win = make_window(sk_root, tmp_path)
    d = tmp_path / "run"
    shutil.copytree(Path(__file__).resolve().parent.parent / "examples" / "dftb_md_water_generated", d)
    win.analysis.set_run_dir(d, ["O", "H"])
    win.analysis.cb_rdf.setChecked(True); win.analysis.cb_msd.setChecked(True)
    res = win.analysis.run()
    assert res is not None and "temperature" in res.figures and "rdf" in res.figures
    assert "温度" in win.analysis.summary.toPlainText()
    assert win.analysis.figs_lay.count() >= 6


def test_structure_source_is_one_dropdown_with_icons(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    src = win.structure.source
    assert src.count() == 9 and src.currentData() == "preset"
    assert [src.itemData(i) for i in range(6, 9)] == ["2d", "cluster", "polymer"]
    assert all(not src.itemIcon(i).isNull() for i in range(src.count())), "すべての項目に絵記号"
    form = win.structure._form
    assert form.isRowVisible(win.structure.preset) and not form.isRowVisible(win.structure.mixture)
    win.structure.set_source("mixture")
    assert win.structure.current_source() == "mixture"
    assert form.isRowVisible(win.structure.mixture) and not form.isRowVisible(win.structure.preset)


def test_toolbar_has_no_menu_duplicates(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    shown = [a for a in win.action_bar.actions() if a.text()]
    assert shown == [], f"操作の列に文字のボタンが残っている: {[a.text() for a in shown]}"
    assert win.act_back not in win.action_bar.actions() and win.act_forward not in win.action_bar.actions()
    assert [b.defaultAction() for b in win.quick_access.buttons] == [win.act_back, win.act_forward]
    assert win.ribbon.page_titles() == ["ファイル", "挿入", "表示", "実行", "設定", "ヘルプ"]
    ribbon_texts = {a.text() for a in win.ribbon.all_actions()}
    for text in ("計算設定 (spec.json) を開く…", "環境設定を再読み込み", "環境設定…", "Draw", "1 つの条件を変えて一括生成…", "生成"):
        assert text in ribbon_texts, text
    assert win.act_scan in win.ribbon.actions_on_page(3) and win.act_generate in win.ribbon.actions_on_page(3)
    assert all(not a.icon().isNull() and a.iconText() for a in win.ribbon.all_actions())
    top = win.menuWidget()
    assert not win.menuBar().isVisibleTo(win) and win.menuWidget() is top is win.top_area


def test_ribbon_collapses_on_second_click(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.show(); app.processEvents()
    rb = win.ribbon
    h_open = rb.height()
    rb.tabs.tabBarClicked.emit(rb.tabs.currentIndex())
    for _ in range(5):
        app.processEvents()
    assert rb.is_collapsed() and not rb.stack.isVisible() and rb.tabs.isVisible()
    assert rb.height() < h_open - 40, (rb.height(), h_open)
    rb.tabs.tabBarClicked.emit(3); rb.tabs.setCurrentIndex(3)
    assert not rb.is_collapsed() and rb.stack.currentIndex() == 3
    rb.toggle.click(); assert rb.is_collapsed()
    rb.toggle.click(); assert not rb.is_collapsed()
    win.ribbon.pages[2].groups[0].buttons[2].click()
    assert win.mode() == win.MODE_ANALYSIS and win.main_stack.currentIndex() == win.MODE_ANALYSIS
    win.structure.set_source("file"); win.structure._rebuild(); win.refresh_preview()
    assert not win.btn_generate.isEnabled() and not win.act_generate.isEnabled()
    win.close()


def test_long_dropdowns_are_capped_and_scroll(app, quiet, sk_root, tmp_path):
    from PySide6.QtWidgets import QComboBox, QStyle, QStyleOptionComboBox
    from adit.gui.style import apply_base_style
    from adit.gui.widgets import POPUP_ROWS
    from adit.mixture import Component
    apply_base_style(app, "light")
    win = make_window(sk_root, tmp_path); win.show(); app.processEvents()
    for combo in (win.structure.preset, win.structure.bulk_el, win.structure.surf_el, win.method.sk_set, win.task.thermostat):
        assert combo.maxVisibleItems() == POPUP_ROWS, combo
    win.structure.set_source("mixture"); win.structure.mixture.add_component(Component(ref="CH4", count=2)); app.processEvents()
    row = win.structure.mixture.ref_cell(win.structure.mixture.table.rowCount() - 1).preset
    assert row.count() > 100 and row.maxVisibleItems() == POPUP_ROWS
    late = QComboBox(win); late.addItems([str(i) for i in range(50)]); late.show(); app.processEvents()
    assert late.maxVisibleItems() == POPUP_ROWS
    opt = QStyleOptionComboBox(); opt.initFrom(row)
    assert not row.style().styleHint(QStyle.StyleHint.SH_ComboBox_Popup, opt, row)
    row.showPopup(); app.processEvents()
    popup = row.view().window()
    assert popup.height() < 20 * row.view().sizeHintForRow(0) and row.view().verticalScrollBar().maximum() > 0
    row.hidePopup(); win.close()


def test_mixture_reference_column_has_real_inputs(app, quiet, sk_root, tmp_path):
    from adit.gui.panels.mixture_editor import MixtureEditor
    from adit.mixture import Component, MixtureSpec
    ed = MixtureEditor()
    cell = ed.ref_cell(0)
    assert cell.kind() == "preset" and cell.preset.currentText() == "H₂O" and cell.text() == "H2O"
    kind = ed.table.cellWidget(0, 0)
    kind.setCurrentIndex(kind.findData("smiles"))
    assert cell.kind() == "smiles" and cell.stack.currentWidget().isAncestorOf(cell.btn_draw) and cell.text() == ""
    ed.set_ref(0, "CCO"); assert ed.spec().components[0].kind == "smiles" and ed.spec().components[0].ref == "CCO"
    kind.setCurrentIndex(kind.findData("file"))
    assert cell.stack.currentWidget().isAncestorOf(cell.btn_browse)
    kind.setCurrentIndex(kind.findData("preset")); assert cell.text() == "H2O"
    m = MixtureSpec(components=[Component(ref="CH4", count=3, label="CH4"), Component(kind="smiles", ref="[Na+]", charge=1),
                                Component(kind="file", ref="/tmp/ion.xyz", count=2)], box_a=12.0, seed=5)
    ed.set_spec(m)
    assert ed.spec() == m
    assert ed.ref_cell(0).preset.currentText() == "CH₄"
    seen = []
    ed.draw_requested.connect(seen.append)
    ed.ref_cell(1).btn_draw.click()
    assert seen == [1] and ed.row_smiles(1) == "[Na+]"
    ed.set_row_smiles(-1, "[Cl-]")
    assert ed.table.rowCount() == 4 and ed.spec().components[3].ref == "[Cl-]"


def test_draw_buttons_have_no_ellipsis(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    assert win.structure.draw_button.text() == win.structure.smiles_draw.text() == win.structure.mixture.btn_draw.text() == "Draw"
    assert win.act_draw.text() == "Draw" and win.structure.mixture.ref_cell(0).btn_draw.text() == "Draw"


def test_combos_have_icons(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    for combo in (win.method.code, win.task.type, win.kpoints.mode, win.runtime.profile):
        assert combo.count() > 0
        assert all(not combo.itemIcon(i).isNull() for i in range(combo.count())), combo.objectName()
    assert all(not win.right_tabs.tabIcon(i).isNull() for i in range(win.right_tabs.count()))


def test_titlebar_buttons_light_up_on_hover(app, quiet, sk_root, tmp_path, monkeypatch):
    from PySide6.QtCore import QEvent, QPoint, QPointF, Qt
    from PySide6.QtGui import QEnterEvent
    from adit.gui.titlebar import BUTTON_COLORS

    monkeypatch.setenv("ADIT_FRAME", "custom")
    win = make_window(sk_root, tmp_path); win.setGeometry(100, 100, 1400, 900); win.show(); app.processEvents()
    assert win.windowFlags() & Qt.WindowType.FramelessWindowHint
    assert win.menuWidget() is win.top_area and win.top_area.isAncestorOf(win.titlebar) and win.top_area.isAncestorOf(win.ribbon)
    assert win.titlebar.leading is win.quick_access and win.titlebar.isAncestorOf(win.quick_access)
    assert "ファイル" in win.ribbon.page_titles()
    bar = win.titlebar
    for b in (bar.btn_close, bar.btn_min, bar.btn_max):
        idle = b.fill_color()
        QApplication.sendEvent(b, QEnterEvent(QPointF(5, 5), QPointF(5, 5), QPointF(5, 5)))
        assert b.fill_color() == BUTTON_COLORS[b.kind][0] != idle
        QApplication.sendEvent(b, QEvent(QEvent.Type.Leave))
        assert b.fill_color() == idle
    fg = win.frameGeometry()
    assert win._edge_resizer.edges_at(QPoint(fg.x() + 2, fg.y() + 300)) == Qt.Edge.LeftEdge
    assert win._edge_resizer.edges_at(QPoint(fg.right() - 1, fg.bottom() - 1)) == Qt.Edge.RightEdge | Qt.Edge.BottomEdge
    assert win._edge_resizer.edges_at(fg.center()) == Qt.Edge(0)
    bar.btn_close.click()
    assert not win.isVisible()


def test_switching_source_to_empty_file_drops_old_structure(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    assert win.btn_generate.isEnabled()
    for src in ("file", "smiles"):
        win.structure.set_source(src); win.structure._rebuild(); win.refresh_preview()
        assert win.structure.structure() is None, src
        assert not win.btn_generate.isEnabled()
        assert win.gen_hint_action.isVisible() and "生成できません" in win.gen_hint.text()


def test_generate_block_reason_is_shown_without_internal_names(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.reload_sets("")
    win.refresh_preview()
    assert not win.btn_generate.isEnabled()
    assert win.method.dftb.sk_browse.isVisibleTo(win.method.dftb) and "dftb.org" in win.method.dftb.sk_info.text()
    assert "Slater-Koster パラメータ" in win.gen_hint.toolTip() and "method.sk_set" not in win.gen_hint.toolTip()
    win.gen_hint.click()
    assert win.mode() == win.MODE_SETTINGS and win.method.sk_set.property("adit_error")


def test_espresso_rows_leave_no_stale_widgets(app, quiet, sk_root, tmp_path):
    from PySide6.QtWidgets import QWidget
    win = make_window(sk_root, tmp_path); win.show(); app.processEvents()
    win.structure.set_source("bulk"); win.structure._rebuild()
    win.method.code.setCurrentIndex(win.method.code.findData("espresso")); win.refresh_preview()
    stray = [c for c in win.method.espresso.findChildren(QWidget) if c.parentWidget() is win.method.espresso and c.isVisible() and c.pos().isNull()]
    assert stray == []


def test_bulk_row_does_not_push_the_method_column_off_screen(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("bulk"); win.structure._rebuild()
    assert win.structure.minimumSizeHint().width() < 600


def test_mixture_does_not_widen_the_structure_group(app, quiet, sk_root, tmp_path):
    from adit.mixture import Component
    win = make_window(sk_root, tmp_path)
    preset_w = win.structure.minimumSizeHint().width()
    win.structure.set_source("mixture")
    win.structure.mixture.add_component(Component(kind="smiles", ref="[Na+]", charge=1))
    win.structure.mixture.add_component(Component(kind="file", ref=""))
    limit = 620 if os.name == "nt" else 480
    width = win.structure.minimumSizeHint().width()
    assert width <= limit, (width, preset_w, limit)


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_left_columns_fit_at_1366(app, quiet, sk_root, tmp_path, lang):
    from adit.gui.i18n import set_language, translate_widgets
    from adit.gui.style import apply_base_style
    apply_base_style(app, "light")
    set_language(lang)
    try:
        win = make_window(sk_root, tmp_path); translate_widgets(win)
        win.set_mode(win.MODE_SETTINGS)
        sc = win._left_scroll

        def settle():
            for _ in range(6):
                app.processEvents()

        def check(tag):
            settle()
            need, have = sc.widget().minimumSizeHint().width(), sc.viewport().width()
            assert need <= have and not sc.horizontalScrollBar().isVisible(), (lang, tag, need, have, win.left_columns())

        win.resize(1366, 768); win.show(); settle()
        for src in ("preset", "bulk", "surface", "mixture", "smiles", "file"):
            win.structure.set_source(src); win.structure._rebuild(); win._on_context(); check(src)
        win.set_mode(win.MODE_STRUCTURE); settle()
        src_combo = win.structure.source
        assert src_combo.width() >= src_combo.sizeHint().width() - 2, (src_combo.width(), src_combo.sizeHint().width())
        win.set_mode(win.MODE_SETTINGS); settle()
        win.structure.set_source("bulk"); win.structure._rebuild(); win._on_context()
        for code in ("espresso", "vasp", "dftbplus"):
            win.method.code.setCurrentIndex(win.method.code.findData(code)); check(code)
        win.task.type.setCurrentIndex(win.task.type.findData("molecular_dynamics")); check("md")
        for width in (1600, 1800):
            win.resize(width, 1000); check(width)
        win.structure.set_source("preset"); win.structure._rebuild(); win._on_context(); check("preset 1800")
        assert win.left_columns() == 2
        win.close()
    finally:
        set_language("ja")


def test_mixture_hides_put_in_cell_row(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    assert win.structure._box_row.isVisibleTo(win.structure)
    win.structure.set_source("mixture")
    assert not win.structure._box_row.isVisibleTo(win.structure)


def test_small_periodic_cell_is_shown_repeated(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.set_source("bulk"); win.structure._rebuild(); win._on_context()
    assert win.structure_view.repeat.currentText() == "2×2×2" and "並べたもの" in win.structure_view.info.text()
    win.structure.set_source("preset"); win.structure._rebuild(); win._on_context()
    assert win.structure_view.repeat.currentText() == "1×1×1"


def test_run_progress_from_output_log(app, tmp_path):
    from adit.gui.panels.run_panel import RunPanel, last_step, progress_text
    log = "***  Geometry step: 0\n...\n***  Geometry step: 41\n"
    assert last_step(log) == 41 and last_step("no steps") is None
    t = progress_text(41, 1000, 84.0)
    assert "42 / 1000" in t and "1 ステップ 2 秒" in t and "残り約 32 分" in t
    assert "/" not in progress_text(3, None, 1.0).split("、")[0]
    (tmp_path / "output.log").write_text(log, encoding="utf-8")
    p = RunPanel(); p.start(tmp_path, 100)
    assert p.progress.value() == 42 and "42 / 100" in p.progress_label.text()
    p.finish(0, "done"); assert p.progress.isHidden()


def test_stale_block_reason_disappears_after_fix(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0")
    win.structure.set_source("smiles"); win.structure._rebuild(); win.refresh_preview()
    assert win.gen_hint_action.isVisible()
    win.structure.set_source("preset"); win.structure._rebuild(); win.refresh_preview()
    assert win.btn_generate.isEnabled() and not win.gen_hint_action.isVisible()


def test_empty_output_dir_blocks_generation(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0")
    win.runtime.outdir.setText(""); win.refresh_preview()
    assert not win.btn_generate.isEnabled() and "出力ディレクトリ" in win.gen_hint.toolTip()


def test_smiles_formal_charge_fills_charge_column(app):
    pytest.importorskip("rdkit")
    from adit.gui.panels.mixture_editor import MixtureEditor, smiles_formal_charge
    from adit.mixture import Component
    assert smiles_formal_charge("[Na+]") == 1 and smiles_formal_charge("[Cl-]") == -1 and smiles_formal_charge("CCO") == 0
    ed = MixtureEditor(); ed.add_component(Component(kind="smiles", ref="C", count=1))
    r = ed.table.rowCount() - 1
    ed.set_ref(r, "[Na+]")
    assert ed.table.cellWidget(r, 3).value() == 1
    cell = ed.ref_cell(r); cell.smiles.setText("[Cl-]"); cell.smiles.editingFinished.emit()
    assert ed.table.cellWidget(r, 3).value() == -1


def test_run_panel_finish_replaces_remaining_time(app, tmp_path):
    from adit.gui.panels.run_panel import RunPanel
    (tmp_path / "output.log").write_text("***  Geometry step: 0\n***  Geometry step: 6\n", encoding="utf-8")
    p = RunPanel(); p.start(tmp_path, 200); p.finish(0, "ok")
    assert "7 ステップで終了しました" in p.progress_label.text() and "残り" not in p.progress_label.text()


def test_unknown_preset_suggests_smiles():
    from adit.structure import StructureError, from_preset
    with pytest.raises(StructureError) as ex:
        from_preset("ethanol")
    assert "SMILES" in str(ex.value) and "CCO" in str(ex.value)


def test_preset_filter_accepts_english_and_japanese_names(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.structure.preset_search.setText("ethanol")
    assert win.structure.preset.currentData() == "CH3CH2OH"
    assert win.structure.preset.count() == 1
    win.structure.preset_search.setText("エタノール")
    assert win.structure.preset.currentData() == "CH3CH2OH"
    win.structure.preset_search.clear()
    assert win.structure.preset.count() == len(win.structure._preset_names)


def test_menu_checks_follow_reloaded_settings(app, quiet, sk_root, tmp_path):
    from adit.config import save_config
    win = make_window(sk_root, tmp_path)
    cfg = win.cfg.model_copy(deep=True); cfg.theme = "dark"; cfg.enable_run = False
    save_config(cfg, tmp_path / "cluster.toml")
    win.reload_config()
    assert [c for c, a in win._theme_actions if a.isChecked()] == ["dark"]


def test_native_frame_keeps_os_titlebar(app, quiet, sk_root, tmp_path, monkeypatch):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMenuBar

    monkeypatch.setenv("ADIT_FRAME", "native")
    win = make_window(sk_root, tmp_path)
    assert not win.windowFlags() & Qt.WindowType.FramelessWindowHint
    assert win.titlebar is None and not isinstance(win.menuWidget(), QMenuBar)
    assert win.menuWidget() is win.top_area and win.ribbon.isAncestorOf(win.quick_access)


def test_header_row_is_aligned(app, quiet, sk_root, tmp_path, monkeypatch):
    monkeypatch.setenv("ADIT_FRAME", "native")
    win = make_window(sk_root, tmp_path); win.resize(1400, 820); win.show()
    for _ in range(5):
        app.processEvents()
    rb = win.ribbon
    heights = {w.height() for w in (win.quick_access, rb.tabs, rb.toggle)}
    tops = {w.mapTo(rb, w.rect().topLeft()).y() for w in (win.quick_access, rb.tabs, rb.toggle)}
    assert heights == {30} and len(tops) == 1, (heights, tops)
    assert all(b.size().width() == b.size().height() == 30 for b in win.quick_access.buttons)
    win.close()


def test_collapse_button_sits_next_to_the_tabs(app, quiet, sk_root, tmp_path, monkeypatch):
    monkeypatch.setenv("ADIT_FRAME", "native")
    win = make_window(sk_root, tmp_path); win.resize(1400, 820); win.show()
    for _ in range(5):
        app.processEvents()
    rb = win.ribbon
    gap = rb.toggle.mapTo(rb, rb.toggle.rect().topLeft()).x() - (rb.tabs.mapTo(rb, rb.tabs.rect().topRight()).x())
    assert 0 <= gap <= 24, gap
    assert rb.toggle.mapTo(rb, rb.toggle.rect().topRight()).x() < rb.width() / 2
    win.close()


def test_main_actions_are_at_the_bottom_right(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.resize(1400, 820); win.show()
    for _ in range(5):
        app.processEvents()
    bar = win.action_bar
    assert bar.mapTo(win, bar.rect().center()).y() > win.height() * 0.8
    assert win.btn_generate.mapTo(win, win.btn_generate.rect().center()).x() > win.width() * 0.6
    assert win.ribbon.isAncestorOf(bar) is False
    win.close()


def test_generate_hint_jumps_to_the_field_and_marks_it_red(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.method.reload_sets("")
    win.refresh_preview(); app.processEvents()
    assert win._error_locations == ["method.sk_set"]
    assert win.gen_hint.click() or True
    win.gen_hint.click()
    for _ in range(4):
        app.processEvents()
    assert win.mode() == win.MODE_SETTINGS
    marked = win._error_marked
    assert marked and all(w.property("adit_error") for w in marked)
    assert any(w.property("adit_key") == "Slater-Koster パラメータ" for w in marked)
    assert win.method.sk_set in marked
    win.method.reload_sets(str(sk_root)); win.method.sk_set.setCurrentIndex(0)
    win.refresh_preview(); app.processEvents()
    assert win._error_marked == [] and not win.method.sk_set.property("adit_error")
    win.close()


def test_generating_points_at_the_terminal(app, quiet, sk_root, tmp_path, monkeypatch):
    """実行ボタンは無い。生成したら、ターミナルで何を打つかを案内する (2026-09-15)。"""
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    win = make_window(sk_root, tmp_path)
    assert not hasattr(win, "btn_run") and not hasattr(win, "act_run")
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    win.generate()
    assert "bash submit.sh" in win.run_hint.text()
    assert win.workspace.root == Path(tmp_path / "out")
    win.runtime.profile.setCurrentText("cluster"); win.refresh_preview(); win.generate()
    assert "transfer_and_submit.sh" in win.run_hint.text()
    win.workspace.close_session()


def test_generate_waits_for_the_pending_preview(app, quiet, sk_root, tmp_path, monkeypatch):
    """出力ディレクトリを空にした直後 (250 ms の待ち中) に Ctrl+G を押しても、カレントディレクトリに書かない。"""
    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    cwd = tmp_path / "cwd"; cwd.mkdir(); monkeypatch.chdir(cwd)
    win.runtime.outdir.setText("")
    assert win._timer.isActive() and win.btn_generate.isEnabled()      # the preview has not caught up yet
    win.generate()
    assert sorted(p.name for p in cwd.iterdir()) == [] and not win.btn_generate.isEnabled()
    win.workspace.close_session()


def test_closing_the_window_asks_about_unsaved_edits(app, quiet, sk_root, tmp_path, monkeypatch):
    win = make_window(sk_root, tmp_path)
    note = tmp_path / "note.txt"; note.write_text("orig", encoding="utf-8")
    assert win.workspace.open_file(note) == ""
    win.workspace.editor.appendPlainText("more")
    assert win.workspace.editor.dirty
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.No))
    assert not win.close()                                             # kept open
    session = win.workspace.terminal.current.session
    assert session is not None and session.alive
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    assert win.close()
    assert note.read_text(encoding="utf-8") == "orig" and not session.alive   # discarded, shell ended


def test_ctrl_s_saves_the_file_in_the_workspace_and_the_spec_elsewhere(app, quiet, sk_root, tmp_path, monkeypatch):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest

    win = make_window(sk_root, tmp_path)
    calls = []
    monkeypatch.setattr(win, "save_spec", lambda: calls.append("spec"))
    monkeypatch.setattr(win.workspace, "save", lambda: calls.append("file") or "")
    win.act_save.triggered.disconnect(); win.act_save.triggered.connect(win.save_spec)
    from PySide6.QtGui import QShortcut
    for sc in win.workspace.findChildren(QShortcut):
        sc.activated.disconnect(); sc.activated.connect(win.workspace.save)
    win.show(); app.processEvents()
    win.set_mode(win.MODE_WORKSPACE); win.workspace.editor.setFocus(); app.processEvents()
    QTest.keyClick(win.workspace.editor, Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier); app.processEvents()
    assert calls == ["file"]
    calls.clear()
    win.set_mode(win.MODE_STRUCTURE); win.structure.setFocus(); app.processEvents()
    QTest.keyClick(win.structure, Qt.Key.Key_S, Qt.KeyboardModifier.ControlModifier); app.processEvents()
    assert calls == ["spec"]
    win.workspace.close_session()
