
from __future__ import annotations

import pytest

from adit.config import default_config, load_config

from tests.conftest import pbs_profile

pytest.importorskip("PySide6")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def quiet(monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: None))


def make_window(sk_root, tmp_path, **cfg_values):
    from adit.gui.main_window import MainWindow
    cfg = default_config(sk_root=str(sk_root)); cfg.profiles["cluster"] = pbs_profile()
    for k, v in cfg_values.items():
        setattr(cfg, k, v)
    win = MainWindow(cfg, tmp_path / "cluster.toml")
    win.method.sk_set.setCurrentText("fake-1-0")
    win.runtime.outdir.setText(str(tmp_path / "out"))
    win.refresh_preview()
    return win


def settle(app, n: int = 5) -> None:
    for _ in range(n):
        app.processEvents()


# ---- proposal 1: the ribbon starts collapsed, opens as one row, and remembers the choice ----

def test_ribbon_starts_collapsed_and_opens_as_one_row(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.resize(1400, 820); win.show(); settle(app)
    rb = win.ribbon
    assert rb.is_collapsed() and not rb.stack.isVisible()
    h_closed = rb.height()
    rb.tabs.tabBarClicked.emit(0); settle(app)
    assert not rb.is_collapsed() and rb.is_peeking()
    assert rb.height() - h_closed <= 48, (rb.height(), h_closed)       # one row of 16 px icons, not two rows of 32 px
    g = rb.pages[0].groups[0]
    assert g.caption.isHidden() and g.toolTip() == "計算設定"
    tops = {b.mapTo(rb, b.rect().topLeft()).y() for b in rb.pages[0].groups[1].buttons}
    assert len(tops) == 1, tops                                          # the group's buttons sit on one line
    win.close()


def test_ribbon_peek_closes_after_a_command_and_the_toggle_pins_it(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.show(); settle(app)
    rb = win.ribbon
    view = rb.page_titles().index("表示")
    rb.tabs.tabBarClicked.emit(view); rb.tabs.setCurrentIndex(view); settle(app)
    assert rb.is_peeking()
    rb.pages[view].groups[0].buttons[2].click()                          # a command on the peeked page
    assert win.mode() == win.MODE_ANALYSIS and rb.is_collapsed()
    rb.tabs.tabBarClicked.emit(view); settle(app)
    rb.toggle.click()                                                    # "keep it open"
    assert not rb.is_collapsed() and not rb.is_peeking()
    rb.pages[view].groups[0].buttons[0].click()
    assert not rb.is_collapsed(), "a pinned ribbon stays open after a command"
    win.close()


def test_ribbon_state_is_saved_in_the_settings(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.show(); settle(app)
    win.ribbon.toggle.click()
    assert not win.ribbon.is_collapsed()
    assert load_config(tmp_path / "cluster.toml").ribbon_collapsed is False
    win.ribbon.toggle.click()
    assert load_config(tmp_path / "cluster.toml").ribbon_collapsed is True
    win.close()
    win2 = make_window(sk_root, tmp_path, ribbon_collapsed=False)
    assert not win2.ribbon.is_collapsed() and win2.ribbon.stack.isVisibleTo(win2)
    win2.close()


# ---- proposal 2: the command palette, shortcut tooltips and the shortcut list ----

def test_palette_search_scores_prefix_over_substring_and_needs_every_word():
    from adit.gui.palette import Item, search
    items = [Item("action", "生成", "Ctrl+G", lambda: None, ("generate",)),
             Item("field", "温度 [K]", "計算条件 › 計算の種類", lambda: None, ("温度 [K]", "Temperature [K]")),
             Item("field", "電子温度 [K]", "計算条件 › 計算手法", lambda: None, ("電子温度 [K]", "Electronic temperature [K]")),
             Item("preset", "H₂O", "プリセット", lambda: None, ("H2O H₂O water 水",))]
    assert [i.title for i in search(items, "温度")] == ["温度 [K]", "電子温度 [K]"]
    assert [i.title for i in search(items, "temp")] == ["温度 [K]", "電子温度 [K]"]     # a word start beats a substring
    assert [i.title for i in search(items, "gen")] == ["生成"]
    assert [i.title for i in search(items, "water")] == ["H₂O"] and [i.title for i in search(items, "水")] == ["H₂O"]
    assert search(items, "電子 温度")[0].title == "電子温度 [K]" and search(items, "電子 圧力") == []
    assert [i.kind for i in search(items, "")] == ["action", "field", "field", "preset"]   # commands first when nothing is typed


def test_palette_jumps_to_a_field_runs_a_command_and_picks_a_preset(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.show(); settle(app)
    win.set_mode(win.MODE_SETTINGS)
    dlg = win.open_palette(); settle(app)
    dlg.search.setText("スピン")
    hits = dlg.shown()
    assert hits and hits[0].kind == "field" and hits[0].title == "スピン多重度" and hits[0].detail.startswith("構造")
    assert dlg.run_current() and win.mode() == win.MODE_STRUCTURE
    assert win.focusWidget() is win.structure.multiplicity
    dlg = win.open_palette(); dlg.search.setText("一括")
    assert dlg.shown()[0].kind == "action" and dlg.shown()[0].detail == "" and "一括生成" in dlg.shown()[0].title
    dlg.reject()
    dlg = win.open_palette(); dlg.search.setText("エタノール")
    assert [i.kind for i in dlg.shown()] == ["preset"]
    assert dlg.run_current() and win.structure.current_source() == "preset" and win.structure.preset.currentData() == "CH3CH2OH"
    dlg = win.open_palette(); dlg.search.setText("undo")
    assert dlg.shown()[0].detail == "Ctrl+Z", "the shortcut is the detail of a command"
    dlg.reject()
    win.close()


def test_tooltips_carry_the_shortcut_once(app, quiet, sk_root, tmp_path):
    from adit.gui.i18n import set_language, translate_widgets
    win = make_window(sk_root, tmp_path)
    assert win.act_generate.toolTip() == "生成 (Ctrl+G)"
    assert win.act_open.toolTip() == "計算設定 (spec.json) を開く… (Ctrl+O)"
    assert win.act_back.toolTip() == "1 つ前の設定に戻す (Ctrl+Z)"
    assert win.ribbon.act_toggle.toolTip().count("Ctrl+F1") == 1
    assert win.act_palette.toolTip().endswith("(Ctrl+K)") and win.act_shortcuts.toolTip().endswith("(Ctrl+/)")
    win.close()
    set_language("en")
    try:
        win = make_window(sk_root, tmp_path); translate_widgets(win)
        assert win.act_open.toolTip() == "Open calculation settings (spec.json)… (Ctrl+O)"
        assert win.act_back.toolTip() == "Back to the previous settings (Ctrl+Z)"
        win.close()
    finally:
        set_language("ja")


def test_shortcut_list_names_every_key(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    rows = dict(win.shortcut_rows())
    for key in ("Ctrl+K", "Ctrl+/", "Ctrl+G", "Ctrl+Z", "Ctrl+Shift+Z", "Ctrl+F1", "Ctrl+S", "Ctrl+O", "Ctrl+Q"):
        assert key in rows, key
    dlg = win.open_shortcuts()
    assert dlg.rows() == win.shortcut_rows()
    dlg.close(); win.close()
