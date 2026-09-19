
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


# ---- proposal 9: one progress strip with the NN/g thresholds, and empty states with the next action ----

def test_progress_strip_follows_the_thresholds(app):
    from adit.gui.progress import EmptyState, ProgressStrip
    strip = ProgressStrip()
    now = [100.0]
    strip.clock = lambda: now[0]
    assert strip.phase() == "idle" and strip.isHidden()
    strip.begin("作っています…")
    strip._tick()
    assert strip.phase() == "hidden" and strip.isHidden(), "nothing under a second"
    now[0] += 1.5; strip._tick()
    assert strip.phase() == "busy" and not strip.isHidden() and strip.bar.maximum() == 0 and strip.btn_cancel.isVisibleTo(strip)
    strip.update(3, 10, "分子を置いています")
    assert strip.phase() == "busy" and strip.label.text() == "分子を置いています  3 / 10"
    now[0] += 10.0; strip._tick()
    assert strip.phase() == "percent" and strip.bar.maximum() == 100 and strip.bar.value() == 30
    assert "30 %" in strip.label.text() and "残り約" in strip.label.text()
    strip.end()
    assert strip.phase() == "idle" and strip.isHidden()
    strip.begin("x", cancellable=False); now[0] += 2; strip._tick()
    assert not strip.btn_cancel.isVisibleTo(strip)
    strip.end()
    empty = EmptyState("見出し", "1 行", "次へ")
    assert empty.button.isVisibleTo(empty) and empty.title.text() == "見出し"
    empty.set_texts(button="")
    assert not empty.button.isVisibleTo(empty)


def test_job_runs_in_the_background_and_drops_a_cancelled_result(app):
    from PySide6.QtCore import QEventLoop, QTimer

    from adit import progress as reports
    from adit.gui.progress import Job

    def spin(ms):
        loop = QEventLoop(); QTimer.singleShot(ms, loop.quit); loop.exec()

    got = []
    job = Job(); job.finished.connect(lambda r, e: got.append((r, e))); job.progress.connect(lambda d, t, w: got.append((d, t, w)))

    def work(n):
        reports.report(1, 2, "half")
        return n * 2

    job.start(work, 21); spin(300)
    assert (1, 2, "half") in got and (42, None) in got
    got.clear()
    job.start(lambda: 1 / 0); spin(300)
    assert len(got) == 1 and got[0][0] is None and isinstance(got[0][1], ZeroDivisionError)
    got.clear()
    import time
    job.start(lambda: time.sleep(0.2) or "late"); job.cancel(); spin(500)
    assert got == [] and not job.is_running()


def test_the_build_shows_the_strip_and_reports_steps(app, quiet, sk_root, tmp_path):
    from PySide6.QtCore import QEventLoop, QTimer

    def spin(ms):
        loop = QEventLoop(); QTimer.singleShot(ms, loop.quit); loop.exec()

    win = make_window(sk_root, tmp_path)
    p = win.structure
    p.set_source("mixture")
    seen = []
    p.job.progress.connect(lambda d, t, w: seen.append((d, t, w)))
    strip = p.recipe.progress
    now = [1000.0]; strip.clock = lambda: now[0]
    p.recipe.btn_build.click()
    assert p.is_building() and strip.phase() == "hidden" and strip.isHidden()
    now[0] += 2.0; strip._tick()
    assert strip.phase() == "busy"
    for _ in range(40):
        spin(100)
        if not p.is_building():
            break
    assert not p.is_building() and strip.phase() == "idle"
    assert any(w == "分子を置いています" for _, _, w in seen), seen[:3]
    assert p.structure() is not None
    win.close()


def test_analysis_runs_in_the_background_with_cancel(app, quiet, tmp_path):
    import shutil
    from pathlib import Path

    from PySide6.QtCore import QEventLoop, QTimer

    from adit.gui.panels.analysis_panel import AnalysisPanel

    def spin(ms):
        loop = QEventLoop(); QTimer.singleShot(ms, loop.quit); loop.exec()

    d = tmp_path / "run"
    shutil.copytree(Path(__file__).resolve().parent.parent / "examples" / "dftb_md_water_generated", d)
    p = AnalysisPanel()
    assert not p.empty.isHidden() and p.empty.button.text() == "解析を実行"
    p.set_run_dir(d, ["O", "H"])
    assert p.run_async()
    assert not p.btn_run.isEnabled() and p.job.is_running()
    p.cancel()
    assert p.btn_run.isEnabled() and not p.job.is_running() and p.summary.toPlainText() == "中止しました"
    assert p.run_async()
    for _ in range(100):
        spin(100)
        if not p.job.is_running():
            break
    assert not p.job.is_running() and p.empty.isHidden() and "温度" in p.summary.toPlainText() and p.figs_lay.count() >= 2
    p.run_dir.setText(str(tmp_path / "nowhere"))
    assert not p.run_async() and "ディレクトリがありません" in p.summary.toPlainText()


def test_empty_states_offer_the_next_action(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    ws = win.workspace
    assert ws.editor_stack.currentWidget() is ws.editor_empty and ws.editor_empty.title.text() == "ファイルを開いていません"
    (tmp_path / "out").mkdir(exist_ok=True); (tmp_path / "out" / "submit.sh").write_text("#!/bin/bash\n", encoding="utf-8")
    ws.set_root(tmp_path / "out")
    assert ws.editor_empty.button.text() == "submit.sh を開く"
    ws.editor_empty.button.click()
    assert ws.editor_stack.currentWidget() is ws.editor and ws.editor.path == tmp_path / "out" / "submit.sh"
    win.method.reload_sets("")
    win.refresh_preview()
    pv = win.preview
    assert pv.tabs.isHidden() and not pv.empty.isHidden() and pv.empty.title.text() == "生成できません"
    assert "Slater-Koster" in pv.empty.line.text() and pv.empty.button.text() == "欄へ移動"
    pv.empty.button.click()
    assert win.mode() == win.MODE_SETTINGS and win.method.sk_set.property("adit_error") is True
    win.method.reload_sets(str(sk_root)); win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    assert not pv.tabs.isHidden() and pv.empty.isHidden()
    dlg = win.scan_dialog()
    assert dlg.empty.button.text() == "例の値を入れる" and dlg.values.text() == ""
    dlg.empty.button.click()
    assert dlg.values.text() == dlg.current_choice().example and dlg.empty.title.text().endswith("個のディレクトリ")
    assert not dlg.empty.button.isVisibleTo(dlg)
    dlg.close(); win.close()
