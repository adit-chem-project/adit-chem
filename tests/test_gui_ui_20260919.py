
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
