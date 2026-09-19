
import os
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from adit import lang  # noqa: E402
from adit.config import load_config, save_config  # noqa: E402
from adit.gui.settings_dialog import SettingsDialog, set_top_level  # noqa: E402
from tests.conftest import cfg_for  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def shown(monkeypatch):
    got = []
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda _p, title, text, *a, **k: got.append((title, text))))
    return got


@pytest.fixture
def ja():
    old = lang.LANGUAGE
    lang.set_language("ja")
    yield
    lang.set_language(old)


def write_cfg(tmp_path: Path, head: str = "") -> Path:
    p = tmp_path / "cluster.toml"
    save_config(cfg_for(tmp_path / "sk"), p)
    if head:
        p.write_text(head + p.read_text(encoding="utf-8"), encoding="utf-8")
    return p


def test_fields_save_and_reload(app, tmp_path, shown, ja):
    p = write_cfg(tmp_path, "# 自分で書いたメモ\n")
    sk, upf = tmp_path / "my_sk", tmp_path / "my_upf"
    sk.mkdir(); upf.mkdir()
    dlg = SettingsDialog(p)
    assert dlg.tabs.currentIndex() == 0
    assert dlg.profiles.rowCount() == 3  # local / cluster / slurm
    kinds = {dlg.profiles.item(r, 0).text().split(" ")[0]: dlg.profiles.item(r, 1).text() for r in range(3)}
    assert "direct" in kinds["local"] and "pbs" in kinds["cluster"] and "slurm" in kinds["slurm"]
    dlg.sk_root.setText(str(sk)); dlg.pseudo_root.setText(str(upf))
    for key, value in (("language", "en"), ("theme", "dark"), ("window_frame", "native")):
        cb = dlg.combos[key]; cb.setCurrentIndex(cb.findData(value))
    dlg.default_profile.setCurrentIndex(dlg.default_profile.findData("cluster"))
    dlg._save()
    assert dlg.result() == SettingsDialog.DialogCode.Accepted and not shown
    cfg = load_config(p)
    assert (cfg.sk_root, cfg.pseudo_root) == (str(sk), str(upf))
    assert (cfg.language, cfg.theme, cfg.window_frame, cfg.default_profile) == ("en", "dark", "native", "cluster")
    assert set(cfg.profiles) == {"local", "cluster", "slurm"}
    assert "# 自分で書いたメモ" in p.read_text(encoding="utf-8")
    assert not p.with_suffix(".toml.tmp").exists()


def test_broken_toml_is_not_saved_and_told_in_japanese(app, tmp_path, shown, ja):
    p = write_cfg(tmp_path)
    before = p.read_text(encoding="utf-8")
    dlg = SettingsDialog(p)
    dlg.tabs.setCurrentIndex(1)
    dlg.editor.setPlainText("sk_root = abc\n")
    dlg._save()
    assert dlg.result() != SettingsDialog.DialogCode.Accepted
    assert p.read_text(encoding="utf-8") == before
    assert not p.with_suffix(".toml.tmp").exists()
    assert len(shown) == 1
    title, text = shown[0]
    assert title == "保存できません"
    assert "環境設定ファイルの 1 行目 11 文字目が TOML の書き方として正しくありません" in text
    assert ".tmp" not in text and "cluster.toml" not in text


def test_valid_toml_with_bad_value_is_not_saved(app, tmp_path, shown, ja):
    p = write_cfg(tmp_path)
    before = p.read_text(encoding="utf-8")
    dlg = SettingsDialog(p)
    dlg.tabs.setCurrentIndex(1)
    dlg.editor.setPlainText('[profiles.x]\nkind = "sge"\n')
    dlg._save()
    assert p.read_text(encoding="utf-8") == before
    assert "profiles.x.kind" in shown[0][1] and "値が使えません" in shown[0][1]


def test_tabs_carry_values_both_ways(app, tmp_path, shown, ja):
    p = write_cfg(tmp_path)
    dlg = SettingsDialog(p)
    dlg.pseudo_root.setText("/data/upf")
    dlg.tabs.setCurrentIndex(1)
    assert 'pseudo_root = "/data/upf"' in dlg.editor.toPlainText()
    dlg.editor.setPlainText(dlg.editor.toPlainText().replace('theme = "auto"', 'theme = "light"'))
    dlg.tabs.setCurrentIndex(0)
    assert dlg.combos["theme"].currentData() == "light" and dlg.pseudo_root.text() == "/data/upf"
    dlg.tabs.setCurrentIndex(1)
    dlg.editor.setPlainText("theme = \n")
    dlg.tabs.setCurrentIndex(0)
    assert dlg.tabs.currentIndex() == 1 and "1 行目" in shown[-1][1]
    assert not p.with_suffix(".toml.tmp").exists()


def test_broken_file_opens_on_advanced_tab(app, tmp_path, shown, ja):
    p = tmp_path / "cluster.toml"
    p.write_text("sk_root = abc\n", encoding="utf-8")
    dlg = SettingsDialog(p)
    assert dlg.tabs.currentIndex() == 1 and "1 行目 11 文字目" in dlg.raw_note.text()


def test_set_top_level_keeps_other_lines():
    text = '# メモ\nsk_root = "/a"  # 行末のメモ\ntheme = "auto"\n\n[profiles.local]\nkind = "direct"\n'
    out = set_top_level(text, {"sk_root": "/a", "theme": "dark", "enable_run": False})
    assert '# メモ' in out and 'sk_root = "/a"  # 行末のメモ' in out
    assert 'theme = "dark"' in out
    assert out.index("enable_run = false") < out.index("[profiles.local]")


def test_english_buttons(app, tmp_path):
    old = lang.LANGUAGE
    lang.set_language("en")
    try:
        dlg = SettingsDialog(write_cfg(tmp_path))
        assert (dlg.btn_save.text(), dlg.btn_cancel.text()) == ("Save", "Cancel")
        assert dlg.tabs.tabText(1).startswith("Advanced")
    finally:
        lang.set_language(old)
    dlg2 = SettingsDialog(write_cfg(tmp_path))
    assert (dlg2.btn_save.text(), dlg2.btn_cancel.text()) == ("保存", "キャンセル")
