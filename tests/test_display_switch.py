import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from adit import lang  # noqa: E402
from adit.config import default_config, load_config, save_config  # noqa: E402
from adit.gui import i18n  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.style import apply_theme, current_theme  # noqa: E402


@pytest.fixture
def window(tmp_path):
    app = QApplication.instance() or QApplication([])
    path = tmp_path / "cluster.toml"
    cfg = default_config(sk_root="")
    cfg.language, cfg.theme = "ja", "light"
    save_config(cfg, path)
    i18n.set_language("ja")
    apply_theme(app, "light")
    win = MainWindow(load_config(path), path)
    win.show()
    app.processEvents()
    yield win
    apply_theme(app, "light")


def _action(items, code):
    return dict(items)[code]


def _background(win) -> str:
    QApplication.instance().processEvents()
    image = win.grab().toImage()
    return image.pixelColor(image.width() // 2, image.height() - 40).name()


def test_the_default_language_is_english():
    assert default_config(sk_root="").language == "en"
    assert i18n.language_from_env() in ("en", os.environ.get("ADIT_LANG", "en"))


def test_a_settings_file_without_language_means_english(tmp_path):
    path = tmp_path / "cluster.toml"
    save_config(default_config(sk_root=""), path)
    text = "\n".join(line for line in path.read_text(encoding="utf-8").splitlines()
                     if not line.startswith("language"))
    path.write_text(text + "\n", encoding="utf-8")
    assert load_config(path).language == "en"


def test_the_language_changes_at_once_and_is_saved(window):
    assert window.btn_generate.text() == "生成"
    _action(window._lang_actions, "en").trigger()
    app = QApplication.instance(); app.processEvents()
    new = MainWindow._current
    assert new is not window and new.isVisible()
    assert new.btn_generate.text() == "Generate"
    assert lang.LANGUAGE == "en" and i18n.LANGUAGE == "en"
    assert load_config(new.cfg_path).language == "en"
    assert _action(new._lang_actions, "en").isChecked()

    _action(new._lang_actions, "ja").trigger()
    app.processEvents()
    back = MainWindow._current
    assert back.btn_generate.text() == "生成" and lang.LANGUAGE == "ja"
    assert load_config(back.cfg_path).language == "ja"


def test_the_theme_changes_at_once_and_is_saved(window):
    light = _background(window)
    _action(window._theme_actions, "dark").trigger()
    app = QApplication.instance(); app.processEvents()
    new = MainWindow._current
    assert current_theme() == "dark"
    assert _background(new) != light
    assert load_config(new.cfg_path).theme == "dark"
    assert _action(new._theme_actions, "dark").isChecked()

    _action(new._theme_actions, "light").trigger()
    app.processEvents()
    assert current_theme() == "light" and _background(MainWindow._current) == light


def test_switching_keeps_the_calculation_settings_and_the_tab(window):
    window.set_mode(window.MODE_SETTINGS)
    window.task.type.setCurrentIndex(window.task.type.findData("molecular_dynamics"))
    window.analysis.run_dir.setText(str(Path("some") / "run"))
    window.runtime.outdir.setText("my_output")
    window._set_language("en")
    new = MainWindow._current
    assert new.task.type.currentData() == "molecular_dynamics"
    assert new.mode() == new.MODE_SETTINGS
    assert new.runtime.outdir.text() == "my_output"
    assert new.analysis.run_dir.text() == str(Path("some") / "run")


def test_a_running_terminal_command_stops_the_switch_unless_confirmed(window, monkeypatch):
    monkeypatch.setattr(window.workspace.terminal, "any_busy", lambda: True)
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.No))
    assert window._change_appearance("theme", "dark") is window
    assert window.cfg.theme == "light" and current_theme() == "light"
    assert _action(window._theme_actions, "light").isChecked()
    assert load_config(window.cfg_path).theme == "light"


def test_the_window_frame_changes_at_once_and_is_saved(window):
    other = "native" if window.titlebar is not None else "custom"
    _action(window._frame_actions, other).trigger()
    new = MainWindow._current
    assert (new.titlebar is None) == (other == "native")
    assert load_config(new.cfg_path).window_frame == other


def test_no_message_is_translated_at_import_time():
    import ast

    root = Path(__file__).resolve().parents[1] / "src" / "adit"
    frozen = []

    def walk(node, path):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                for default in [*child.args.defaults, *[d for d in child.args.kw_defaults if d]]:
                    walk(ast.Expr(default), path)
                continue
            if isinstance(child, ast.Call) and getattr(child.func, "id", "") in ("L", "tr", "pick"):
                frozen.append(f"{path}:{child.lineno}")
            walk(child, path)

    for path in root.rglob("*.py"):
        walk(ast.parse(path.read_text(encoding="utf-8")), path)
    assert frozen == []


def test_choosing_the_current_value_does_not_rebuild(window):
    assert window._change_appearance("language", "ja") is window


def test_the_settings_dialog_applies_a_change(window):
    cfg = load_config(window.cfg_path)
    cfg.theme = "dark"
    save_config(cfg, window.cfg_path)
    window.reload_config()
    assert current_theme() == "dark" and MainWindow._current is not window
