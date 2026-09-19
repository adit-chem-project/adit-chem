
from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication, QMessageBox

from adit.gui.style import apply_theme
from adit.config import ConfigError, config_path, default_config, load_config, save_config


def _ensure_japanese_font(app: QApplication) -> None:
    from adit.gui.style import apply_font

    apply_font(app)


def _prefer_xcb_on_wsl() -> None:
    import os

    if os.environ.get("QT_QPA_PLATFORM") or sys.platform != "linux":
        return
    try:
        is_wsl = "microsoft" in Path("/proc/version").read_text().lower()
    except OSError:
        return
    if is_wsl and os.environ.get("DISPLAY"):
        os.environ["QT_QPA_PLATFORM"] = "xcb"


def main() -> int:
    _prefer_xcb_on_wsl()
    app = QApplication(sys.argv)
    _ensure_japanese_font(app)
    from adit.config import ensure_config, first_run_message
    from adit.lang import L
    try:
        cfg, cfg_file, created = ensure_config(config_path())
    except (ConfigError, OSError) as ex:
        QMessageBox.critical(None, L("環境設定ファイルを読めません", "Cannot read the settings file"),
                             str(ex) + L("\n\nファイルは書き換えていません。直してから、もう一度起動してください。",
                                         "\n\nThe file was not changed. Fix it and start ADIT again."))
        return 2
    if created:
        QMessageBox.information(None, L("環境設定ファイルを作りました", "Settings file created"), first_run_message(cfg_file))
    from adit.config import env_var
    print("theme:", apply_theme(app, env_var("THEME") or cfg.theme), file=sys.stderr)
    from adit.gui.main_window import MainWindow

    from adit.gui.i18n import language_from_env, set_language, translate_widgets

    set_language(language_from_env(cfg.language))
    win = MainWindow(cfg, config_path())
    translate_widgets(win)
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
