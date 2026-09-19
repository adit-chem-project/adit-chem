
from __future__ import annotations

from adit.config import env_var
import os
from dataclasses import dataclass

from PySide6.QtGui import QFont, QFontDatabase
from PySide6.QtWidgets import QApplication

LABEL_WIDTH = 120       # floor; the widest label of a group sets its column
NARROW_FIELD = 120      # every numeric field
LABEL_GAP = 12          # between the label column and the fields
GROUP_SPACING = 14
ROW_SPACING = 8
PANEL_MARGIN = 14
RADIUS = 12
RADIUS_SMALL = 10
FONT_FAMILIES = ["Inter", "Helvetica Neue", "SF Pro Text", "Noto Sans CJK JP", "Noto Sans JP", "Segoe UI", "Hiragino Sans", "sans-serif"]


@dataclass(frozen=True)
class Tokens:

    bg: str
    card: str
    field: str
    line: str
    fg: str
    muted: str
    accent: str
    accent_hover: str
    ok: str
    ng: str
    bg_top: str = ""
    bg_bottom: str = ""
    # Surfaces are (nearly) opaque so that text contrast does not depend on what lies underneath;
    # only the ground keeps its gradient. The names are historical.
    glass: str = ""
    glass_field: str = ""
    glass_hover: str = ""
    glass_edge: str = ""
    accent_top: str = ""
    pill_bg: str = ""
    pill_fg: str = ""


LIGHT = Tokens(bg="#EFEFF1", card="#FFFFFF", field="#F4F4F6", line="rgba(0, 0, 0, 30)", fg="#14161B", muted="#4B5260",
               accent="#0A7AFF", accent_hover="#0866D6", ok="#1F7A38", ng="#C8001A",
               bg_top="#E9E9EC", bg_bottom="#F5F5F7", glass="rgba(255, 255, 255, 242)", glass_field="rgba(250, 250, 252, 250)",
               glass_hover="rgba(255, 255, 255, 255)", glass_edge="rgba(255, 255, 255, 230)", accent_top="#3A95FF",
               pill_bg="#E1E3E8", pill_fg="#3B4150")
DARK = Tokens(bg="#1C1C1E", card="#2C2C2E", field="#3A3A3C", line="rgba(255, 255, 255, 34)", fg="#F0F1F4", muted="#B6B8C2",
              accent="#3D8CFF", accent_hover="#5C9EFF", ok="#5FD07A", ng="#FF6961",
              bg_top="#1C1C1E", bg_bottom="#2A2A2D", glass="rgba(48, 48, 51, 242)", glass_field="rgba(62, 62, 66, 250)",
              glass_hover="rgba(78, 78, 83, 255)", glass_edge="rgba(255, 255, 255, 38)", accent_top="#6AA8FF",
              pill_bg="#4A4B52", pill_fg="#E8E9EE")

GROUP_COLORS = {"structure": "#0FA3A3", "task": "#F0883E", "method": "#8E5AD6", "kpoints": "#34A853", "runtime": "#5B7DB1", "analysis": "#D6567A"}

THEME_COLORS = {
    "[light]": {"background": LIGHT.bg, "border": "#C9CDD8", "input.background": LIGHT.field, "primary": LIGHT.accent, "foreground": LIGHT.fg},
    "[dark]": {"background": DARK.bg, "border": "#454852", "input.background": DARK.field, "primary": DARK.accent, "foreground": DARK.fg},
}


def qss(t: Tokens) -> str:
    r, rs = RADIUS, RADIUS_SMALL
    group_css = "\n".join(f"QGroupBox#{k} {{ border-left: 3px solid {v}; }}" for k, v in GROUP_COLORS.items())
    ground = f"qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 {t.bg_top}, stop:1 {t.bg_bottom})"
    primary = f"qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 {t.accent_top}, stop:1 {t.accent})"
    return f"""
QMainWindow, QDialog {{ background: {ground}; }}
/* 地のグラデーションを透かすため、土台の部品 (中央の部品・分割・スクロールの中・タブの中身) は塗らない。
   既製テーマ (qdarktheme) はこれらを不透明に塗るので、ここで打ち消す */
QMainWindow > QWidget, QSplitter, QSplitter > QWidget, QStackedWidget, QStackedWidget > QWidget,
QScrollArea > QWidget#qt_scrollarea_viewport, QAbstractScrollArea#qt_scrollarea_viewport {{ background: transparent; }}
/* 本体の窓を親にしたダイアログも「QMainWindow > QWidget」に当たって透明になる (> は親オブジェクトで判定する)。地を塗り直す */
QMainWindow > QDialog {{ background: {ground}; }}
QScrollArea {{ border: none; background: transparent; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}
QSplitter::handle {{ background: transparent; width: {GROUP_SPACING}px; }}

QGroupBox {{
    font-weight: 700;
    font-size: 12.5pt;
    background: {t.glass};
    border: 1px solid {t.glass_edge};
    border-left: 3px solid {t.line};
    border-radius: {r}px;
    margin-top: 0px;
    padding: 34px {PANEL_MARGIN}px 12px {PANEL_MARGIN}px;
}}
/* 見出しは面 (カード) の中の左上に置く。境界線の上にまたがらせない */
QGroupBox::title {{
    subcontrol-origin: padding;
    subcontrol-position: top left;
    left: {PANEL_MARGIN}px;
    top: 10px;
    padding: 0;
    color: {t.fg};
    background: transparent;
}}
/* 群のタイトルだけ太字にし、中の部品は通常の字に戻す */
QGroupBox QWidget {{ font-weight: 400; font-size: 11pt; color: {t.fg}; }}
{group_css}
QWidget#rowbox {{ background: transparent; }}
QLabel {{ background: transparent; }}
QLabel#hint {{ color: {t.muted}; font-size: 10pt; }}
QLabel#status_ok {{ color: {t.ok}; font-weight: 600; }}
QLabel#status_ng {{ color: {t.ng}; font-weight: 600; }}
QLabel#title {{ font-weight: 700; font-size: 12.5pt; }}
/* Sub-headings inside a card are set apart by space above, not by weight */
QLabel#subtitle {{ font-weight: 500; font-size: 11pt; color: {t.muted}; margin-top: 16px; }}
QLabel#unit {{ color: {t.muted}; }}
QLabel#required {{ color: {t.accent}; font-weight: 600; }}
/* The word next to a required label; the color alone is not enough (also used by links) */
QLabel#required_pill {{ background: {t.pill_bg}; color: {t.pill_fg}; border-radius: 9px; font-size: 8.5pt; font-weight: 600; }}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QPlainTextEdit, QTextEdit {{
    background: {t.glass_field};
    border: 1px solid {t.line};
    border-radius: {rs}px;
    padding: 3px 8px;
    min-height: 24px;
    selection-background-color: {t.accent};
}}
QPlainTextEdit, QTextEdit {{ padding: 8px 10px; }}
QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover {{ background: {t.glass_hover}; }}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QPlainTextEdit:focus {{ border: 1px solid {t.accent}; background: {t.glass_hover}; }}
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled {{ color: {t.muted}; background: transparent; }}
QComboBox::drop-down {{ border: none; width: 18px; }}
/* プルダウンの一覧は別の窓なので、読みやすさを優先して不透明。
   combobox-popup: 0 で「選んだ項目の位置に画面いっぱいの一覧を重ねる」出し方を止め、一度に見える行数 (widgets.POPUP_ROWS) で切って
   スクロールさせる (選択肢の多い一覧が画面の上下で見切れていた。2026-09-12 の方針) */
QComboBox {{ combobox-popup: 0; }}
QComboBox QAbstractItemView {{ background: {t.card}; border: 1px solid {t.line}; border-radius: 8px; selection-background-color: {t.accent}; }}
QComboBox QAbstractItemView QScrollBar:vertical {{ background: transparent; width: 12px; margin: 4px 2px 4px 0; border: none; }}
QComboBox QAbstractItemView QScrollBar::handle:vertical {{ background: {t.muted}; border-radius: 4px; min-height: 28px; margin: 0 2px; }}
QComboBox QAbstractItemView QScrollBar::add-line:vertical, QComboBox QAbstractItemView QScrollBar::sub-line:vertical {{ height: 0; border: none; background: none; }}
QComboBox QAbstractItemView QScrollBar::add-page:vertical, QComboBox QAbstractItemView QScrollBar::sub-page:vertical {{ background: none; }}
QSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{ border: none; width: 16px; }}

QPushButton {{
    background: {t.glass_field};
    border: 1px solid {t.line};
    border-radius: {rs}px;
    padding: 3px 14px;
    min-height: 24px;
    font-weight: 500;
}}
QPushButton:hover {{ background: {t.glass_hover}; border: 1px solid {t.glass_edge}; }}
QPushButton:pressed {{ background: {t.line}; }}
QPushButton:disabled {{ color: {t.muted}; background: transparent; }}
QPushButton#primary {{
    background: {primary};
    color: #FFFFFF;
    border: 1px solid {t.accent};
    padding: 5px 20px;
    font-weight: 600;
}}
QPushButton#primary:hover {{ background: {t.accent_hover}; }}
QPushButton#primary:disabled {{ background: {t.glass_field}; color: {t.muted}; border: 1px solid {t.line}; }}
QPushButton#cell_button {{ padding: 2px 6px; min-height: 20px; }}
QTableWidget QSpinBox, QTableWidget QComboBox, QTableWidget QLineEdit {{ padding: 1px 6px; min-height: 18px; }}
QPushButton#link {{ border: none; background: transparent; color: {t.accent}; padding: 2px 6px; }}
/* Underlined so that a link-style button is not mistaken for a required label */
QPushButton#link {{ text-decoration: underline; }}
QPushButton#link:hover {{ text-decoration: underline; background: transparent; border: none; }}
/* Number of problems, left of Generate; the full list opens from it */
QPushButton#error_badge {{ background: {t.ng}; color: #FFFFFF; border: none; border-radius: 11px; padding: 2px 12px; min-height: 18px; font-weight: 600; }}
QPushButton#error_badge:hover {{ background: {t.ng}; color: #FFFFFF; border: none; }}
/* The reason, in red, under the field it belongs to; the row is hidden when there is no error */
QLabel#field_error {{ color: {t.ng}; font-size: 10pt; }}
/* Value differs from the code default (labels paint the same bar themselves, see widgets.FieldLabel) */
QCheckBox[adit_changed="user"] {{ border-left: 2px solid {t.accent}; padding-left: 6px; }}
QCheckBox[adit_changed="template"] {{ border-left: 2px solid {t.muted}; padding-left: 6px; }}

QToolBar {{ background: transparent; border: none; spacing: 8px; padding: 8px {PANEL_MARGIN}px 4px {PANEL_MARGIN}px; }}
QToolBar QToolButton {{ border-radius: {rs}px; padding: 5px 12px; background: transparent; }}
QToolBar QToolButton:hover {{ background: {t.glass_hover}; }}
QToolBar QToolButton:checked {{ background: {t.glass}; border: 1px solid {t.line}; }}
QToolBar QLabel {{ color: {t.muted}; font-size: 10pt; }}
/* 自前の題名の帯 (gui/titlebar.py)。地のグラデーションを見せる。枠を外した窓には 1 px の縁を付け、隣の窓と見分けられるようにする */
QMainWindow#adit_frameless {{ border: 1px solid {t.line}; }}
QWidget#titlebar {{ background: transparent; border-bottom: 1px solid {t.line}; }}
QWidget#titlebar QMenuBar {{ background: transparent; border: none; }}
/* リボン (gui/ribbon.py)。見出しは右側のタブと同じ丸い札、中身は群と同じすりガラスの面 */
QWidget#ribbon, QWidget#qat, QWidget#ribbon_group, QStackedWidget#ribbon_stack {{ background: transparent; }}
/* メインの面の札 (gui/ribbon.py の ModeBar)。作業の順に左から並ぶ。選ばれている札だけ色を持つ */
QWidget#modebar {{ background: transparent; }}
QToolButton#mode_button {{
    background: transparent; border: 1px solid transparent; border-radius: 8px;
    padding: 4px 14px; font-weight: 600; color: {t.muted};
}}
QToolButton#mode_button:hover {{ background: {t.glass_hover}; }}
QToolButton#mode_button:checked {{ background: {t.glass}; border: 1px solid {t.glass_edge}; color: {t.fg}; }}
QTabBar#ribbon_tabs::tab {{ padding: 4px 11px; margin: 2px 2px 2px 0; font-weight: 600; min-height: 18px; }}
QWidget#ribbon_page {{ background: {t.glass}; border: 1px solid {t.glass_edge}; border-radius: {r}px; }}
QFrame#ribbon_sep {{ background: {t.line}; border: none; margin: 6px 0 6px 0; }}
QLabel#ribbon_caption {{ color: {t.muted}; font-size: 8.5pt; }}
QToolButton#ribbon_large, QToolButton#ribbon_small, QToolButton#qat_button, QToolButton#ribbon_toggle {{
    background: transparent; border: 1px solid transparent; border-radius: 8px;
}}
QToolButton#ribbon_large {{ padding: 2px 6px 1px 6px; min-width: 48px; }}
QToolButton#ribbon_small {{ padding: 1px 6px; }}
QToolButton#qat_button, QToolButton#ribbon_toggle {{ padding: 4px; border-radius: 6px; min-width: 18px; min-height: 18px; }}
QToolButton#ribbon_large:hover, QToolButton#ribbon_small:hover, QToolButton#qat_button:hover, QToolButton#ribbon_toggle:hover {{
    background: {t.glass_hover}; border: 1px solid {t.glass_edge};
}}
QToolButton#ribbon_large:pressed, QToolButton#ribbon_small:pressed, QToolButton#qat_button:pressed {{ background: {t.line}; }}
QToolButton#ribbon_large:checked, QToolButton#ribbon_small:checked {{ background: {t.glass_hover}; border: 1px solid {t.accent}; }}
QToolButton#ribbon_large:disabled, QToolButton#ribbon_small:disabled, QToolButton#qat_button:disabled {{ color: {t.muted}; }}
QToolBar#ribbon_actions {{ padding: 0 0 0 8px; spacing: 8px; }}
QMenuBar {{ background: transparent; }}
QMenuBar::item {{ background: transparent; padding: 4px 10px; border-radius: 6px; }}
QMenuBar::item:selected {{ background: {t.glass_hover}; }}
QMenu {{ background: {t.card}; border: 1px solid {t.line}; border-radius: 10px; padding: 6px; }}
QMenu::item {{ padding: 5px 22px 5px 12px; border-radius: 6px; }}
QMenu::item:selected {{ background: {t.accent}; color: #FFFFFF; }}
QLabel#window_title {{ color: {t.muted}; font-size: 10pt; }}
QStatusBar {{ background: transparent; color: {t.muted}; border: none; }}
QStatusBar::item {{ border: none; }}

QTabWidget::pane {{ border: 1px solid {t.glass_edge}; border-radius: {r}px; background: {t.glass}; top: -1px; }}
QTabWidget > QWidget {{ background: transparent; }}
QTabBar {{ qproperty-drawBase: 0; background: transparent; }}
QTabBar::tab {{
    background: transparent;
    border: none;
    border-radius: {rs}px;
    padding: 6px 14px;
    margin: 0 4px 6px 0;
    color: {t.muted};
    font-weight: 500;
}}
QTabBar::tab:selected {{ background: {t.glass_hover}; color: {t.fg}; border: 1px solid {t.glass_edge}; }}
QTabBar::tab:hover:!selected {{ background: {t.glass}; }}
QTableWidget, QTableView {{ background: {t.glass_field}; border: 1px solid {t.line}; border-radius: 8px; gridline-color: {t.line}; }}
QListWidget {{ background: {t.glass_field}; border: 1px solid {t.line}; border-radius: 8px; padding: 4px; }}
QListWidget::item {{ padding: 4px 6px; border-radius: 6px; }}
/* 選んだ行は、窓の前後 (active / inactive) に依らず地の文字色のまま、縁をアクセント色に (ライトで白い字が淡い地に乗り、読めなかった) */
QListWidget::item:selected, QListWidget::item:selected:!active {{ background: {t.glass_hover}; color: {t.fg}; border: 1px solid {t.accent}; }}
QHeaderView::section {{ background: {t.glass}; border: none; border-bottom: 1px solid {t.line}; padding: 4px 6px; color: {t.muted}; }}
QProgressBar {{ background: {t.glass_field}; border: 1px solid {t.line}; border-radius: 6px; min-height: 10px; max-height: 10px; }}
QProgressBar::chunk {{ background: {primary}; border-radius: 5px; }}
QCheckBox, QRadioButton {{ spacing: 8px; background: transparent; }}
QCheckBox::indicator, QRadioButton::indicator {{ width: 16px; height: 16px; }}
QToolTip {{ background: {t.card}; color: {t.fg}; border: 1px solid {t.line}; border-radius: 6px; padding: 4px 8px; }}

/* 生成できない欄の印 (2026-09-13 の方針)。理由の行を押すと、その欄へ飛んでここが赤くなる。
   **このまとまりは一番最後に置く。**前に書くと、:focus などあとの規則に上書きされて色が出ない */
QComboBox[adit_error="true"], QLineEdit[adit_error="true"], QSpinBox[adit_error="true"],
QDoubleSpinBox[adit_error="true"], QPlainTextEdit[adit_error="true"], QPushButton[adit_error="true"],
QComboBox[adit_error="true"]:focus, QLineEdit[adit_error="true"]:focus, QSpinBox[adit_error="true"]:focus,
QDoubleSpinBox[adit_error="true"]:focus, QPlainTextEdit[adit_error="true"]:focus {{
    border: 2px solid {t.ng};
}}
QLabel[adit_error="true"], QLabel#required[adit_error="true"] {{ color: {t.ng}; font-weight: 700; }}
QWidget[adit_error="true"] {{ border: 2px solid {t.ng}; border-radius: {rs}px; }}
"""


def resolve_theme(theme: str | None = None) -> str:
    theme = (theme or env_var("THEME", "auto")).lower()
    if theme in ("light", "dark"):
        return theme
    try:
        import darkdetect

        return "dark" if (darkdetect.theme() or "").lower() == "dark" else "light"
    except Exception:
        return "light"


def tokens_for(theme: str) -> Tokens:
    return DARK if resolve_theme(theme) == "dark" else LIGHT


def apply_font(app: QApplication) -> str:
    import sys
    from pathlib import Path

    cand = Path(sys.prefix) / "fonts" / "NotoSansCJKjp-VF.ttf"
    if cand.is_file():
        QFontDatabase.addApplicationFont(str(cand))
    available = set(QFontDatabase.families())
    families = [f for f in FONT_FAMILIES if f in available] or [app.font().family()]
    font = QFont(); font.setFamilies(families); font.setPointSize(11)
    app.setFont(font)
    return families[0]


_APPLIED = ""


def current_theme() -> str:
    return _APPLIED or resolve_theme()


def _tell_icons(theme: str) -> None:
    from adit.gui import icons

    icons.set_theme(theme)


def apply_base_style(app: QApplication, theme: str | None = None) -> None:
    global _APPLIED
    resolved = resolve_theme(theme)
    _APPLIED = resolved
    _tell_icons(resolved)
    app.setStyleSheet(qss(tokens_for(resolved)))


def apply_theme(app: QApplication, theme: str | None = None) -> str:
    global _APPLIED
    theme = theme or env_var("THEME", "auto")
    resolved = resolve_theme(theme)
    _APPLIED = resolved
    _tell_icons(resolved)
    try:
        import qdarktheme

        qdarktheme.setup_theme(resolved, custom_colors=THEME_COLORS, additional_qss=qss(tokens_for(resolved)), default_theme="light")
        return f"qdarktheme:{resolved}"
    except Exception as ex:
        apply_base_style(app, resolved)
        return f"base:{resolved} (qdarktheme を使えない: {type(ex).__name__}: {ex})"
