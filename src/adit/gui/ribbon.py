
from __future__ import annotations

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (QFrame, QHBoxLayout, QLabel, QSizePolicy, QStackedWidget, QTabBar, QToolButton,
                               QVBoxLayout, QWidget)

from adit.gui import icons
from adit.lang import L

LARGE_ICON = 16
SMALL_ICON = 16
HEADER_HEIGHT = 30
PAGE_PADDING = 5


def _bind(button: QToolButton, action: QAction) -> QToolButton:
    button.setDefaultAction(action)
    button.setVisible(action.isVisible())
    action.changed.connect(lambda b=button, a=action: b.setVisible(a.isVisible()))
    return button


def _row_button(action: QAction, name: str) -> QToolButton:
    # One row of 16 px icons with the label beside: the ribbon body stays about 40 px tall.
    b = QToolButton(); b.setObjectName(name)
    b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
    b.setIconSize(QSize(SMALL_ICON, SMALL_ICON)); b.setAutoRaise(True)
    b.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
    return _bind(b, action)


def large_button(action: QAction) -> QToolButton:
    return _row_button(action, "ribbon_large")


def small_button(action: QAction) -> QToolButton:
    return _row_button(action, "ribbon_small")


class RibbonGroup(QWidget):

    triggered = Signal()

    def __init__(self, title: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("ribbon_group")
        self.body = QHBoxLayout(); self.body.setContentsMargins(0, 0, 0, 0); self.body.setSpacing(2)
        # The caption is kept for the group name (tooltip, translation) but no longer takes a row.
        self.caption = QLabel(title); self.caption.setObjectName("ribbon_caption"); self.caption.hide()
        self.setToolTip(title)
        lay = QVBoxLayout(self); lay.setContentsMargins(4, 0, 4, 0); lay.setSpacing(0)
        lay.addLayout(self.body, 1); lay.addWidget(self.caption)
        self.buttons: list[QToolButton] = []
        self.separator: QFrame | None = None

    def _track(self, b: QToolButton) -> None:
        self.buttons.append(b)
        b.defaultAction().changed.connect(self.sync_visibility)
        b.defaultAction().triggered.connect(lambda *_: self.triggered.emit())

    def sync_visibility(self) -> None:
        show = any(b.defaultAction().isVisible() for b in self.buttons)
        self.setVisible(show)
        if self.separator is not None:
            self.separator.setVisible(show)

    def add_large(self, action: QAction) -> QToolButton:
        b = large_button(action); self.body.addWidget(b); self._track(b)
        return b

    def add_small(self, actions: list[QAction]) -> list[QToolButton]:
        made = [small_button(a) for a in actions]
        for b in made:
            self.body.addWidget(b); self._track(b)
        return made


class RibbonPage(QWidget):

    triggered = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("ribbon_page")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._lay = QHBoxLayout(self); self._lay.setContentsMargins(8, PAGE_PADDING, 8, PAGE_PADDING); self._lay.setSpacing(6)
        self._lay.addStretch(1)
        self.groups: list[RibbonGroup] = []

    def add_group(self, title: str) -> RibbonGroup:
        g = RibbonGroup(title)
        if self.groups:
            sep = QFrame(); sep.setObjectName("ribbon_sep"); sep.setFixedWidth(1)
            self._lay.insertWidget(self._lay.count() - 1, sep); g.separator = sep
        self._lay.insertWidget(self._lay.count() - 1, g)
        g.triggered.connect(self.triggered.emit)
        self.groups.append(g)
        return g


class Ribbon(QWidget):

    collapsed_changed = Signal(bool)

    def __init__(self, parent: QWidget | None = None, *, collapsed: bool = True):
        super().__init__(parent)
        self.setObjectName("ribbon")
        self.tabs = QTabBar(); self.tabs.setObjectName("ribbon_tabs")
        self.tabs.setDrawBase(False); self.tabs.setExpanding(False); self.tabs.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.stack = QStackedWidget(); self.stack.setObjectName("ribbon_stack")
        self.stack.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        self.toggle = QToolButton(); self.toggle.setObjectName("ribbon_toggle"); self.toggle.setAutoRaise(True)
        self.toggle.setIconSize(QSize(SMALL_ICON, SMALL_ICON)); self.toggle.clicked.connect(self._on_toggle)
        self.toggle.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.toggle.setFixedSize(HEADER_HEIGHT, HEADER_HEIGHT)
        self.act_toggle = QAction(L("リボンを畳む / 開く", "Collapse / expand the ribbon"), self)
        self.act_toggle.setShortcut("Ctrl+F1"); self.act_toggle.triggered.connect(self._on_toggle)
        self.addAction(self.act_toggle)
        self._header = QHBoxLayout(); self._header.setContentsMargins(8, 2, 8, 1); self._header.setSpacing(6)
        self._lead = QHBoxLayout(); self._lead.setContentsMargins(0, 0, 0, 0)
        self._trail = QHBoxLayout(); self._trail.setContentsMargins(0, 0, 0, 0); self._trail.setSpacing(6)
        self._header.addLayout(self._lead)
        self._header.addWidget(self.tabs, 0, Qt.AlignmentFlag.AlignVCenter)
        self._header.addWidget(self.toggle, 0, Qt.AlignmentFlag.AlignVCenter)
        self._header.addStretch(1); self._header.addLayout(self._trail)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)
        lay.addLayout(self._header)
        pad = QHBoxLayout(); pad.setContentsMargins(8, 0, 8, 4); pad.addWidget(self.stack); lay.addLayout(pad)
        self.pages: list[RibbonPage] = []
        self._collapsed = collapsed
        self._peek = False          # opened from the collapsed state by a tab click: closes again after a command
        self.stack.setVisible(not collapsed)
        self.tabs.tabBarClicked.connect(self._on_tab_clicked)
        self.tabs.currentChanged.connect(self.stack.setCurrentIndex)
        self._update_toggle()

    def add_page(self, title: str) -> RibbonPage:
        page = RibbonPage()
        self.stack.addWidget(page); self.tabs.addTab(title)
        page.triggered.connect(self._on_command)
        self.pages.append(page)
        self._align_leading()
        return page

    def set_leading(self, w: QWidget) -> None:
        self._lead.addWidget(w, 0, Qt.AlignmentFlag.AlignVCenter)
        self._leading = w
        self._align_leading()

    def _align_leading(self) -> None:
        self.tabs.setFixedHeight(HEADER_HEIGHT)
        self.toggle.setFixedSize(HEADER_HEIGHT, HEADER_HEIGHT)
        w = getattr(self, "_leading", None)
        if w is None:
            return
        w.setFixedHeight(HEADER_HEIGHT)
        for b in getattr(w, "buttons", []):
            b.setFixedSize(HEADER_HEIGHT, HEADER_HEIGHT)

    def add_trailing(self, w: QWidget) -> None:
        self._trail.addWidget(w)

    def page_titles(self) -> list[str]:
        return [self.tabs.tabText(i) for i in range(self.tabs.count())]

    def actions_on_page(self, i: int) -> list[QAction]:
        return [b.defaultAction() for g in self.pages[i].groups for b in g.buttons]

    def all_actions(self) -> list[QAction]:
        return [a for i in range(len(self.pages)) for a in self.actions_on_page(i)]

    def is_collapsed(self) -> bool:
        return self._collapsed

    def is_peeking(self) -> bool:
        return self._peek

    def set_collapsed(self, collapsed: bool, *, remember: bool = True) -> None:
        # remember=False opens the ribbon only until the next command (a peek); it is not saved.
        self._peek = not collapsed and not remember
        if collapsed != self._collapsed:
            self._collapsed = collapsed
            self.stack.setVisible(not collapsed)
        self._update_toggle()
        if remember:
            self.collapsed_changed.emit(collapsed)

    def show_page(self, i: int) -> None:
        self.tabs.setCurrentIndex(i); self.set_collapsed(False, remember=False)

    def _on_toggle(self) -> None:
        # While peeking the button means "keep it open"; otherwise it flips the state.
        self.set_collapsed(False if self._peek else not self._collapsed)

    def _on_tab_clicked(self, i: int) -> None:
        if i < 0:
            return
        if self._collapsed:
            self.set_collapsed(False, remember=False)
        elif i == self.tabs.currentIndex():
            self.set_collapsed(True)

    def _on_command(self) -> None:
        if self._peek:
            self.set_collapsed(True, remember=False)

    def _update_toggle(self) -> None:
        if self._collapsed or self._peek:
            self.toggle.setIcon(icons.icon("chevron_down", SMALL_ICON))
            self.toggle.setToolTip(L("リボンを開いたままにする (Ctrl+F1)", "Keep the ribbon open (Ctrl+F1)"))
        else:
            self.toggle.setIcon(icons.icon("chevron_up", SMALL_ICON))
            self.toggle.setToolTip(L("リボンを畳む。見出しだけを残します。見出しを押すと開きます (Ctrl+F1)",
                                     "Collapse the ribbon to its tab names; click a tab to open it again (Ctrl+F1)"))


class ModeBar(QWidget):

    changed = Signal(int)

    def __init__(self, labels: list[tuple[str, str]], parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("modebar")
        lay = QHBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(4)
        self.buttons: list[QToolButton] = []
        for index, (icon_name, text) in enumerate(labels):
            b = QToolButton(); b.setObjectName("mode_button"); b.setCheckable(True); b.setAutoExclusive(True)
            b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
            b.setIcon(icons.icon(icon_name, SMALL_ICON)); b.setIconSize(QSize(SMALL_ICON, SMALL_ICON))
            b.setText(text); b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            b.clicked.connect(lambda _c=False, i=index: self.changed.emit(i))
            lay.addWidget(b); self.buttons.append(b)
        lay.addStretch()
        if self.buttons:
            self.buttons[0].setChecked(True)

    def set_current(self, index: int) -> None:
        if 0 <= index < len(self.buttons) and not self.buttons[index].isChecked():
            self.buttons[index].setChecked(True)

    def current(self) -> int:
        return next((i for i, b in enumerate(self.buttons) if b.isChecked()), 0)

    def set_text(self, index: int, text: str) -> None:
        if 0 <= index < len(self.buttons):
            self.buttons[index].setText(text)

    def add_trailing(self, widget: QWidget) -> None:
        self.layout().addWidget(widget)


class QuickAccessBar(QWidget):

    def __init__(self, actions: list[QAction], parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("qat")
        lay = QHBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(2)
        self.buttons: list[QToolButton] = []
        for a in actions:
            b = QToolButton(); b.setObjectName("qat_button"); b.setAutoRaise(True)
            b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly); b.setIconSize(QSize(SMALL_ICON, SMALL_ICON))
            b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            _bind(b, a); lay.addWidget(b); self.buttons.append(b)


__all__ = ["Ribbon", "RibbonPage", "RibbonGroup", "QuickAccessBar", "ModeBar", "large_button", "small_button"]
