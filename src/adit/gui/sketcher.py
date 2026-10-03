
from __future__ import annotations

import copy
import math
import re

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QKeyEvent, QKeySequence, QMouseEvent, QPainter, QPen, QShortcut
from PySide6.QtWidgets import (QButtonGroup, QDialog, QDialogButtonBox, QHBoxLayout, QInputDialog, QLabel, QMessageBox, QPushButton, QToolButton,
                               QVBoxLayout, QWidget)

from adit.gui.style import RADIUS_SMALL, current_theme, tokens_for
from adit.sketch import (BOND_PX, BOND_SYMBOLS, CLICK_PX, ELEMENT_COLORS, ELEMENT_COLORS_DARK, ELEMENTS,
                          HISTORY_MAX, SAtom, SBond, Sketch, TEMPLATE_FORMULAS, TEMPLATES, WINDOW_TITLE)
from adit.lang import L
from adit.structure import pretty_formula

def _qcolor(s: str) -> QColor:
    m = re.fullmatch(r"\s*rgba\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*", s)
    if m:
        return QColor(*(int(v) for v in m.groups()))
    return QColor(s)


def _palette() -> dict:
    theme = current_theme()
    t = tokens_for(theme)
    dark = theme == "dark"
    return {"face": _qcolor(t.card), "edge": _qcolor(t.muted if dark else "#B8C0D0"), "bond": _qcolor(t.fg), "muted": _qcolor(t.muted),
            "accent": _qcolor(t.accent), "ng": _qcolor(t.ng), "colors": {**ELEMENT_COLORS, **(ELEMENT_COLORS_DARK if dark else {}), "C": t.fg}}


class SketchCanvas(QWidget):
    changed = Signal()
    history_changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setMinimumSize(520, 360)
        self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.sketch = Sketch()
        self.element = "C"
        self.tool = "bond"  # "bond" / "ring" / "eraser"
        self.bond_order = 1
        self.template: str | None = None
        self.selected: int | None = None
        self._hover: int | SBond | None = None
        self._drag_from: int | None = None
        self._drag_to: QPointF | None = None
        self._offset = QPointF(0, 0)
        self._undo: list[Sketch] = []
        self._redo: list[Sketch] = []
        self.setMouseTracking(True)

    def _record(self) -> None:
        self._undo.append(copy.deepcopy(self.sketch))
        del self._undo[:-HISTORY_MAX]
        self._redo.clear()
        self.history_changed.emit()

    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo(self) -> bool:
        if not self._undo:
            return False
        self._redo.append(self.sketch); self.sketch = self._undo.pop()
        self._after_history(); return True

    def redo(self) -> bool:
        if not self._redo:
            return False
        self._undo.append(self.sketch); self.sketch = self._redo.pop()
        self._after_history(); return True

    def reset_history(self) -> None:
        self._undo.clear(); self._redo.clear(); self.history_changed.emit()

    def _after_history(self) -> None:
        self.selected = None; self._drag_from = None; self._drag_to = None; self._hover = None
        self.update(); self.changed.emit(); self.history_changed.emit()

    def _to_screen(self, a: SAtom) -> QPointF:
        return QPointF(a.x + self.width() / 2 + self._offset.x(), a.y + self.height() / 2 + self._offset.y())

    def _to_model(self, p: QPointF) -> tuple[float, float]:
        return p.x() - self.width() / 2 - self._offset.x(), p.y() - self.height() / 2 - self._offset.y()

    def _atom_at(self, p: QPointF) -> int | None:
        for i, a in enumerate(self.sketch.atoms):
            q = self._to_screen(a)
            if (q.x() - p.x()) ** 2 + (q.y() - p.y()) ** 2 <= 12 ** 2:
                return i
        return None

    def _bond_at(self, p: QPointF) -> SBond | None:
        for b in self.sketch.bonds:
            pa, pb = self._to_screen(self.sketch.atoms[b.a]), self._to_screen(self.sketch.atoms[b.b])
            dx, dy = pb.x() - pa.x(), pb.y() - pa.y(); ln = math.hypot(dx, dy) or 1
            t = ((p.x() - pa.x()) * dx + (p.y() - pa.y()) * dy) / (ln * ln)
            if 0.1 < t < 0.9:
                d = abs((p.x() - pa.x()) * dy - (p.y() - pa.y()) * dx) / ln
                if d <= 6:
                    return b
        return None

    def add_template(self, name: str, x: float, y: float) -> None:
        n, aromatic = TEMPLATES[name]
        r = BOND_PX / (2 * math.sin(math.pi / n))
        start = len(self.sketch.atoms)
        for k in range(n):
            t = 2 * math.pi * k / n - math.pi / 2
            self.sketch.atoms.append(SAtom(x + r * math.cos(t), y + r * math.sin(t), "C"))
        for k in range(n):
            order = 2 if aromatic and k % 2 == 0 else 1
            self.sketch.bonds.append(SBond(start + k, start + (k + 1) % n, order))

    def load_smiles(self, smiles: str) -> None:
        new = Sketch.from_smiles(smiles)
        self._record()
        self.sketch = new; self.selected = None; self._offset = QPointF(0, 0); self.update(); self.changed.emit()

    def clear(self) -> None:
        if not self.sketch.atoms:
            return
        self._record()
        self.sketch = Sketch(); self.selected = None; self.update(); self.changed.emit()

    def set_charge(self, delta: int) -> bool:
        if self.selected is None:
            return False
        self._record()
        self.sketch.atoms[self.selected].charge += delta; self.update(); self.changed.emit()
        return True

    def delete_selected(self) -> bool:
        if self.selected is None:
            return False
        self._record()
        self.sketch.remove_atom(self.selected); self.selected = None; self.update(); self.changed.emit()
        return True

    def set_tool(self, tool: str, *, order: int | None = None, template: str | None = None) -> None:
        self.tool = tool
        if order is not None:
            self.bond_order = order
        self.template = template if tool == "ring" else None
        if tool == "eraser":
            self.selected = None
        self._hover = None
        self.setCursor(Qt.CursorShape.PointingHandCursor if tool == "eraser" else Qt.CursorShape.ArrowCursor)
        self.update()

    def erase_at(self, p: QPointF) -> bool:
        i = self._atom_at(p)
        b = None if i is not None else self._bond_at(p)
        if i is None and b is None:
            return False
        self._record()
        if i is not None:
            self.sketch.remove_atom(i)
        else:
            self.sketch.bonds.remove(b)
        self.selected = None; self._hover = None; self.update(); self.changed.emit()
        return True

    def mousePressEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        p = ev.position(); i = self._atom_at(p)
        if ev.button() == Qt.MouseButton.RightButton:
            self.erase_at(p); return
        if ev.button() != Qt.MouseButton.LeftButton:
            return
        if self.tool == "eraser":
            self.erase_at(p); return
        if i is not None:
            self.selected = i; self._drag_from = i; self._drag_to = p
            self.update(); return
        b = self._bond_at(p)
        if b is not None:
            self._record(); b.order = self.bond_order if b.order != self.bond_order else b.order % 3 + 1
            self.update(); self.changed.emit(); return
        x, y = self._to_model(p)
        self._record()
        if self.tool == "ring" and self.template:
            self.add_template(self.template, x, y)
        else:
            self.sketch.atoms.append(SAtom(x, y, self.element)); self.selected = len(self.sketch.atoms) - 1
        self.update(); self.changed.emit()

    def _hover_target(self, p: QPointF) -> int | SBond | None:
        i = self._atom_at(p)
        return i if i is not None else self._bond_at(p)

    def mouseMoveEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        if self._drag_from is not None:
            self._drag_to = ev.position(); self.update()
        elif self.tool == "eraser":
            h = self._hover_target(ev.position())
            if h is not self._hover and h != self._hover:
                self._hover = h; self.update()

    def leaveEvent(self, ev) -> None:  # noqa: N802
        if self._hover is not None:
            self._hover = None; self.update()
        super().leaveEvent(ev)

    def mouseReleaseEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        if self._drag_from is None:
            return
        i, p = self._drag_from, ev.position(); self._drag_from = None; self._drag_to = None
        src = self._to_screen(self.sketch.atoms[i])
        if math.hypot(p.x() - src.x(), p.y() - src.y()) < CLICK_PX:
            if self.tool == "bond" and self.sketch.atoms[i].elem != self.element:
                self._record(); self.sketch.atoms[i].elem = self.element; self.changed.emit()
            self.update(); return
        j = self._atom_at(p)
        if j == i:
            self.update(); return
        if j is not None:
            old = self.sketch.bond_between(i, j)
            if old is not None:
                if old.order != self.bond_order:
                    self._record(); old.order = self.bond_order; self.changed.emit()
                self.selected = j; self.update(); return
        self._record()
        if j is None:
            ang = math.atan2(p.y() - src.y(), p.x() - src.x())
            ang = math.radians(round(math.degrees(ang) / 30) * 30)
            a = self.sketch.atoms[i]
            self.sketch.atoms.append(SAtom(a.x + BOND_PX * math.cos(ang), a.y + BOND_PX * math.sin(ang), self.element))
            j = len(self.sketch.atoms) - 1
        self.sketch.bonds.append(SBond(i, j, self.bond_order))
        self.selected = j; self.update(); self.changed.emit()

    def keyPressEvent(self, ev: QKeyEvent) -> None:  # noqa: N802
        if ev.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace) and self.delete_selected():
            return
        super().keyPressEvent(ev)

    def paintEvent(self, _ev) -> None:  # noqa: N802
        pal = _palette()
        p = QPainter(self); p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(QPen(pal["edge"], 1.0)); p.setBrush(pal["face"])
        p.drawRoundedRect(QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5), RADIUS_SMALL, RADIUS_SMALL)
        sk = self.sketch
        if not sk.atoms:
            p.setPen(pal["muted"])
            p.drawText(QRectF(self.rect()), Qt.AlignmentFlag.AlignCenter,
                       L("ここをクリックすると、選んでいる元素の原子を置きます\n原子からドラッグすると結合が伸びます",
                         "Click here to place an atom of the selected element\nDrag from an atom to draw a bond"))
        for b in sk.bonds:
            a1, a2 = sk.atoms[b.a], sk.atoms[b.b]; p1, p2 = self._to_screen(a1), self._to_screen(a2)
            self._draw_bond(p, p1, p2, b.order, a1.elem != "C" or a1.charge != 0, a2.elem != "C" or a2.charge != 0, pal["bond"])
        if self._drag_from is not None and self._drag_to is not None:
            pen = QPen(pal["muted"]); pen.setStyle(Qt.PenStyle.DashLine); pen.setWidthF(1.5); p.setPen(pen)
            p.drawLine(self._to_screen(sk.atoms[self._drag_from]), self._drag_to)
        if self.tool == "eraser" and isinstance(self._hover, SBond) and self._hover in sk.bonds:
            hb = self._hover; pen = QPen(pal["ng"]); pen.setWidthF(9.0); pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            c = QColor(pal["ng"]); c.setAlpha(70); pen.setColor(c); p.setPen(pen)
            p.drawLine(self._to_screen(sk.atoms[hb.a]), self._to_screen(sk.atoms[hb.b]))
        font = QFont(self.font()); font.setPointSize(12); font.setBold(True); p.setFont(font)
        hs = sk.h_counts()
        sel_fill = QColor(pal["accent"]); sel_fill.setAlpha(40)
        ng_fill = QColor(pal["ng"]); ng_fill.setAlpha(50)
        for i, a in enumerate(sk.atoms):
            q = self._to_screen(a)
            if i == self.selected:
                p.setPen(QPen(pal["accent"], 1.5)); p.setBrush(sel_fill); p.drawEllipse(q, 13, 13)
            if self.tool == "eraser" and isinstance(self._hover, int) and not isinstance(self._hover, bool) and i == self._hover:
                p.setPen(QPen(pal["ng"], 1.5)); p.setBrush(ng_fill); p.drawEllipse(q, 13, 13)
            if a.elem != "C" or a.charge != 0 or sk.degree(i) == 0:
                h = hs[i] if i < len(hs) else 0
                text = a.elem + (("H" + (str(h) if h > 1 else "")) if h else "")
                text += (("+" if a.charge > 0 else "−") + (str(abs(a.charge)) if abs(a.charge) > 1 else "")) if a.charge else ""
                w = 14 + 9 * len(text)
                p.setPen(Qt.PenStyle.NoPen); p.setBrush(pal["face"]); p.drawEllipse(q, w / 2, 11)
                p.setPen(QColor(pal["colors"].get(a.elem, pal["colors"]["C"]))); p.drawText(QRectF(q.x() - w / 2, q.y() - 12, w, 24), Qt.AlignmentFlag.AlignCenter, text)
        p.end()

    def _draw_bond(self, p: QPainter, p1: QPointF, p2: QPointF, order: int, trim1: bool, trim2: bool, color: QColor | None = None) -> None:
        dx, dy = p2.x() - p1.x(), p2.y() - p1.y(); ln = math.hypot(dx, dy) or 1
        ux, uy = dx / ln, dy / ln
        if trim1:
            p1 = QPointF(p1.x() + ux * 9, p1.y() + uy * 9)
        if trim2:
            p2 = QPointF(p2.x() - ux * 9, p2.y() - uy * 9)
        pen = QPen(color or QColor("#202020")); pen.setWidthF(2.0); pen.setCapStyle(Qt.PenCapStyle.RoundCap); p.setPen(pen)
        nx, ny = -uy * 4, ux * 4
        offsets = {1: [0], 2: [-1, 1], 3: [-1, 0, 1]}[order]
        for k in offsets:
            p.drawLine(QPointF(p1.x() + nx * k, p1.y() + ny * k), QPointF(p2.x() + nx * k, p2.y() + ny * k))


def _dialog_qss() -> str:
    t = tokens_for(current_theme()); rs = RADIUS_SMALL
    return f"""
QDialog {{ background: {t.bg}; color: {t.fg}; }}
QLabel {{ color: {t.fg}; background: transparent; }}
QLabel#hint {{ color: {t.muted}; font-size: 10pt; }}
QLabel#status_ng {{ color: {t.ng}; font-weight: 600; }}
QLabel#smiles {{ color: {t.fg}; font-weight: 600; }}
QLabel#section {{ color: {t.muted}; font-weight: 700; font-size: 10pt; }}
QToolButton#bond {{ font-size: 13pt; font-weight: 700; min-width: 22px; padding: 2px 9px; }}
QToolButton {{ color: {t.fg}; background: {t.field}; border: 1px solid {t.line}; border-radius: {rs}px; padding: 4px 9px; min-height: 22px; }}
QToolButton:hover {{ border: 1px solid {t.accent}; }}
QToolButton:checked {{ background: {t.accent}; color: #FFFFFF; border: 1px solid {t.accent}; font-weight: 600; }}
QToolButton:disabled, QPushButton:disabled {{ color: {t.muted}; background: transparent; border: 1px solid {t.line}; }}
"""


class _WrapLabel(QLabel):

    def resizeEvent(self, ev) -> None:  # noqa: N802
        super().resizeEvent(ev)
        h = self.heightForWidth(self.width())
        if h > 0 and h != self.minimumHeight():
            self.setMinimumHeight(h)


class SketchDialog(QDialog):

    def __init__(self, smiles: str = "", parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle(WINDOW_TITLE); self.resize(900, 620)
        self.setStyleSheet(_dialog_qss())
        self.canvas = SketchCanvas()
        self.smiles_label = QLabel(""); self.smiles_label.setObjectName("smiles"); self.smiles_label.setWordWrap(True)
        self.smiles_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self._elem_group = QButtonGroup(self); self._elem_group.setExclusive(True)
        self._elem_buttons: list[QToolButton] = []
        for e in ELEMENTS:
            b = QToolButton(); b.setText(e); b.setCheckable(True); b.setChecked(e == "C"); b.clicked.connect(lambda _c, e=e: self._pick_element(e))
            b.setToolTip(L(f"{e} を置く。既存の原子をクリックすると {e} に置き換わります", f"Place {e}. Clicking an existing atom changes it to {e}"))
            self._elem_group.addButton(b); self._elem_buttons.append(b)

        self._tool_group = QButtonGroup(self); self._tool_group.setExclusive(True)
        self._bond_buttons: dict[int, QToolButton] = {}
        bond_names = {1: (L("単結合", "Single bond"), "single"), 2: (L("二重結合", "Double bond"), "double"), 3: (L("三重結合", "Triple bond"), "triple")}
        for order, (name, en) in bond_names.items():
            b = QToolButton(); b.setObjectName("bond"); b.setText(BOND_SYMBOLS[order]); b.setCheckable(True); b.setChecked(order == 1)
            b.setAccessibleName(name); b.clicked.connect(lambda _c, o=order: self._pick_bond(o))
            b.setToolTip(L(f"{name}。原子からドラッグするとこの結合で原子が伸び、別の原子の上で離すとその 2 つをつなぎます。"
                           f"結合をクリックすると{name}に変わります",
                           f"{name}. Drag from an atom to add a new atom with a {en} bond, or release on another atom to join the two. "
                           f"Click a bond to make it {en}"))
            self._tool_group.addButton(b); self._bond_buttons[order] = b
        self._tmpl_buttons: list[QToolButton] = []
        for name in TEMPLATES:
            f = pretty_formula(TEMPLATE_FORMULAS[name])
            b = QToolButton(); b.setText(f); b.setCheckable(True); b.clicked.connect(lambda _c, n=name: self._pick_template(n))
            b.setToolTip(L(f"{f} の環。選んでから空いた所をクリックすると置きます", f"{f} ring. Select it, then click an empty spot to place it"))
            self._tool_group.addButton(b); self._tmpl_buttons.append(b)
        self.btn_eraser = QToolButton(); self.btn_eraser.setText(L("消しゴム", "Eraser")); self.btn_eraser.setCheckable(True)
        self.btn_eraser.setToolTip(L("選んでいる間は、クリックした原子か結合を消します (原子を消すと、つながった結合も消えます)。Ctrl+Z で戻せます",
                                     "While selected, clicking an atom or bond removes it (removing an atom also removes its bonds). Ctrl+Z brings it back"))
        self.btn_eraser.clicked.connect(self._pick_eraser); self._tool_group.addButton(self.btn_eraser)

        plus = QPushButton("+"); minus = QPushButton("−")
        plus.setToolTip(L("選んだ原子 (青い丸) の電荷を 1 増やします", "Raise the charge of the selected atom (blue circle) by 1"))
        minus.setToolTip(L("選んだ原子 (青い丸) の電荷を 1 減らします", "Lower the charge of the selected atom (blue circle) by 1"))
        plus.clicked.connect(lambda: self._charge(1)); minus.clicked.connect(lambda: self._charge(-1))
        self.btn_undo = QPushButton(L("元に戻す", "Undo")); self.btn_undo.setToolTip("Ctrl+Z"); self.btn_undo.clicked.connect(self.canvas.undo)
        self.btn_redo = QPushButton(L("やり直す", "Redo")); self.btn_redo.setToolTip("Ctrl+Shift+Z / Ctrl+Y"); self.btn_redo.clicked.connect(self.canvas.redo)
        load = QPushButton(L("SMILES から…", "From SMILES…")); load.clicked.connect(self._load)
        load.setToolTip(L("SMILES を入力して 2 次元に展開し、続けて編集します", "Enter a SMILES, lay it out in 2D and keep editing"))
        self.btn_clear = QPushButton(L("すべて消去", "Clear all")); self.btn_clear.clicked.connect(self.canvas.clear)
        self.btn_clear.setToolTip(L("描いたものをすべて消します。Ctrl+Z で戻せます", "Erase everything. Ctrl+Z brings it back"))

        def section(title: str, widgets) -> QVBoxLayout:
            box = QVBoxLayout(); box.setSpacing(4)
            head = QLabel(title); head.setObjectName("section"); box.addWidget(head)
            row = QHBoxLayout(); row.setSpacing(4)
            for w in widgets:
                row.addWidget(w)
            box.addLayout(row)
            return box

        row1 = QHBoxLayout(); row1.setSpacing(18)
        row1.addLayout(section(L("元素", "Element"), self._elem_buttons))
        row1.addLayout(section(L("結合", "Bond"), self._bond_buttons.values()))
        row1.addLayout(section(L("環", "Ring"), self._tmpl_buttons))
        row1.addStretch()
        row2 = QHBoxLayout(); row2.setSpacing(18)
        row2.addLayout(section(L("ツール", "Tool"), [self.btn_eraser]))
        row2.addLayout(section(L("編集", "Edit"), [self.btn_undo, self.btn_redo, self.btn_clear]))
        row2.addLayout(section(L("電荷", "Charge"), [plus, minus]))
        row2.addStretch()
        row2.addLayout(section(L("読み込み", "Import"), [load]))

        for seq in ("Ctrl+Z",):
            QShortcut(QKeySequence(seq), self, activated=self.canvas.undo)
        for seq in ("Ctrl+Shift+Z", "Ctrl+Y"):
            QShortcut(QKeySequence(seq), self, activated=self.canvas.redo)

        hint = _WrapLabel(L("空いた所をクリック: 選んだ元素の原子を置く ・ 原子をクリック: 選んだ元素に置き換える ・ 水素は自動で補われます\n"
                        "原子からドラッグ: 選んだ結合で新しい原子をつなぐ (別の原子の上で離すと、その 2 つをつなぐ) ・ "
                        "結合をクリック: 選んだ結合に変える (同じ結合なら 単 → 二重 → 三重 と進む)\n"
                        "消しゴム: クリックした原子・結合を消す (右クリックや、原子を選んで Delete でも消せます) ・ Ctrl+Z: 元に戻す",
                        "Click empty space: place an atom of the selected element · Click an atom: change it to the selected element · "
                        "Hydrogens are added automatically\n"
                        "Drag from an atom: attach a new atom with the selected bond (release on another atom to join the two) · "
                        "Click a bond: change it to the selected bond (if it already is, it cycles single → double → triple)\n"
                        "Eraser: click an atom or bond to remove it (right-click, or select an atom and press Delete, also works) · Ctrl+Z: undo"))
        hint.setObjectName("hint"); hint.setWordWrap(True); self._hint = hint
        self.note = QLabel(""); self.note.setObjectName("hint")
        from adit.gui.copy_save import copy_button, save_button

        self.btn_copy = copy_button(self)
        self.btn_save = save_button(self)
        self.btn_copy.clicked.connect(self._copy_image)
        self.btn_save.clicked.connect(self._save_image)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText(L("この分子を使う", "Use this molecule"))
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText(L("キャンセル", "Cancel"))
        buttons.accepted.connect(self._accept); buttons.rejected.connect(self.reject)
        buttons.addButton(self.btn_copy, QDialogButtonBox.ButtonRole.ActionRole)
        buttons.addButton(self.btn_save, QDialogButtonBox.ButtonRole.ActionRole)
        lay = QVBoxLayout(self); lay.setContentsMargins(14, 14, 14, 14); lay.setSpacing(10)
        lay.addLayout(row1); lay.addLayout(row2); lay.addWidget(self.canvas, 1); lay.addWidget(hint); lay.addWidget(self.note)
        lay.addWidget(self.smiles_label); lay.addWidget(buttons)
        self.canvas.changed.connect(self._refresh)
        self.canvas.history_changed.connect(self._refresh_history)
        if smiles:
            try:
                self.canvas.load_smiles(smiles)
            except ValueError:
                pass
        self.canvas.reset_history()
        self._refresh()

    def _pick_element(self, e: str) -> None:
        self.canvas.element = e
        for b in self._elem_buttons:
            b.setChecked(b.text() == e)
        if self.canvas.tool != "bond":
            self._pick_bond(self.canvas.bond_order)

    def _pick_bond(self, order: int) -> None:
        self.canvas.set_tool("bond", order=order)
        self._bond_buttons[order].setChecked(True)

    def _pick_template(self, name: str) -> None:
        self.canvas.set_tool("ring", template=name)
        for b, n in zip(self._tmpl_buttons, TEMPLATES):
            if n == name:
                b.setChecked(True)

    def _pick_eraser(self) -> None:
        self.canvas.set_tool("eraser")
        self.btn_eraser.setChecked(True)

    def _charge(self, delta: int) -> None:
        if not self.canvas.set_charge(delta):
            self.note.setText(L("電荷を変えるには、先に原子をクリックして選んでください (青い丸が付きます)",
                                "Click an atom first to select it (a blue circle appears), then change its charge"))
            self.note.setVisible(True)

    def _load(self) -> None:
        text, ok = QInputDialog.getText(self, "SMILES", "SMILES:")
        if ok and text.strip():
            try:
                self.canvas.load_smiles(text.strip())
            except ValueError as ex:
                QMessageBox.warning(self, "SMILES", str(ex))

    def _refresh_history(self) -> None:
        self.btn_undo.setEnabled(self.canvas.can_undo()); self.btn_redo.setEnabled(self.canvas.can_redo())

    def _refresh(self) -> None:
        self.note.setText(""); self.note.setVisible(False)
        try:
            self._smiles = self.canvas.sketch.to_smiles()
            if self._smiles:
                components = self.canvas.sketch.component_count()
                separated = (L(f"   分子が {components} つあります (互いにつながっていません。SMILES の . は分離を表します)",
                               f"   There are {components} separate molecules (they are not connected; a dot in SMILES denotes separation)")
                             if components > 1 else "")
                self.smiles_label.setText(L(f"SMILES: {self._smiles}   分子式: {self.canvas.sketch.formula()}",
                                            f"SMILES: {self._smiles}   Formula: {self.canvas.sketch.formula()}") + separated)
            else:
                self.smiles_label.setText(L("(まだ何も描いていません)", "(nothing drawn yet)"))
            self.smiles_label.setObjectName("status_ng" if self.canvas.sketch.component_count() > 1 else "smiles")
        except Exception as ex:
            self._smiles = ""
            self.smiles_label.setText(L(f"この形は分子になりません (Ctrl+Z で 1 つ前に戻せます): {ex}", f"not a valid molecule (Ctrl+Z to step back): {ex}"))
            self.smiles_label.setObjectName("status_ng")
        self.smiles_label.style().unpolish(self.smiles_label); self.smiles_label.style().polish(self.smiles_label)
        self._refresh_history()

    def _accept(self) -> None:
        if not self._smiles:
            QMessageBox.warning(self, WINDOW_TITLE, L("分子として成り立つ形にしてから「この分子を使う」を押してください", "Make a valid molecule before pressing \"Use this molecule\"")); return
        self.accept()

    def smiles(self) -> str:
        return self._smiles

    def _copy_image(self) -> None:
        from adit.gui.copy_save import copy_widget

        copy_widget(self.canvas)

    def _save_image(self) -> None:
        from adit.export_image import sketch_to_mol, sketch_to_svg
        from adit.gui.copy_save import save_dialog

        sketch = self.canvas.sketch
        if not sketch.atoms:
            return
        try:
            svg = sketch_to_svg(sketch)
        except ValueError:
            svg = None
        try:
            mol = sketch_to_mol(sketch)
        except Exception:
            mol = None
        save_dialog(self, "structure.svg", svg_text=svg, widget=self.canvas, mol_text=mol)
