
import os

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QPoint, Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from adit import lang  # noqa: E402


@pytest.fixture()
def dlg():
    app = QApplication.instance() or QApplication([])
    lang.set_language("ja")
    from adit.gui.sketcher import SketchDialog

    d = SketchDialog("")
    d.resize(900, 600); d.show(); d.activateWindow(); app.processEvents()
    yield d
    d.close()


def _pt(c, i):
    q = c._to_screen(c.sketch.atoms[i])
    return QPoint(round(q.x()), round(q.y()))


def _drag(c, a: QPoint, b: QPoint):
    QTest.mousePress(c, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, a)
    QTest.mouseMove(c, QPoint((a.x() + b.x()) // 2, (a.y() + b.y()) // 2))
    QTest.mouseMove(c, b)
    QTest.mouseRelease(c, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, b)


def _click(c, p: QPoint):
    QTest.mouseClick(c, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, p)


def test_drag_with_other_element_keeps_root(dlg):
    c = dlg.canvas
    _click(c, QPoint(350, 300))
    _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(44, -25))  # C-C
    dlg._pick_element("O")
    _drag(c, _pt(c, 1), _pt(c, 1) + QPoint(44, 25))
    assert [a.elem for a in c.sketch.atoms] == ["C", "C", "O"]
    assert c.sketch.to_smiles() == "CCO"
    assert dlg.smiles() == "CCO"


def test_press_alone_does_not_change_element(dlg):
    c = dlg.canvas
    _click(c, QPoint(350, 300)); dlg._pick_element("O")
    QTest.mousePress(c, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, _pt(c, 0))
    assert c.sketch.atoms[0].elem == "C"
    QTest.mouseRelease(c, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, _pt(c, 0))
    assert c.sketch.atoms[0].elem == "O"


def test_click_changes_element_and_undo_redo(dlg):
    c = dlg.canvas
    _click(c, QPoint(350, 300))
    _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(44, 0))
    dlg._pick_element("N")
    _click(c, _pt(c, 1))
    assert c.sketch.to_smiles() == "CN"
    assert c.undo() and c.sketch.to_smiles() == "CC"
    assert c.redo() and c.sketch.to_smiles() == "CN"


def test_clear_all_is_undoable_with_keys(dlg):
    c = dlg.canvas
    assert dlg.btn_clear.text() == "すべて消去"
    _click(c, QPoint(350, 300)); _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(44, 0))
    assert c.sketch.to_smiles() == "CC"
    dlg.btn_clear.click()
    assert c.sketch.atoms == [] and dlg.btn_undo.isEnabled()
    QTest.keyClick(dlg, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
    assert c.sketch.to_smiles() == "CC"
    QTest.keyClick(dlg, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier)
    assert c.sketch.atoms == []
    QTest.keyClick(dlg, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
    assert c.sketch.to_smiles() == "CC"
    QTest.keyClick(dlg, Qt.Key.Key_Y, Qt.KeyboardModifier.ControlModifier)
    assert c.sketch.atoms == []


def test_history_starts_at_opened_molecule():
    QApplication.instance() or QApplication([])
    from adit.gui.sketcher import SketchDialog

    d = SketchDialog("CCO")
    assert not d.canvas.can_undo() and not d.btn_undo.isEnabled()
    assert d.smiles() == "CCO"


def test_delete_key_and_right_click_remove_one(dlg):
    c = dlg.canvas
    _click(c, QPoint(350, 300)); _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(44, 0))
    _drag(c, _pt(c, 1), _pt(c, 1) + QPoint(44, 0))
    assert c.sketch.to_smiles() == "CCC"
    QTest.mouseClick(c, Qt.MouseButton.RightButton, Qt.KeyboardModifier.NoModifier, _pt(c, 2))
    assert c.sketch.to_smiles() == "CC"
    _click(c, _pt(c, 1)); c.setFocus()
    QTest.keyClick(c, Qt.Key.Key_Delete)
    assert c.sketch.to_smiles() == "C"
    assert c.undo() and c.sketch.to_smiles() == "CC"


def test_charge_without_selection_tells_what_to_do(dlg):
    dlg.canvas.selected = None
    dlg._charge(1)
    assert "原子をクリック" in dlg.note.text() and dlg.note.isVisibleTo(dlg)


def test_smiles_is_centered():
    from adit.gui.sketcher import Sketch

    sk = Sketch.from_smiles("CC(=O)[O-].c1ccncc1")
    xs, ys = [a.x for a in sk.atoms], [a.y for a in sk.atoms]
    assert abs(min(xs) + max(xs)) < 1e-6 and abs(min(ys) + max(ys)) < 1e-6


def _mid(c, bi):
    b = c.sketch.bonds[bi]
    p, q = _pt(c, b.a), _pt(c, b.b)
    return QPoint((p.x() + q.x()) // 2, (p.y() + q.y()) // 2)


@pytest.mark.parametrize("language", ["ja", "en"])
def test_window_title_is_ascii(language):
    QApplication.instance() or QApplication([])
    lang.set_language(language)
    try:
        from adit.gui.sketcher import SketchDialog

        t = SketchDialog("").windowTitle()
        assert t == "ADIT Draw" and t.isascii()
    finally:
        lang.set_language("ja")


def test_ring_buttons_show_formulas(dlg):
    assert [b.text() for b in dlg._tmpl_buttons] == ["C₆H₆", "C₆H₁₂", "C₅H₁₀"]
    assert [b.text() for b in dlg._bond_buttons.values()] == ["—", "=", "≡"]


def test_eraser_removes_atom_with_bonds_and_undo(dlg):
    c = dlg.canvas
    _click(c, QPoint(350, 300)); _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(44, 0)); _drag(c, _pt(c, 1), _pt(c, 1) + QPoint(44, 0))
    assert c.sketch.to_smiles() == "CCC"
    _click(dlg.btn_eraser, QPoint(5, 5))
    assert dlg.btn_eraser.isChecked() and not dlg._bond_buttons[1].isChecked() and c.tool == "eraser"
    _click(c, _pt(c, 1))
    assert len(c.sketch.atoms) == 2 and c.sketch.bonds == [] and c.sketch.to_smiles() == "C.C"
    QTest.keyClick(dlg, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
    assert c.sketch.to_smiles() == "CCC"
    _click(c, _mid(c, 0))
    assert len(c.sketch.atoms) == 3 and len(c.sketch.bonds) == 1
    _click(c, QPoint(700, 100))
    assert len(c.sketch.atoms) == 3
    dlg._pick_element("N")
    assert c.tool == "bond" and dlg._bond_buttons[1].isChecked() and not dlg.btn_eraser.isChecked()


def test_bond_tool_joins_ring_and_atom_to_toluene(dlg):
    c = dlg.canvas
    _click(dlg._tmpl_buttons[0], QPoint(5, 5)); _click(c, QPoint(300, 300))
    _click(dlg._bond_buttons[1], QPoint(5, 5)); _click(c, QPoint(650, 300))
    assert c.sketch.to_smiles() == "C.c1ccccc1"
    _drag(c, _pt(c, 1), _pt(c, 6))
    assert c.sketch.to_smiles() == "Cc1ccccc1" and dlg.smiles() == "Cc1ccccc1"
    assert c.undo() and c.sketch.to_smiles() == "C.c1ccccc1"


def test_separate_molecules_are_clearly_reported(dlg):
    c = dlg.canvas
    _click(dlg._tmpl_buttons[0], QPoint(5, 5)); _click(c, QPoint(300, 300)); _click(c, QPoint(390, 300))
    assert c.sketch.component_count() == 2 and "." in dlg.smiles()
    assert "分子が 2 つあります" in dlg.smiles_label.text()
    assert dlg.smiles_label.objectName() == "status_ng"


def test_double_bond_tool_between_existing_atoms_and_new_atom(dlg):
    c = dlg.canvas
    _click(c, QPoint(300, 300)); _click(c, QPoint(500, 300))
    _click(dlg._bond_buttons[2], QPoint(5, 5))
    _drag(c, _pt(c, 0), _pt(c, 1))
    assert c.sketch.to_smiles() == "C=C"
    dlg._pick_element("O"); assert c.tool == "bond" and c.bond_order == 2
    _drag(c, _pt(c, 1), _pt(c, 1) + QPoint(0, -44))
    assert c.sketch.to_smiles() == "C=C=O"
    _click(dlg._bond_buttons[1], QPoint(5, 5))
    _drag(c, _pt(c, 0), _pt(c, 1))
    assert c.sketch.to_smiles() == "CC=O" and len(c.sketch.bonds) == 2


def test_bond_click_sets_selected_order(dlg):
    c = dlg.canvas
    _click(c, QPoint(350, 300)); _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(88, 0))
    _click(dlg._bond_buttons[3], QPoint(5, 5))
    _click(c, _mid(c, 0)); assert c.sketch.to_smiles() == "C#C"
    _click(dlg._bond_buttons[2], QPoint(5, 5))
    _click(c, _mid(c, 0)); assert c.sketch.to_smiles() == "C=C"
    _click(c, _mid(c, 0)); assert c.sketch.to_smiles() == "C#C"
    assert c.undo() and c.sketch.to_smiles() == "C=C"


def test_overvalent_bond_is_reported_not_blocked(dlg):
    c = dlg.canvas
    dlg._pick_element("O"); _click(c, QPoint(350, 300))
    _click(dlg._bond_buttons[3], QPoint(5, 5)); _drag(c, _pt(c, 0), _pt(c, 0) + QPoint(88, 0))
    assert c.sketch.bonds[0].order == 3 and dlg.smiles() == ""
    assert "分子になりません" in dlg.smiles_label.text()


def test_dark_theme_colors_follow_tokens():
    QApplication.instance() or QApplication([])
    from adit.gui import style
    from adit.gui.sketcher import SketchDialog, _palette

    saved_theme = style._APPLIED
    try:
        style._APPLIED = "dark"
        d = SketchDialog("")
        assert style.DARK.fg in d.styleSheet() and style.DARK.bg in d.styleSheet()
        assert _palette()["face"].name().upper() == style.DARK.card.upper()
        style._APPLIED = "light"
        assert _palette()["face"].name().upper() == style.LIGHT.card.upper()
    finally:
        style._APPLIED = saved_theme
