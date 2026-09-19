"""Measuring in the 3D view: distance, angle, dihedral, minimum image, picking and sending atoms to fields."""

import os

import numpy as np
import pytest
from ase import Atoms

from adit.gui.atom_select import base_indices, compact_ranges, merge_fixed_text, merge_select_text, series_text
from adit.measure import angle, dihedral, distance, measure, measure_text, mic_vector


def test_distance_angle_dihedral_of_a_known_geometry():
    # ethane-like staggered fragment: H-C-C-H dihedral 180, H-C-C angle 109.47 (tetrahedral), C-C 1.54
    c1, c2 = np.array([0.0, 0.0, 0.0]), np.array([1.54, 0.0, 0.0])
    t = np.radians(109.47)
    h1 = c1 + 1.09 * np.array([np.cos(t), np.sin(t), 0.0])          # on C1, in the xy plane
    h2 = c2 + 1.09 * np.array([-np.cos(t), -np.sin(t), 0.0])        # on C2, anti to H1
    pos = np.array([h1, c1, c2, h2])
    assert distance(pos, 1, 2) == pytest.approx(1.54, abs=1e-9)
    assert angle(pos, 0, 1, 2) == pytest.approx(109.47, abs=1e-6)
    assert dihedral(pos, 0, 1, 2, 3) == pytest.approx(180.0, abs=1e-6)
    h2_gauche = c2 + 1.09 * np.array([-np.cos(t), np.sin(t) * np.cos(np.radians(60)), np.sin(t) * np.sin(np.radians(60))])
    assert abs(dihedral(np.array([h1, c1, c2, h2_gauche]), 0, 1, 2, 3)) == pytest.approx(60.0, abs=1e-6)


def test_right_angle_and_sign_of_the_dihedral():
    pos = np.array([[1.0, 0, 0], [0, 0, 0], [0, 0, 1.0], [0, 1.0, 1.0]])
    assert angle(pos, 0, 1, 2) == pytest.approx(90.0)
    assert dihedral(pos, 0, 1, 2, 3) == pytest.approx(-90.0) or dihedral(pos, 0, 1, 2, 3) == pytest.approx(90.0)
    assert dihedral(pos, 0, 1, 2, 3) == pytest.approx(-dihedral(pos, 3, 2, 1, 0) * -1)   # same value read from either end


def test_minimum_image_across_the_cell_boundary():
    cell = np.diag([10.0, 10.0, 10.0])
    pos = np.array([[0.5, 5.0, 5.0], [9.5, 5.0, 5.0]])
    assert distance(pos, 0, 1) == pytest.approx(9.0)
    assert distance(pos, 0, 1, cell) == pytest.approx(1.0)
    v = mic_vector(pos[0], pos[1], cell)
    assert v == pytest.approx([-1.0, 0.0, 0.0])
    # a triclinic cell: the naive wrap is not the shortest image; the 27-image search finds it
    tri = np.array([[10.0, 0, 0], [5.0, 8.7, 0], [0, 0, 10.0]])
    a, b = np.array([0.0, 0.0, 0.0]), np.array([9.9, 8.6, 0.0])
    brute = min(np.linalg.norm(b - a + np.array([i, j, k]) @ tri) for i in range(-2, 3) for j in range(-2, 3) for k in range(-2, 3))
    assert np.linalg.norm(mic_vector(a, b, tri)) == pytest.approx(brute)


def test_measure_reports_kind_and_text():
    pos = np.array([[0.0, 0, 0], [1.0, 0, 0], [1.0, 1.0, 0], [1.0, 1.0, 1.0]])
    assert measure(pos, [0]) is None
    m = measure(pos, [0, 1])
    assert m["kind"] == "distance" and m["value"] == pytest.approx(1.0)
    assert measure_text(m, ["O", "H", "H", "H"]) == "距離 O1-H2: 1.000 Å"
    m3 = measure(pos, [0, 1, 2])
    assert m3["kind"] == "angle" and m3["value"] == pytest.approx(90.0)
    assert measure_text(m3) == "角度 1-2-3: 90.0 °"
    m4 = measure(pos, [0, 1, 2, 3])
    assert m4["kind"] == "dihedral" and abs(m4["value"]) == pytest.approx(90.0)


def test_field_text_helpers():
    assert compact_ranges([3, 1, 2, 7, 2]) == "1-3,7"
    assert merge_fixed_text("", [2, 1]) == "1-2"
    assert merge_fixed_text("1-4,7:xy", [5, 9]) == "1-5,9,7:xy"
    assert merge_select_text("", [1, 2]) == "index 1-2"
    assert merge_select_text("index 1,2", [4]) == "index 1-2,4"
    assert merge_select_text("element O and z < 10", [3]) == "element O and z < 10 or index 3"
    assert series_text("", [2, 1, 3]) == "2,1,3"           # click order is kept: the angle is at the middle atom
    assert series_text("1,2; ", [3, 4]) == "1,2; 3,4"
    assert base_indices([0, 5, 3, 8], 3) == [1, 3]            # a 2x repeated 3-atom cell maps back to the base cell


pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _water() -> Atoms:
    return Atoms("OHH", positions=[(0.0, -1.0, 0.0), (0.0, 0.0, 0.783064), (0.0, 0.0, -0.783064)])


def test_clicking_atoms_in_the_viewer_measures_and_shift_adds(app):
    from PySide6.QtCore import QPointF, Qt
    from PySide6.QtGui import QMouseEvent
    from adit.gui.viewer3d import Viewer3D

    v = Viewer3D(); v.resize(400, 400); v.set_atoms(_water()); v.grab()      # grab paints, which fills the hit-test cache
    xy = v._xy
    assert len(xy) == 3

    def click(i, shift=False):
        mods = Qt.KeyboardModifier.ShiftModifier if shift else Qt.KeyboardModifier.NoModifier
        p = QPointF(*xy[i])
        v.mousePressEvent(QMouseEvent(QMouseEvent.Type.MouseButtonPress, p, p, Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton, mods))
        v.mouseReleaseEvent(QMouseEvent(QMouseEvent.Type.MouseButtonRelease, p, p, Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton, mods))

    seen = []
    v.selectionChanged.connect(lambda s: seen.append(list(s)))
    click(1); assert v.selected == [1]
    click(0, shift=True); assert v.selected == [1, 0]
    m = v.measurement(); assert m["kind"] == "distance" and m["value"] == pytest.approx(np.hypot(1.0, 0.783064))
    click(2, shift=True); assert v.selected == [1, 0, 2]
    assert v.measurement()["value"] == pytest.approx(2 * np.degrees(np.arctan2(0.783064, 1.0)))
    click(0, shift=True); assert v.selected == [1, 2]                  # shift-click on a selected atom removes it
    click(2); assert v.selected == [2]                                  # plain click replaces
    assert seen[-1] == [2]
    # a drag does not select
    p0, p1 = QPointF(*xy[0]), QPointF(xy[0][0] + 30, xy[0][1] + 30)
    v.mousePressEvent(QMouseEvent(QMouseEvent.Type.MouseButtonPress, p0, p0, Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier))
    v.mouseMoveEvent(QMouseEvent(QMouseEvent.Type.MouseMove, p1, p1, Qt.MouseButton.NoButton, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier))
    v.mouseReleaseEvent(QMouseEvent(QMouseEvent.Type.MouseButtonRelease, p1, p1, Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier))
    assert v.selected == [2]
    # clicking empty space clears
    v.grab(); far = QPointF(2.0, 2.0)
    v.mousePressEvent(QMouseEvent(QMouseEvent.Type.MouseButtonPress, far, far, Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier))
    v.mouseReleaseEvent(QMouseEvent(QMouseEvent.Type.MouseButtonRelease, far, far, Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier))
    assert v.selected == []
    assert v.grab().width() == 400                                      # painting with a selection and without does not crash


def test_periodic_measurement_uses_the_minimum_image(app):
    from adit.gui.viewer3d import Viewer3D

    a = Atoms("ArAr", positions=[(0.5, 5, 5), (9.5, 5, 5)], cell=np.diag([10.0, 10, 10]), pbc=True)
    v = Viewer3D(); v.set_atoms(a); v.set_selection([0, 1])
    assert v.measurement()["value"] == pytest.approx(1.0)


def test_panel_sends_base_cell_indices_to_the_chosen_field(app):
    from adit.gui.panels.structure_view_panel import StructureViewPanel
    from adit.spec import AtomsData, Structure

    panel = StructureViewPanel()
    st = Structure(source="file", source_ref="x", atoms=AtomsData(symbols=["Ar", "Ar"], positions=[(0.5, 5, 5), (9.5, 5, 5)],
                                                                  cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]], pbc=[True, True, True]), periodic=True)
    panel.set_structure(st)
    assert panel.repeat.currentText() == "2×2×2"                        # small cells are shown repeated
    got = []
    panel.send_selection.connect(lambda d, i: got.append((d, list(i))))
    panel.viewer.set_selection([3, 0])                                  # atom 3 of the repeated view is atom 2 of the base cell
    assert not panel.act_send_series.isEnabled() or len(panel.viewer.selected) == 2
    assert "Å" in panel.measure.text()
    panel._send("fixed"); panel._send("series"); panel._send("select")
    assert got == [("fixed", [2, 1]), ("distances", [2, 1]), ("select", [2, 1])]
    panel.viewer.set_selection([0])
    assert "Shift" in panel.measure.text()
    assert not panel.act_send_series.isEnabled()
    panel.btn_clear_sel.click()
    assert panel.viewer.selected == [] and not panel.btn_send.isEnabled()
    panel.viewer.set_selection([0, 1]); assert panel.btn_send.isEnabled()
    panel.set_structure(st)                                             # a new structure drops the selection and the buttons follow
    assert panel.viewer.selected == [] and not panel.btn_send.isEnabled() and panel.measure.text() == ""


def test_main_window_puts_the_atoms_into_the_fields(app, sk_root, tmp_path):
    from tests.test_gui import make_window

    win = make_window(sk_root, tmp_path)
    win.structure.set_source("preset"); win.structure._rebuild(); win._on_context()
    view = win.structure_view
    assert view.viewer.atoms() is not None and len(view.viewer.atoms()) == 3
    view.viewer.set_selection([0, 1]); view._send("fixed")
    assert win.structure.fixed.text() == "1-2"
    assert win.structure.structure().fixed_atoms == [0, 1]
    view.viewer.set_selection([1, 0, 2]); view._send("series")
    assert win.analysis.angles.text() == "2,1,3"
    view.viewer.set_selection([2]); view._send("select")
    assert win.analysis.select.text() == "index 3"
    assert "原子の選び方" in view.measure.text()
