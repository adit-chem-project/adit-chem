
from __future__ import annotations

import os

import numpy as np
import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from ase import Atoms  # noqa: E402
from ase.build import molecule  # noqa: E402
from ase.io import write  # noqa: E402
from PySide6.QtCore import QEventLoop, QTimer  # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from tests.conftest import water_spec  # noqa: E402
from adit.builder import Recipe, file_base, recipe_structure  # noqa: E402
from adit.builder.model import (Adsorb, Fix, Remove, Selection, Slab, SolventLayer, Solvate, Substitute, Supercell,  # noqa: E402
                               Vacuum)
from adit.config import default_config  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.panels.recipe_editor import INTERFACE, failed_step  # noqa: E402
from adit.mixture import Component  # noqa: E402
from adit.spec import AtomsData, Structure  # noqa: E402

SKEW_AB = [[11.0, 0.0, 0.0], [-4.4, 16.5, 0.0]]


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def win(app, sk_root, tmp_path, monkeypatch):
    for name in ("information", "critical"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(lambda *a, **k: None))
    w = MainWindow(default_config(sk_root=str(sk_root)), tmp_path / "cluster.toml")
    w.method.sk_set.setCurrentText("fake-1-0")
    w.runtime.outdir.setText(str(tmp_path / "out"))
    yield w
    w.close()


def add(panel, op: str):
    rc = panel.recipe
    rc.add_kind.setCurrentIndex(rc.add_kind.findData(op)); rc.btn_add.click()
    return rc.editors[-1] if op != INTERFACE else None


def build_and_wait(panel, timeout_ms: int = 120_000) -> None:
    loop = QEventLoop()
    panel.built.connect(loop.quit)
    QTimer.singleShot(timeout_ms, loop.quit)
    panel.recipe.btn_build.click()
    if panel.is_building():
        loop.exec()
    panel.built.disconnect(loop.quit)
    assert not panel.is_building(), "時間内に作り終わらなかった"


def spin_events(ms: int) -> None:
    loop = QEventLoop(); QTimer.singleShot(ms, loop.quit); loop.exec()


def rows(panel) -> list[str]:
    return [panel.recipe.list.item(i).text() for i in range(panel.recipe.list.count())]


def test_add_move_remove(win):
    p = win.structure
    p.set_source("bulk")
    assert p.recipe.count() == 0 and p.structure().source == "bulk"
    for op in ("supercell", "slab", "fix"):
        add(p, op)
    assert [r.split("  —  ")[0] for r in rows(p)] == ["1. 超格子", "2. 面で切る", "3. 固定"]
    assert p.recipe.stack.currentWidget() is p.recipe.editors[2]
    p.recipe.btn_up.click()
    assert [e.op for e in p.recipe.editors] == ["supercell", "fix", "slab"] and p.recipe.list.currentRow() == 1
    assert p.recipe.stack.currentWidget() is p.recipe.editors[1] and rows(p)[1].startswith("2. 固定")
    p.recipe.list.setCurrentRow(0); p.recipe.btn_down.click()
    assert [e.op for e in p.recipe.editors] == ["fix", "supercell", "slab"]
    p.recipe.btn_del.click()
    assert [e.op for e in p.recipe.editors] == ["fix", "slab"] and len(rows(p)) == 2
    assert p.recipe.stack.count() == 2
    assert not p.fixed.isEnabled()
    assert p.structure() is None and "作る" in p.error()
    while p.recipe.count():
        p.recipe.list.setCurrentRow(0); p.recipe.btn_del.click()
    assert p.fixed.isEnabled() and p.structure() is not None and p.structure().source == "bulk"


def test_each_step_fields_to_json(win):
    p = win.structure
    p.set_source("bulk"); p.bulk_cubic.setChecked(True)
    ed = add(p, "supercell"); ed.mode.setCurrentIndex(1); ed.mat[0][1].setValue(1); ed.mat[1][0].setValue(-1)
    assert ed.step() == Supercell(matrix=[[1, 1, 0], [-1, 1, 0], [0, 0, 1]])
    ed.mode.setCurrentIndex(0); ed.rep[1].setValue(3)
    assert ed.step() == Supercell(repeat=(2, 3, 1))
    ed = add(p, "slab")
    for w, v in zip(ed.hkl, (1, 1, 1)):
        w.setValue(v)
    ed.layers.setValue(2); ed.vacuum.setValue(8.0)
    assert ed.step() == Slab(miller=(1, 1, 1), layers=2, vacuum=8.0, termination=0)
    ed = add(p, "vacuum"); ed.thickness.setValue(20.0)
    assert ed.step() == Vacuum(axis=2, thickness=20.0)
    ed = add(p, "remove"); ed.sel.elements.setText("Si"); ed.sel.z_min.setText("1.5"); ed.count.setValue(2); ed.seed.setValue(4)
    assert ed.step() == Remove(where=Selection(elements=["Si"], z_min=1.5), count=2, seed=4)
    ed.sel.z_min.setText("abc")
    with pytest.raises(ValueError):
        ed.step()
    with pytest.raises(ValueError, match="手順 4"):
        p.recipe.steps()
    ed.sel.z_min.setText("")
    ed = add(p, "substitute"); ed.mode.setCurrentIndex(1); ed.fraction.setValue(0.25); ed.to.setCurrentText("Ge")
    assert ed.step() == Substitute(fraction=0.25, to="Ge")
    ed = add(p, "adsorb")
    ed.kind.setCurrentIndex(1); ed.ref.smiles.setText("[C-]#[O+]")
    ed.place.setCurrentIndex(ed.place.findData("above_atom")); ed.atom.setValue(3); ed.height.setValue(1.9); ed.down.setValue(1)
    assert ed.step() == Adsorb(molecule={"kind": "smiles", "ref": "[C-]#[O+]"}, above_atom=2, height=1.9, down_atom=0)
    ed = add(p, "solvent_layer")
    ed.mix.add_component(Component(kind="smiles", ref="[Li+]", count=2, charge=1, label="Li+"))
    ed.mode.setCurrentIndex(1); ed.thickness.setValue(15.0); ed.gap.setValue(2.5); ed.vacuum.setValue(0.0)
    s = ed.step()
    assert isinstance(s, SolventLayer) and s.thickness == 15.0 and s.gap == 2.5 and [c.ref for c in s.components] == ["H2O", "[Li+]"]
    ed = add(p, "solvate")
    assert ed.step() == Solvate(components=[Component(ref="H2O", count=0, label="H2O")])
    assert ed.mix.table.cellWidget(0, 2).text() == "自動"
    ed = add(p, "fix"); ed.layers.setValue(2); ed.sel.elements.setText("Si")
    assert ed.step() == Fix(bottom_layers=2, where=Selection(elements=["Si"]))
    add(p, "box")
    rec = p.current_recipe()
    assert [s.op for s in rec.steps] == ["supercell", "slab", "vacuum", "remove", "substitute", "adsorb", "solvent_layer", "solvate", "fix", "box"]
    assert rec.base.source == "bulk" and rec.base.ref == "Si cubic"
    assert Recipe.from_ref(rec.to_ref()) == rec
    assert p.charge.value() == 2
    assert "1. 超格子  —  2 × 3 × 1" in rows(p)[0]


def test_recipe_spec_roundtrip(win, tmp_path):
    slab = tmp_path / "slab.xyz"; write(slab, molecule("H2O"))
    rec = Recipe.model_validate({"base": file_base(slab).model_dump(), "steps": [
        {"op": "supercell", "matrix": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]},
        {"op": "slab", "miller": [1, 1, 0], "layers": 4, "vacuum": 12.5, "termination": 1},
        {"op": "vacuum", "axis": 1, "thickness": 7.25},
        {"op": "box", "padding": 3.5, "max_multiple": 4},
        {"op": "adsorb", "molecule": {"kind": "preset", "ref": "CO"}, "site": "fcc", "height": 1.25, "down_atom": 1},
        {"op": "remove", "where": {"elements": ["Li"], "z_min": 4.5, "indices": [3, 5]}, "fraction": 0.5, "seed": 9},
        {"op": "substitute", "where": {"z_max": 2.0}, "count": 2, "to": "Mg"},
        {"op": "solvent_layer", "components": [{"kind": "smiles", "ref": "O=C1OCCO1", "count": 32, "label": "C3H4O3"},
                                               {"kind": "smiles", "ref": "[Li+]", "count": 3, "charge": 1}],
         "density_g_cm3": 1.194, "gap": 2.2, "vacuum": 0.0, "min_distance": 1.9, "seed": 3, "max_tries": 900},
        {"op": "solvate", "components": [{"kind": "preset", "ref": "H2O", "count": 0}], "padding": 6.0},
        {"op": "fix", "bottom_layers": 2, "where": {"elements": ["Ni"]}, "layer_tolerance": 0.3}]})
    st = Structure(source="recipe", source_ref=rec.to_ref(), atoms=AtomsData.from_ase(molecule("H2O")), charge=3, fixed_atoms=[0])
    win.apply_spec(water_spec(structure=st))
    p = win.structure
    assert p.current_source() == "file" and p.recipe.count() == 10
    assert p.current_recipe() == rec
    assert win.current_spec().structure == st
    assert not p.fixed.isEnabled() and p.charge.value() == 3
    assert p.recipe.editors[1].term.currentData() == 1
    for base in ({"source": "2d", "ref": {"kind": "mx2", "formula": "WS2", "size": [3, 2, 1], "thickness": 3.1}},
                 {"source": "cluster", "ref": {"kind": "decahedron", "symbol": "Au", "p": 3, "q": 2, "r": 1, "lattice_constant": 4.08}},
                 {"source": "polymer", "ref": {"unit": "*CC(*)c1ccccc1", "n": 4, "seed": 2}}):
        r2 = Recipe.model_validate({"base": base, "steps": [{"op": "box", "padding": 4.0}]})
        s2 = Structure(source="recipe", source_ref=r2.to_ref(), atoms=AtomsData.from_ase(molecule("H2O")))
        win.apply_spec(water_spec(structure=s2))
        assert p.current_source() == base["source"]
        assert p.current_recipe() == r2 and win.current_spec().structure == s2


def test_build_in_thread_and_cancel(win):
    p = win.structure
    p.set_source("bulk"); p.bulk_cubic.setChecked(True)
    ed = add(p, "supercell"); ed.rep[0].setValue(2); ed.rep[1].setValue(1)
    ed = add(p, "remove"); ed.sel.elements.setText("Si"); ed.count.setValue(1); ed.seed.setValue(3)
    p.recipe.changed.emit()
    assert p.structure() is None
    build_and_wait(p)
    st = p.structure()
    assert st is not None and st.source == "recipe" and len(st.atoms.symbols) == 15
    lines = p.recipe.log.text().splitlines()
    assert len(lines) == 3 and "15 原子" in lines[2] and lines[0].startswith("0. 土台")
    assert win.structure_view._structure is not None
    assert Recipe.from_ref(st.source_ref) == p.current_recipe()
    win.refresh_preview()
    assert win.current_spec().structure.source == "recipe"
    ed.count.setValue(2)
    assert p.structure() is None and "作り直" in p.error()
    ed.count.setValue(1)
    assert p.structure() is not None
    ed.count.setValue(2)
    p.recipe.btn_build.click(); assert p.is_building() and p.recipe.btn_cancel.isVisible() is not None
    p.recipe.btn_cancel.click()
    spin_events(1500)
    assert not p.is_building() and p.structure() is None and "中止" in p.recipe.status.text()


def test_error_selects_the_failing_step(win):
    p = win.structure
    p.set_source("bulk"); p.bulk_cubic.setChecked(True)
    add(p, "supercell")
    ed = add(p, "remove"); ed.count.setValue(100)
    add(p, "fix")
    p.recipe.list.setCurrentRow(0)
    build_and_wait(p)
    assert p.structure() is None and p.error().startswith("手順 2 (原子を抜く)")
    assert p.recipe.list.currentRow() == 1 and p.recipe.stack.currentWidget() is ed
    assert p.recipe.status.objectName() == "status_ng"
    assert failed_step("after step 3 (slab), atoms 1 (Ni) and 2 (O) are 0.20 Å apart") == 3 and failed_step("手順 0 (土台) のあとで") is None


def _oblique_slab(tmp_path) -> str:
    cell = np.array(SKEW_AB + [[0.0, 0.0, 20.0]])
    sym, pos = [], []
    for layer, (el, z) in enumerate([("Ni", 0.0), ("O", 2.1), ("Ni", 4.2)]):
        for i in range(4):
            for j in range(6):
                sym.append(el); pos.append((i / 4 + (layer % 2) * 0.125) * cell[0] + (j / 6) * cell[1] + [0, 0, z])
    path = tmp_path / "oblique_slab.extxyz"
    write(path, Atoms(sym, positions=pos, cell=cell, pbc=True))
    return str(path)


def test_skew_cell_interface_from_the_screen(win, tmp_path):
    pytest.importorskip("rdkit")
    p = win.structure
    p.set_source("file"); p.file.setText(_oblique_slab(tmp_path)); p.file.editingFinished.emit()
    assert p.structure() is not None and len(p.structure().atoms.symbols) == 72
    add(p, INTERFACE)
    assert [e.op for e in p.recipe.editors] == ["supercell", "solvent_layer", "fix"]
    assert p.recipe.editors[0]._model.fit_components
    sl = p.recipe.editors[1]
    labels = [sl.mix.table.item(r, 4).text() for r in range(sl.mix.table.rowCount())]
    assert labels == ["C3H4O3", "Li+", "PF6-"]
    assert sl.vacuum.value() == 0.0
    sl.density.setValue(1.2)
    sl.conc.setValue(1.0)
    got = []
    for r in (1, 2, 1):
        sl.mix.table.selectRow(r)
        msg = sl.fill_count_from_concentration()
        assert "mol/L" in msg
        got.append(sl.mix.table.cellWidget(r, 2).value())
    assert got == [2, 3, 3], got
    counts = [sl.mix.table.cellWidget(r, 2).value() for r in range(3)]
    assert counts == [32, 3, 3], counts
    assert p.charge.value() == 0
    build_and_wait(p)
    st = p.structure()
    assert st is not None, p.error()
    a = st.atoms.to_ase()
    assert len(a) == 72 + 32 * 10 + 3 + 3 * 7
    s = a.get_chemical_symbols()
    assert (s.count("Li"), s.count("P"), s.count("F"), s.count("Ni")) == (3, 3, 18, 48)
    assert a.cell.lengths()[0] == pytest.approx(11.0000) and a.cell.angles()[2] == pytest.approx(104.9, abs=0.1)
    assert len(st.fixed_atoms) == 24 and set(np.array(s)[st.fixed_atoms]) == {"Ni"}
    assert st.charge == 0
    lines = p.recipe.log.text().splitlines()
    assert len(lines) == 4 and "416" in lines[2]
    rec = Recipe.from_ref(st.source_ref)
    assert rec.base.source == "file" and len(rec.base.sha256) == 64
    spec = win.current_spec()
    win.apply_spec(spec)
    assert win.current_spec().structure == spec.structure and p.current_recipe() == rec


def test_new_bases(win):
    p = win.structure
    p.set_source("2d")
    st = p.structure()
    assert st is not None and st.source == "recipe" and len(st.atoms.symbols) == 8
    p.set_source("cluster"); p.cl_shells.setValue(2)
    assert p.structure() is not None and len(p.structure().atoms.symbols) == 13 and not p.structure().periodic
    ed = add(p, "box"); ed.padding.setValue(4.0)
    build_and_wait(p)
    assert p.structure().periodic
    pytest.importorskip("rdkit")
    p.recipe.btn_del.click()
    p.set_source("polymer"); p.pl_unit.setText("*CC*"); p.pl_n.setValue(3); p.pl_unit.editingFinished.emit()
    assert p.structure() is None
    build_and_wait(p)
    assert p.structure() is not None and len(p.structure().atoms.symbols) == 20  # C6H14


def test_english_labels(app, sk_root, tmp_path):
    from adit.gui.i18n import set_language, translate_widgets
    set_language("en")
    try:
        w = MainWindow(default_config(sk_root=str(sk_root)), tmp_path / "c.toml"); translate_widgets(w)
        p = w.structure
        assert p.source.itemText(p.source.findData("cluster")) == "Nanoparticle"
        p.set_source("bulk")
        ed = add(p, "slab")
        texts = [ed.form.itemAt(i, ed.form.ItemRole.LabelRole).widget().text() for i in range(ed.form.rowCount())
                 if ed.form.itemAt(i, ed.form.ItemRole.LabelRole) is not None]
        assert "Miller indices (h k l)" in texts and "Vacuum per side [Å]" in texts
        assert rows(p)[0].startswith("1. slab")
        w.close()
    finally:
        set_language("ja")


def test_unticked_step_is_skipped_without_deleting_it(win):
    from PySide6.QtCore import Qt

    p = win.structure
    p.set_source("bulk"); p.bulk_cubic.setChecked(True)
    ed = add(p, "supercell"); ed.rep[0].setValue(2); ed.rep[1].setValue(1); ed.rep[2].setValue(1)
    add(p, "fix")
    assert not p.fixed.isEnabled()
    item = p.recipe.list.item(1)
    assert item.checkState() == Qt.CheckState.Checked
    item.setCheckState(Qt.CheckState.Unchecked)
    assert not p.recipe.editors[1].enabled and p.fixed.isEnabled() and "(無効" in rows(p)[1]
    assert p.recipe.count() == 2 and not p.current_recipe().steps[1].enabled
    build_and_wait(p)
    st = p.structure()
    assert st is not None and len(st.atoms.symbols) == 16 and st.fixed_atoms == []
    assert "2. 固定: 無効なので飛ばしました" in p.recipe.log.text()
    assert Recipe.from_ref(st.source_ref).steps[1].enabled is False
    p.recipe.set_enabled(2, True)
    assert item.checkState() == Qt.CheckState.Checked and not p.fixed.isEnabled() and p.structure() is None
    p.set_structure(st)
    assert p.recipe.list.item(1).checkState() == Qt.CheckState.Unchecked and p.fixed.isEnabled()
