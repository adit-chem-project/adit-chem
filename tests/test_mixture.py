
from __future__ import annotations

import numpy as np
import pytest

from adit.mixture import Component, MixtureError, MixtureSpec, build_mixture, concentration_mol_per_l, mixture_text, parse_mixture_text
from adit.structure import build_structure


@pytest.mark.parametrize("language,expected", [
    ("ja", "スラブとの隙間 (gap) は帯の厚みを変えません"),
    ("en", "The slab gap does not change the band thickness"),
])
def test_band_packing_failure_distinguishes_gap_from_band_thickness(monkeypatch, language, expected):
    from adit import lang
    from adit.mixture import Packed, pack_failure_message

    monkeypatch.setattr(lang, "LANGUAGE", language)
    packed = Packed([], [], [], [], False, 1.83, (12.1, 87.3), cell=np.diag([10.0, 10.0, 100.0]))
    assert expected in pack_failure_message("C3H4O3 × 32", packed, 2.0)
    packed.band = None
    assert expected not in pack_failure_message("C3H4O3 × 32", packed, 2.0)


def _spec(**kw) -> MixtureSpec:
    comps = [Component(ref="H2O", count=20, label="水"), Component(ref="NH3", count=2)]
    return MixtureSpec(components=comps, density_g_cm3=0.8, seed=1, **kw)


def test_counts_and_min_distance():
    a = build_mixture(_spec())
    assert len(a) == 20 * 3 + 2 * 4 and all(a.pbc)
    tags = a.get_tags(); d = a.get_all_distances(mic=True)
    mol_id = a.get_array("molecule")
    other = mol_id[:, None] != mol_id[None, :]
    assert d[other].min() >= 2.0 - 1e-9
    assert set(tags) == {0, 1}


def test_seed_reproduces_and_changes():
    a, b, c = build_mixture(_spec()), build_mixture(_spec()), build_mixture(_spec(seed=2) if False else MixtureSpec(components=_spec().components, density_g_cm3=0.8, seed=2))
    assert np.allclose(a.get_positions(), b.get_positions()) and not np.allclose(a.get_positions(), c.get_positions())


def test_edge_box_and_concentration():
    m = MixtureSpec(components=[Component(ref="H2O", count=10)], box_a=12.0, seed=0)
    a = build_mixture(m)
    assert abs(a.cell.lengths()[0] - 12.0) < 1e-9
    assert abs(concentration_mol_per_l(10, 12.0) - 10 / (6.02214076e23 * 12.0 ** 3 * 1e-27)) < 1e-9


def test_total_charge_and_spec_roundtrip():
    m = MixtureSpec(components=[Component(ref="H2O", count=8), Component(kind="smiles", ref="[Na+]", count=2, charge=1, label="Na+"),
                                Component(kind="smiles", ref="[Cl-]", count=1, charge=-1, label="Cl-")], box_a=10.0)
    pytest.importorskip("rdkit")
    assert m.total_charge() == 1
    st = build_structure("mixture", m.to_ref(), charge=m.total_charge())
    assert st.source == "mixture" and st.periodic and st.charge == 1 and len(st.atoms.symbols) == 8 * 3 + 3
    st2 = build_structure("mixture", MixtureSpec.from_ref(st.source_ref).to_ref(), charge=1)
    assert st2.atoms.positions == st.atoms.positions


def test_text_roundtrip():
    comps = [Component(ref="H2O", count=32), Component(kind="smiles", ref="[Na+]", count=1, charge=1, label="Na+"), Component(kind="file", ref="/tmp/x.xyz", count=3)]
    assert parse_mixture_text(mixture_text(comps)) == comps
    assert parse_mixture_text("  # comment\nH2O * 4\n") == [Component(ref="H2O", count=4)]
    with pytest.raises(MixtureError):
        parse_mixture_text("H2O * many")


def test_cannot_place_gives_readable_error():
    with pytest.raises(MixtureError) as ex:
        build_mixture(MixtureSpec(components=[Component(ref="H2O", count=50)], box_a=5.0, max_tries=50))
    assert "H2O" in str(ex.value)


def test_sketch_to_smiles():
    pytest.importorskip("rdkit")
    pytest.importorskip("PySide6")
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from adit.gui.sketcher import SAtom, SBond, Sketch, SketchCanvas

    c = SketchCanvas(); c.add_template("benzene", 0, 0); c.sketch.atoms[0].elem = "N"
    assert c.sketch.to_smiles() == "c1ccncc1"
    sk = Sketch(atoms=[SAtom(0, 0, "C"), SAtom(44, 0, "O")], bonds=[SBond(0, 1, 2)])
    assert sk.to_smiles() == "C=O"
    back = Sketch.from_smiles("CC(=O)[O-]")
    assert len(back.atoms) == 4 and back.to_smiles() == "CC(=O)[O-]"


def test_packs_at_liquid_density():
    m = MixtureSpec(components=[Component(ref="H2O", count=32)], density_g_cm3=1.0, seed=0)
    a = build_mixture(m)
    assert len(a) == 96
    ids = a.get_array("molecule")
    d = a.get_all_distances(mic=True)
    assert d[ids[:, None] != ids[None, :]].min() >= 2.0 - 1e-9


def test_molecules_are_not_split_by_the_boundary():
    m = MixtureSpec(components=[Component(ref="H2O", count=24)], density_g_cm3=1.0, seed=1)
    a = build_mixture(m)
    box = float(a.cell.lengths()[0])
    pos, ids = a.get_positions(), a.get_array("molecule")
    for k in set(ids.tolist()):
        p = pos[ids == k]
        assert float(np.abs(p - p.mean(axis=0)).max()) < box / 2
