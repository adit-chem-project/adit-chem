
from __future__ import annotations

import hashlib

import numpy as np
import pytest
from ase import Atoms
from ase.build import fcc111, molecule
from ase.io import write
from ase.neighborlist import neighbor_list

from adit.builder import (Recipe, RecipeError, build_recipe, file_base, has_op, interface_steps, list_terminations, recipe_structure,
                         salt_count, split_molecules)
from adit.builder.ops import planes
from adit.mixture import Component, MixtureError, MixtureSpec, build_mixture
from adit.structure import StructureError, build_structure

SKEW_AB = [[11.0, 0.0, 0.0], [-4.4, 16.5, 0.0]]
EC = "O=C1OCCO1"


def R(base: dict, *steps: dict) -> Recipe:
    return Recipe.model_validate({"base": base, "steps": list(steps)})


def run(base: dict, *steps: dict) -> Atoms:
    return build_recipe(R(base, *steps))[0]


def min_inter(atoms: Atoms, cutoff: float = 3.5) -> float:
    i, j, d = neighbor_list("ijd", atoms, cutoff)
    mol = atoms.get_array("molecule")
    s = mol[i] != mol[j]
    return float(d[s].min()) if s.any() else float("inf")


def sorted_frac(atoms: Atoms) -> np.ndarray:
    f = np.round(atoms.get_scaled_positions(wrap=True), 6) % 1.0
    return f[np.lexsort(f.T[::-1])]


def density(atoms: Atoms) -> float:
    return float(atoms.get_masses().sum()) / atoms.get_volume()


def test_supercell_repeat_and_matrix_agree():
    si = {"source": "bulk", "ref": "Si diamond 5.43 cubic"}
    a8 = run(si)
    a = run(si, {"op": "supercell", "repeat": [2, 2, 2]})
    b = run(si, {"op": "supercell", "matrix": [[2, 0, 0], [0, 2, 0], [0, 0, 2]]})
    assert len(a8) == 8 and len(a) == len(b) == 64
    assert abs(density(a) - density(a8)) < 1e-9
    assert np.allclose(np.asarray(a.cell), np.asarray(b.cell)) and np.allclose(sorted_frac(a), sorted_frac(b), atol=1e-5)
    c = run(si, {"op": "supercell", "matrix": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]})
    assert len(c) == 16 and abs(density(c) - density(a8)) < 1e-9
    with pytest.raises(RecipeError):
        run(si, {"op": "supercell", "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, -1]]})


@pytest.mark.parametrize("step", [{"op": "supercell", "repeat": [2, 1, 1]},
                                  {"op": "supercell", "matrix": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]}, {"op": "box"}])
def test_supercell_keeps_molecules_separate(step):
    mix = MixtureSpec(components=[Component(ref="H2O", count=12)], cell=[[10.0, 0, 0], [5.0, 8.66, 0], [0, 0, 10.0]], seed=0)
    base = {"source": "mixture", "ref": mix.to_ref()}
    a, logs = build_recipe(R(base, step))
    mol = a.get_array("molecule")
    assert len(set(mol.tolist())) == len(a) // 3 and np.bincount(mol).tolist() == [3] * (len(a) // 3)
    for k in set(mol.tolist()):
        p = a.positions[mol == k]
        assert np.linalg.norm(p - p.mean(axis=0), axis=1).max() < 1.2
    assert split_molecules(a) == [] and logs[-1].note == ""
    assert logs[-1].min_distance >= 2.0 - 1e-6


def test_mixture_cell_none_keeps_the_old_coordinates():
    def fp(spec):
        return hashlib.sha256(np.round(build_mixture(spec).get_positions(), 6).tobytes()).hexdigest()[:16]

    assert fp(MixtureSpec(components=[Component(ref="H2O", count=20, label="水"), Component(ref="NH3", count=2)],
                          density_g_cm3=0.8, seed=1)) == "715e611ac8daafd3"
    assert fp(MixtureSpec(components=[Component(ref="H2O", count=32)], density_g_cm3=1.0, seed=0)) == "eef3f0e56f296fd9"
    assert "cell" not in MixtureSpec(components=[Component(ref="H2O")]).to_ref()


@pytest.mark.parametrize("cell", [
    [[12.0, 0, 0], [0, 14.0, 0], [0, 0, 16.0]],
    SKEW_AB + [[0.0, 0.0, 20.0]],
])
def test_mixture_in_any_cell(cell):
    m = MixtureSpec(components=[Component(ref="H2O", count=80)], cell=cell, seed=0)
    a = build_mixture(m)
    assert len(a) == 240 and np.allclose(np.asarray(a.cell), cell)
    assert min_inter(a) >= 2.0 - 1e-6
    assert split_molecules(a) == []
    b = build_mixture(m)
    assert np.allclose(a.get_positions(), b.get_positions())


def test_mixture_z_band():
    cell = [[15.0, 0, 0], [0, 15.0, 0], [0, 0, 30.0]]
    a = build_mixture(MixtureSpec(components=[Component(ref="H2O", count=40)], cell=cell, z_band=(5.0, 17.0), seed=2))
    z = a.positions[:, 2]
    assert z.min() >= 5.0 - 1e-9 and z.max() <= 17.0 + 1e-9 and min_inter(a) >= 2.0 - 1e-6


AU = {"source": "surface", "ref": "fcc111 Au 4x4x4 vacuum=10"}


def _water_layer(**kw) -> dict:
    return {"op": "solvent_layer", "components": [{"kind": "preset", "ref": "H2O", "count": 30}], "gap": 2.5, **kw}


def test_interface_au111_water():
    from adit.structure import from_surface
    slab = from_surface(AU["ref"])
    a = run(AU, _water_layer())
    b = run(AU, _water_layer(vacuum=10.0))
    assert len(a) == 64 + 90
    assert min_inter(a) >= 2.0 - 1e-6 and min_inter(b) >= 2.0 - 1e-6
    lo, hi = a.info["solvent_band"]
    zs = a.positions[64:, 2]
    assert zs.min() >= lo - 1e-9 and zs.max() <= hi + 1e-9
    assert lo == pytest.approx(a.positions[:64, 2].max() + 2.5)
    assert np.allclose(np.asarray(a.cell)[:2], np.asarray(slab.cell)[:2])
    assert np.allclose(np.asarray(a.cell)[2, :2], 0)
    assert a.cell.lengths()[2] + 10.0 == pytest.approx(b.cell.lengths()[2])
    assert a.cell.lengths()[2] - hi == pytest.approx(2.5)
    assert np.allclose(a.get_positions(), run(AU, _water_layer()).get_positions())
    assert not np.allclose(a.get_positions(), run(AU, _water_layer(seed=1)).get_positions())
    shift = slab.positions[0] - a.positions[0]
    assert np.allclose(slab.positions - shift, a.positions[:64], atol=1e-9)


def _oblique_slab(tmp_path) -> str:
    cell = np.array(SKEW_AB + [[0.0, 0.0, 20.0]])
    sym, pos = [], []
    for layer, (el, z) in enumerate([("Ni", 0.0), ("O", 2.1), ("Ni", 4.2)]):
        for i in range(4):
            for j in range(6):
                sym.append(el); pos.append((i / 4 + (layer % 2) * 0.125) * cell[0] + (j / 6) * cell[1] + [0, 0, z])
    p = tmp_path / "oblique_slab.extxyz"
    write(p, Atoms(sym, positions=pos, cell=cell, pbc=True))
    return str(p)


def test_skew_cell_interface_both_sides_solution(tmp_path):
    pytest.importorskip("rdkit")
    path = _oblique_slab(tmp_path)
    base = file_base(path).model_dump()
    comps = [{"kind": "smiles", "ref": EC, "count": 32}, {"kind": "smiles", "ref": "[Li+]", "count": 3, "charge": 1},
             {"kind": "smiles", "ref": "F[P-](F)(F)(F)(F)F", "count": 3, "charge": -1}]
    rec = R(base, {"op": "solvent_layer", "components": comps, "density_g_cm3": 1.2, "gap": 2.0, "vacuum": 0.0})
    a, logs = build_recipe(rec)
    assert len(a) == 72 + 32 * 10 + 3 + 3 * 7
    s = a.get_chemical_symbols()
    assert (s.count("Li"), s.count("P"), s.count("F"), s.count("C")) == (3, 3, 18, 96)
    assert rec.total_charge() == 0
    lengths, angles = a.cell.lengths(), a.cell.angles()
    assert lengths[0] == pytest.approx(11.0000) and lengths[1] == pytest.approx(17.0766, abs=1e-3)
    assert angles[2] == pytest.approx(104.9, abs=0.1) and angles[0] == pytest.approx(90) and angles[1] == pytest.approx(90)
    assert min_inter(a) >= 2.0 - 1e-6
    assert split_molecules(a) == []
    lo, hi = a.info["solvent_band"]
    assert a.cell.lengths()[2] - hi == pytest.approx(2.0)
    assert all(np.isfinite(l.seconds) for l in logs) and logs[-1].n_atoms == len(a)


def test_interface_helper_cuts_before_fitting_cross_section():
    pytest.importorskip("rdkit")
    components = [Component(kind="smiles", ref=EC, count=1, label="C3H4O3")]
    steps = interface_steps(False, components)
    assert [s.op for s in steps] == ["slab", "supercell", "solvent_layer", "fix"]
    atoms, logs = build_recipe(Recipe(base={"source": "bulk", "ref": "Al fcc 4.05 cubic"}, steps=steps))
    assert len(atoms) > 10 and "面内" in logs[2].note
    assert atoms.cell.lengths()[0] > 2.86 and atoms.cell.lengths()[1] > 2.86


def test_supercell_before_slab_records_miller_index_warning():
    rec = R({"source": "bulk", "ref": "Al fcc 4.05 cubic"},
            {"op": "supercell", "repeat": [2, 1, 1]},
            {"op": "slab", "miller": [1, 0, 0], "layers": 2})
    _atoms, logs = build_recipe(rec)
    assert "ミラー指数" in logs[2].note


def test_solvate_keeps_solute_fixed():
    a = run({"source": "preset", "ref": "C6H6"},
            {"op": "solvate", "components": [{"kind": "preset", "ref": "H2O", "count": 0}], "padding": 6.0, "density_g_cm3": 1.0})
    ref = molecule("C6H6")
    ref.set_cell(np.diag(np.ptp(ref.positions, axis=0) + 12.0)); ref.center()
    assert np.allclose(a.positions[:12], ref.positions)
    n_water = (len(a) - 12) // 3
    assert n_water > 10 and a.get_chemical_symbols()[:12] == ref.get_chemical_symbols()
    assert min_inter(a) >= 2.0 - 1e-6


def test_slab_al111_matches_fcc111_area_density():
    a = run({"source": "bulk", "ref": "Al fcc 4.05 cubic"}, {"op": "slab", "miller": [1, 1, 1], "layers": 2, "vacuum": 8.0})
    g = planes(a.positions[:, 2], 0.2)
    area = np.linalg.norm(np.cross(a.cell[0], a.cell[1]))
    ref = fcc111("Al", (1, 1, 3), a=4.05, vacuum=8.0)
    ref_area = np.linalg.norm(np.cross(ref.cell[0], ref.cell[1]))
    assert len(g[0]) / area == pytest.approx(1 / ref_area)
    assert np.allclose(np.asarray(a.cell)[2, :2], 0) and np.allclose(np.asarray(a.cell)[:2, 2], 0)
    assert all(a.pbc)


def test_slab_nio100_is_stoichiometric_and_terminations():
    nio = {"source": "bulk", "ref": "NiO rocksalt 4.17 cubic"}
    a = run(nio, {"op": "slab", "miller": [1, 0, 0], "layers": 2, "vacuum": 8.0})
    s = np.array(a.get_chemical_symbols())
    assert (s == "Ni").sum() == (s == "O").sum()
    for g in planes(a.positions[:, 2], 0.2):
        assert (s[g] == "Ni").sum() == (s[g] == "O").sum()
    from adit.structure import from_bulk
    terms = list_terminations(from_bulk("NiO rocksalt 4.17 cubic"), (1, 1, 1))
    assert len(terms) >= 2 and all(set(Atoms(t).get_chemical_symbols()) in ({"Ni"}, {"O"}) for t in terms)
    t0 = run(nio, {"op": "slab", "miller": [1, 1, 1], "layers": 2, "vacuum": 8.0, "termination": 0})
    t1 = run(nio, {"op": "slab", "miller": [1, 1, 1], "layers": 2, "vacuum": 8.0, "termination": 1})
    top = lambda x: x[planes(x.positions[:, 2], 0.2)[-1]].get_chemical_symbols()[0]
    assert top(t0) != top(t1) and t0.get_chemical_formula() == t1.get_chemical_formula()
    with pytest.raises(RecipeError):
        run(nio, {"op": "slab", "miller": [1, 1, 1], "layers": 2, "termination": 99})
    one = run(nio, {"op": "slab", "miller": [1, 1, 1], "layers": 1, "vacuum": 8.0})
    g = planes(one.positions[:, 2], 0.2)
    assert [one[x].get_chemical_formula() for x in g] in (["Ni4", "O4"], ["O4", "Ni4"])
    assert np.ptp(one.positions[:, 2]) == pytest.approx(4.17 / np.sqrt(3) / 2, abs=1e-6)


@pytest.mark.parametrize("where", [{"site": "ontop"}, {"above_atom": "top"}])
def test_adsorb_c_on_cu_at_1_9(where):
    cu = {"source": "surface", "ref": "fcc111 Cu 3x3x3 vacuum=8"}
    from adit.structure import from_surface
    slab = from_surface(cu["ref"])
    if "above_atom" in where:
        where = {"above_atom": int(np.argmax(slab.positions[:, 2]))}
    co = molecule("CO")
    c_index = co.get_chemical_symbols().index("C")
    a = run(cu, {"op": "adsorb", "molecule": {"kind": "preset", "ref": "CO"}, "height": 1.9, "down_atom": c_index, **where})
    c = len(slab) + c_index
    o = len(slab) + 1 - c_index
    d = a.get_distances(c, list(range(len(slab))), mic=True)
    assert d.min() == pytest.approx(1.9, abs=1e-6)
    assert a.positions[o, 2] > a.positions[c, 2]
    below = int(np.argmin(d))
    assert np.allclose(a.positions[below, :2], a.positions[c, :2]) or np.allclose(
        a.get_distance(c, below, mic=True, vector=True)[:2], 0, atol=1e-6)


def test_adsorb_overlap_is_stopped_with_reason():
    with pytest.raises(RecipeError):
        run({"source": "surface", "ref": "fcc111 Cu 2x2x2 vacuum=8"},
            {"op": "adsorb", "molecule": {"kind": "preset", "ref": "H2O"}, "above_atom": 7, "height": 0.1})


def test_remove_and_substitute():
    si64 = [{"source": "bulk", "ref": "Si diamond 5.43 cubic"}, {"op": "supercell", "repeat": [2, 2, 2]}]
    a = run(si64[0], si64[1], {"op": "remove", "where": {"elements": ["Si"]}, "count": 1, "seed": 3})
    b = run(si64[0], si64[1], {"op": "remove", "where": {"elements": ["Si"]}, "count": 1, "seed": 3})
    c = run(si64[0], si64[1], {"op": "remove", "where": {"elements": ["Si"]}, "count": 1, "seed": 4})
    assert len(a) == 63 and np.allclose(a.positions, b.positions) and not np.allclose(a.positions, c.positions)
    with pytest.raises(RecipeError):
        run(si64[0], si64[1], {"op": "remove", "where": {"elements": ["Si"]}, "count": 100})
    with pytest.raises(RecipeError):
        run(si64[0], si64[1], {"op": "remove", "where": {"elements": ["Li"]}, "count": 1})
    d = run(si64[0], si64[1], {"op": "substitute", "where": {"elements": ["Si"]}, "fraction": 0.25, "to": "Ge", "seed": 0})
    assert d.get_chemical_symbols().count("Ge") == 16 and len(d) == 64


def test_vacuum_box_and_fix():
    au = run(AU, {"op": "vacuum", "axis": 2, "thickness": 20.0})
    zspan = np.ptp(au.positions[:, 2])
    assert au.cell.lengths()[2] == pytest.approx(zspan + 20.0)
    rect = run({"source": "surface", "ref": "fcc111 Al 2x2x3 vacuum=5"}, {"op": "box"})
    assert np.allclose(rect.cell.angles(), 90) and len(rect) == 24
    assert abs(density(rect) - density(run({"source": "surface", "ref": "fcc111 Al 2x2x3 vacuum=5"}))) < 1e-9
    mol = run({"source": "preset", "ref": "H2O"}, {"op": "box", "padding": 5.0})
    assert all(mol.pbc) and np.allclose(mol.cell.angles(), 90)
    st, _ = recipe_structure(R(AU, {"op": "fix", "bottom_layers": 2}))
    assert len(st.fixed_atoms) == 32 and max(st.fixed_atoms) < 64
    z = np.array(st.atoms.positions)[:, 2]
    assert z[st.fixed_atoms].max() < np.sort(np.unique(np.round(z, 3)))[2]


def test_2d_nanotube_cluster():
    g = run({"source": "2d", "ref": {"kind": "graphene", "size": [2, 2, 1]}})
    assert len(g) == 8 and all(g.pbc) and g.cell.lengths()[2] == pytest.approx(20.0)
    assert len(run({"source": "2d", "ref": {"kind": "mx2", "size": [2, 2, 1]}})) == 12
    t = run({"source": "2d", "ref": {"kind": "nanotube", "n": 6, "m": 0, "length": 2, "vacuum": 5.0}})
    assert len(t) == 48 and all(t.pbc)
    r = run({"source": "2d", "ref": {"kind": "nanoribbon", "n": 3, "m": 2, "ribbon_type": "armchair"}})
    assert len(r) == 32 and all(r.pbc)
    ico = run({"source": "cluster", "ref": {"kind": "icosahedron", "symbol": "Cu", "shells": 4}})
    assert len(ico) == 147 and not any(ico.pbc)
    assert len(run({"source": "cluster", "ref": {"kind": "octahedron", "symbol": "Au", "length": 5, "cutoff": 2}})) == 55
    w = run({"source": "cluster", "ref": {"kind": "wulff", "symbol": "Cu", "size": 80}})
    assert 40 < len(w) < 160


def test_polymer_polyethylene():
    pytest.importorskip("rdkit")
    a = run({"source": "polymer", "ref": {"unit": "*CC*", "n": 10}})
    assert a.get_chemical_formula(mode="hill") == "C20H42"
    with pytest.raises(RecipeError):
        run({"source": "polymer", "ref": {"unit": "CC*", "n": 3}})
    with pytest.raises(RecipeError):
        run({"source": "polymer", "ref": {"unit": "*CC*", "n": 100000}})
    import time
    t = time.perf_counter()
    with pytest.raises(RecipeError):
        run({"source": "polymer", "ref": {"unit": "*CC*", "n": 100}})
    assert time.perf_counter() - t < 5.0
    ps = run({"source": "polymer", "ref": {"unit": "*CC(*)c1ccccc1", "n": 3}})
    assert ps.get_chemical_formula(mode="hill") == "C24H26"


def test_recipe_json_roundtrip_and_build_structure():
    from tests.conftest import water_spec
    from adit.spec import CalculationSpec

    rec = R(AU, _water_layer(seed=5), {"op": "fix", "bottom_layers": 1},
            {"op": "remove", "where": {"elements": ["Au"], "z_min": 5.0}, "count": 2, "seed": 1})
    ref = rec.to_ref()
    assert Recipe.from_ref(ref) == rec
    atoms, _ = build_recipe(rec)
    st = build_structure("recipe", ref)
    assert st.source == "recipe" and np.allclose(np.array(st.atoms.positions), atoms.positions)
    assert len(st.fixed_atoms) == 16
    assert build_structure("recipe", st.source_ref).atoms == st.atoms
    spec = water_spec(structure=st)
    back = CalculationSpec.from_json(spec.to_json())
    assert back.structure == st
    with pytest.raises(StructureError):
        Recipe.from_ref('{"base": {"source": "preset", "ref": "H2O"}, "steps": [{"op": "twist"}]}')
    with pytest.raises(StructureError):
        Recipe.from_ref('{"base": {"source": "preset", "ref": "H2O"}, "steps": [{"op": "vacuum", "thicknes": 3}]}')


def test_file_base_detects_changed_file(tmp_path):
    p = tmp_path / "m.xyz"
    write(p, molecule("H2O"))
    base = file_base(p)
    assert len(run(base.model_dump())) == 3
    write(p, molecule("NH3"))
    with pytest.raises(RecipeError):
        run(base.model_dump())


def test_salt_count():
    n, c = salt_count(1.0, 1000.0)
    assert n == 1 and c == pytest.approx(1 / (6.02214076e23 * 1000.0 * 1e-27))
    n, c = salt_count(1.0, 204.25 * 22.3)
    assert n == 3 and c == pytest.approx(3 / (6.02214076e23 * 204.25 * 22.3 * 1e-27))


def test_limits_stop_before_building():
    with pytest.raises(RecipeError):
        run({"source": "bulk", "ref": "Si diamond 5.43 cubic"}, {"op": "supercell", "repeat": [20, 20, 20]})
    with pytest.raises(RecipeError):
        run({"source": "cluster", "ref": {"kind": "icosahedron", "symbol": "Cu", "shells": 60}})
    with pytest.raises(MixtureError):
        build_mixture(MixtureSpec(components=[Component(ref="H2O", count=20000)]))
    with pytest.raises(RecipeError):
        run(AU, {"op": "solvent_layer", "components": [{"kind": "preset", "ref": "H2O", "count": 20000}]})


def test_solvent_layer_needs_a_slab_cell():
    with pytest.raises(RecipeError):
        run({"source": "bulk", "ref": "Si"}, _water_layer())
    with pytest.raises(RecipeError):
        run({"source": "preset", "ref": "H2O"}, _water_layer())


def test_disabled_step_is_skipped_and_only_false_is_written():
    rec = R(AU, {"op": "vacuum", "axis": 2, "thickness": 20.0, "enabled": False}, {"op": "fix", "bottom_layers": 2})
    assert not rec.steps[0].enabled and rec.steps[1].enabled
    ref = rec.to_ref()
    assert '"enabled":false' in ref and '"enabled":true' not in ref
    assert Recipe.from_ref(ref) == rec
    atoms, logs = build_recipe(rec)
    plain = run(AU, {"op": "fix", "bottom_layers": 2})
    assert np.allclose(atoms.cell, plain.cell) and len(atoms) == len(plain)
    assert logs[1].skipped and "飛ばしました" in logs[1].line() and not logs[2].skipped
    st, _ = recipe_structure(rec)
    assert len(st.fixed_atoms) == 32
    assert R(AU, _water_layer(enabled=False)).total_charge() == 0
    assert not has_op(rec.steps, "vacuum") and has_op(rec.steps, "fix")
