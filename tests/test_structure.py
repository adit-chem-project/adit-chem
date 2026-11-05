from pathlib import Path

import numpy as np
import pytest
from ase.build import molecule
from ase.io import write

from adit.structure import (StructureError, build_structure, from_file, from_preset,
                             from_smiles, has_rdkit, preset_matches, preset_names, preset_search_text)

REPO = Path(__file__).resolve().parent.parent


def test_preset_list_contains_water_and_is_sorted():
    names = preset_names()
    assert "H2O" in names and names == sorted(names)
    assert "ethanol" in preset_search_text("CH3CH2OH") and "エタノール" in preset_search_text("CH3CH2OH")
    assert preset_matches("CH3CH2OH", "ethanol") and not preset_matches("CH3OH", "ethanol")


def test_preset_water():
    atoms = from_preset("H2O")
    assert sorted(atoms.get_chemical_symbols()) == ["H", "H", "O"]
    assert not any(atoms.pbc)


def test_preset_unknown():
    with pytest.raises(StructureError):
        from_preset("unobtainium")


@pytest.mark.parametrize("fmt", ["xyz", "cif", "proteindatabank", "gen"])
def test_file_formats(tmp_path, fmt):
    ref = molecule("CH4")
    if fmt == "cif":
        ref.set_cell([8, 8, 8]); ref.pbc = True
    p = tmp_path / f"m.{fmt}"
    write(p, ref, format=fmt)
    atoms = from_file(p)
    assert sorted(atoms.get_chemical_symbols()) == ["C", "H", "H", "H", "H"]


def test_file_tutorial_gen():
    atoms = from_file(REPO / "examples" / "water" / "geometry.gen")
    assert atoms.get_chemical_symbols() == ["O", "H", "H"]
    assert np.isclose(atoms.get_distance(1, 2), 2 * 0.783064)


def test_file_missing_or_broken(tmp_path):
    with pytest.raises(StructureError):
        from_file(tmp_path / "nope.xyz")
    p = tmp_path / "broken.xyz"; p.write_text("not an xyz\n", encoding="utf-8")
    with pytest.raises(StructureError):
        from_file(p)


def test_file_errors_are_bilingual_and_empty_path_is_named(tmp_path):
    from adit.lang import LANGUAGE, set_language
    from adit.mixture import Component, MixtureSpec, build_mixture
    from adit.structure import from_bulk, from_surface

    before = LANGUAGE
    try:
        for lang, empty in (("ja", "構造ファイルを指定してください"), ("en", "choose a structure file")):
            set_language(lang)
            for blank in ("", "   ", None):
                with pytest.raises(StructureError) as ex:
                    from_file(blank)
                assert str(ex.value) == empty
            with pytest.raises(StructureError) as ex:
                build_mixture(MixtureSpec(components=[Component(kind="file", ref="", count=1)]))
            assert str(ex.value) == empty
            missing = tmp_path / "nope.xyz"
            with pytest.raises(StructureError) as ex:
                from_file(missing)
            assert str(missing) in str(ex.value)
            for bad in (lambda: from_bulk(""), lambda: from_bulk("Si diamond abc"), lambda: from_surface("fcc111 Al"),
                        lambda: build_structure("magic", "x")):
                with pytest.raises(StructureError) as ex:
                    bad()
                has_ja = any("぀" <= ch <= "ヿ" for ch in str(ex.value))
                assert has_ja == (lang == "ja"), str(ex.value)
    finally:
        set_language(before)


def test_file_multi_frame_takes_last(tmp_path):
    a, b = molecule("H2"), molecule("H2"); b.positions[1, 2] += 0.5
    p = tmp_path / "traj.xyz"; write(p, [a, b])
    assert np.isclose(from_file(p).get_distance(0, 1), b.get_distance(0, 1))


@pytest.mark.skipif(not has_rdkit(), reason="RDKit is missing")
def test_smiles_ethanol():
    atoms = from_smiles("CCO")
    assert sorted(atoms.get_chemical_symbols()) == sorted(["C", "C", "O"] + ["H"] * 6)
    d = atoms.get_all_distances(); np.fill_diagonal(d, 9.0)
    assert d.min() > 0.9


@pytest.mark.skipif(not has_rdkit(), reason="RDKit is missing")
def test_smiles_invalid():
    with pytest.raises(StructureError):
        from_smiles("C(C")


def test_build_structure_preset_and_bad_source():
    s = build_structure("preset", "H2O", charge=1, multiplicity=2)
    assert s.source == "preset" and s.charge == 1 and len(s.atoms.symbols) == 3
    with pytest.raises(StructureError):
        build_structure("magic", "x")


def test_bulk_defaults_from_reference_state():
    from adit.structure import default_bulk, from_bulk
    assert default_bulk("Si") == ("diamond", 5.43) and default_bulk("Al") == ("fcc", 4.05)
    a = from_bulk("Si")
    assert len(a) == 2 and all(a.pbc)
    a8 = from_bulk("Si diamond 5.43 cubic")
    assert len(a8) == 8 and abs(a8.cell.lengths()[0] - 5.43) < 1e-9
    with pytest.raises(StructureError):
        from_bulk("Si diamond abc")
    with pytest.raises(StructureError):
        from_bulk("Xx")


def test_surface_slab():
    from adit.structure import from_surface
    s = from_surface("fcc111 Al 2x2x3 vacuum=10")
    assert len(s) == 12 and tuple(s.pbc) == (True, True, True)
    assert s.cell.lengths()[2] > 20
    with pytest.raises(StructureError):
        from_surface("fcc111 Al")
    with pytest.raises(StructureError):
        from_surface("fcc999 Al 2x2x3")
    with pytest.raises(StructureError):
        from_surface("fcc111 Al 2x2x3 foo=1")


def test_build_structure_bulk_and_surface():
    b = build_structure("bulk", "Si")
    assert b.source == "bulk" and b.atoms.pbc == (True, True, True)
    s = build_structure("surface", "fcc111 Al 2x2x3")
    assert len(s.atoms.symbols) == 12
