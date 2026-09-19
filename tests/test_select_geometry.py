
import numpy as np
import pytest
from ase.build import bulk, molecule

from adit.analysis.geometry_series import (GeometryError, angle_series, dihedral_series, distance_series,
                                            kabsch_rmsd, parse_atom_list, rmsd_series, rmsf)
from adit.analysis.select import SelectionError, describe, parse_indices, select


def test_element_and_index_and_coordinate():
    w = molecule("H2O")
    assert select(w, "element O").tolist() == [0]
    assert select(w, "O").tolist() == [0]
    assert select(w, "element O H").tolist() == [0, 1, 2]
    assert select(w, "index 1-2").tolist() == [0, 1]
    assert select(w, None).tolist() == [0, 1, 2]
    assert select(w, "all").tolist() == [0, 1, 2]
    below = select(w, "z < 0")
    assert set(below.tolist()) <= {0, 1, 2}


def test_and_or_not_and_parentheses():
    w = molecule("H2O")
    assert select(w, "element H and z > -10").tolist() == [1, 2]
    assert select(w, "element O or element H").tolist() == [0, 1, 2]
    assert select(w, "not element H").tolist() == [0]
    assert select(w, "(element O or element H) and index 1").tolist() == [0]


def test_within_uses_the_minimum_image():
    si = bulk("Si", cubic=True)
    near = select(si, "within 2.5 of index 1")
    assert len(near) > 1 and 0 in near.tolist()
    free = si.copy(); free.pbc = False
    assert len(select(free, "within 2.5 of index 1")) <= len(near)


def test_selection_errors_say_what_is_wrong():
    w = molecule("H2O")
    with pytest.raises(SelectionError, match="分からない言葉|unknown word"):
        select(w, "resname SOL")
    with pytest.raises(SelectionError, match="原子は 3 個|only 3 atoms"):
        select(w, "index 9")
    with pytest.raises(SelectionError, match="1 から|start at 1"):
        select(w, "index 0")
    with pytest.raises(SelectionError, match="括弧|parenthesis"):
        select(w, "(element O")
    assert "3" in describe(w, None)


def test_parse_indices_is_one_based():
    assert parse_indices("1-3,5") == [0, 1, 2, 4]


def test_distance_angle_dihedral_match_ase():
    w = molecule("H2O")
    assert distance_series([w], 0, 1)[0] == pytest.approx(w.get_distance(0, 1))
    assert angle_series([w], 1, 0, 2)[0] == pytest.approx(w.get_angle(1, 0, 2))
    e = molecule("C2H6")
    assert dihedral_series([e], 2, 0, 1, 5)[0] == pytest.approx(e.get_dihedral(2, 0, 1, 5), abs=1e-6) or \
           abs(abs(dihedral_series([e], 2, 0, 1, 5)[0]) - 180.0) < 1e-6


def test_distance_crosses_the_periodic_boundary():
    from ase import Atoms

    cell = np.diag([10.0, 10.0, 10.0])
    a = Atoms("H2", positions=[[0.5, 0.5, 0.5], [9.5, 0.5, 0.5]], cell=cell, pbc=True)
    assert distance_series([a], 0, 1)[0] == pytest.approx(1.0)
    free = a.copy(); free.pbc = False
    assert distance_series([free], 0, 1)[0] == pytest.approx(9.0)


def test_series_follow_the_motion():
    w = molecule("H2O")
    frames = []
    for i in range(4):
        fr = w.copy(); fr.positions[1, 2] += 0.1 * i
        frames.append(fr)
    d = distance_series(frames, 0, 1)
    assert len(d) == 4 and np.all(np.diff(d) < 0)


# ---- RMSD, RMSF ----
def test_rmsd_removes_rotation_when_superposing():
    w = molecule("H2O")
    rotated = w.copy(); rotated.rotate(37, "z")
    assert kabsch_rmsd(w.get_positions(), rotated.get_positions())[0] == pytest.approx(0.0, abs=1e-9)
    assert kabsch_rmsd(w.get_positions(), rotated.get_positions(), superpose=False)[0] > 0.3


def test_rmsd_series_starts_at_zero():
    w = molecule("H2O")
    frames = [w.copy() for _ in range(3)]
    frames[1].positions[1, 0] += 0.2
    frames[2].positions[1, 0] += 0.4
    series = rmsd_series(frames)
    assert series[0] == pytest.approx(0.0, abs=1e-12) and series[2] > series[1] > 0


def test_rmsf_reports_per_atom_values():
    w = molecule("H2O")
    frames = []
    for i in range(5):
        fr = w.copy(); fr.positions[1, 2] += 0.05 * i
        frames.append(fr)
    got = rmsf(frames)
    assert len(got["rmsf_A"]) == 3 and got["index"] == [1, 2, 3] and got["symbol"][0] == "O"
    assert got["superposed"] and "判定" in got["note"] or "judgement" in got["note"]
    with pytest.raises(GeometryError, match="2 つ以上|at least two"):
        rmsf([w])


def test_parse_atom_list_checks_the_count():
    assert parse_atom_list("1,2,3", 3) == [0, 1, 2]
    with pytest.raises(GeometryError, match="4 個|give 4"):
        parse_atom_list("1,2", 4)


def test_selection_rewritten_for_vmd_and_ovito():
    from adit.analysis.select import to_ovito, to_vmd

    assert to_vmd("element O") == "name O" and to_ovito("element O") == 'ParticleType == "O"'
    assert to_vmd("index 1-10,12") == "index 0 to 9 11"
    assert to_ovito("index 1-10,12") == "((ParticleIndex >= 0 && ParticleIndex <= 9) || ParticleIndex == 11)"
    assert to_vmd("z < 10 and O H") == "(z < 10 and name O H)"
    assert to_ovito("z < 10 and O H") == '(Position.Z < 10 && (ParticleType == "O" || ParticleType == "H"))'
    assert to_vmd("not (element O or index 3)") == "not (name O or index 2)"
    assert to_ovito("not element O") == '(ParticleType == "O") == 0'
    assert to_vmd("within 5 of element O") == "(within 5 of name O)"
    with pytest.raises(SelectionError):
        to_ovito("within 5 of element O")
    with pytest.raises(SelectionError):
        to_vmd("element O extra")
