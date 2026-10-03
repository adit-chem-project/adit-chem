

import pytest
from ase.build import molecule
from ase.io import write

from adit.analysis.readers import _INPUTS, detect_code, load_run
from adit.analysis.readers_generic import find_readable
from adit.codes import GENERATORS

GENERIC = ("gaussian", "gamess", "qchem", "openmx", "amber", "namd", "grrm", "dcdftbmd")


def test_every_generator_has_an_analysis_entry():
    covered = {code for _, code in _INPUTS}
    assert set(GENERATORS) <= covered, sorted(set(GENERATORS) - covered)


@pytest.mark.parametrize("name,code", [(n, c) for n, c in _INPUTS if c in GENERIC])
def test_generic_codes_are_detected_by_their_real_input_name(tmp_path, name, code):
    (tmp_path / name).write_text("dummy\n", encoding="utf-8")
    assert detect_code(tmp_path) == code


def test_dftbplus_wins_over_dcdftbmd_when_both_names_exist(tmp_path):
    (tmp_path / "dftb_in.hsd").write_text("Geometry = {}\n", encoding="utf-8")
    (tmp_path / "dftb.inp").write_text("KEYWORD\n", encoding="utf-8")
    assert detect_code(tmp_path) == "dftbplus"


def test_generic_reader_uses_an_ase_readable_trajectory(tmp_path):
    (tmp_path / "gaussian.gjf").write_text("#p HF/STO-3G\n", encoding="utf-8")
    frames = [molecule("H2O") for _ in range(3)]
    for i, fr in enumerate(frames):
        fr.positions[0, 2] += 0.01 * i
    write(tmp_path / "trajectory.xyz", frames)
    data = load_run(tmp_path)
    assert data.code == "gaussian" and len(data.frames) == 3
    assert data.frame_source == "trajectory.xyz"
    assert any("ASE" in n for n in data.notes)
    assert data.energies_ev == []
    assert any("エネルギーの推移" in n or "no energies" in n for n in data.notes)


def test_generic_reader_reads_energies_from_extxyz(tmp_path):
    from ase.calculators.singlepoint import SinglePointCalculator

    (tmp_path / "qchem.in").write_text("$molecule\n$end\n", encoding="utf-8")
    frames = []
    for i in range(4):
        fr = molecule("H2O")
        fr.calc = SinglePointCalculator(fr, energy=-76.0 + 0.01 * i)
        frames.append(fr)
    write(tmp_path / "traj.extxyz", frames)
    data = load_run(tmp_path)
    assert len(data.energies_ev) == 4 and data.energies_ev[0] == pytest.approx(-76.0)


def test_generic_reader_says_what_is_missing(tmp_path):
    (tmp_path / "namd.conf").write_text("numsteps 10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="ASE が読める|ASE can read"):
        load_run(tmp_path)


def test_find_readable_prefers_the_trajectory(tmp_path):
    write(tmp_path / "final.xyz", molecule("H2O"))
    write(tmp_path / "traj.extxyz", [molecule("H2O"), molecule("H2O")])
    frames, kind, path = find_readable(tmp_path)
    assert path.name == "traj.extxyz" and len(frames) == 2
