
from pathlib import Path

import pytest

from adit.analysis import trajectory_ext as ext

HAS_LIB = bool(ext.available())
GROMACS = Path("examples/gromacs_spce_nvt_generated")


def test_external_formats_are_recognised():
    assert ext.needs_external(Path("a.xtc")) and ext.needs_external(Path("b.trr")) and ext.needs_external(Path("c.dcd"))
    assert not ext.needs_external(Path("d.xyz")) and not ext.needs_external(Path("e.extxyz"))


def test_missing_library_says_how_to_convert(monkeypatch):
    monkeypatch.setattr(ext, "available", lambda: [])
    with pytest.raises(ext.ExternalTrajectoryError) as got:
        ext.read_external(GROMACS / "adit.xtc")
    message = str(got.value)
    assert "MDAnalysis" in message and "gmx trjconv" in message


def test_gromacs_reader_does_not_pretend_to_have_the_trajectory(monkeypatch):
    from adit.analysis import readers_extra

    monkeypatch.setattr(ext, "available", lambda: [])
    data = readers_extra.read_gromacs(GROMACS)
    assert len(data.frames) == 1 and data.frame_source.endswith(".gro")
    assert any("読めません" in n or "cannot be read" in n for n in data.notes)


@pytest.mark.skipif(not HAS_LIB, reason="MDAnalysis / MDTraj / chemfiles がない")
def test_reads_the_real_xtc_and_matches_gromacs():
    frames, lib = ext.read_external(GROMACS / "adit.xtc", ext.find_topology(GROMACS, GROMACS / "adit.xtc"))
    assert len(frames) == 11 and len(frames[0]) == 648
    assert frames[0].get_chemical_symbols()[:3] == ["O", "H", "H"]
    assert frames[0].cell.lengths()[0] == pytest.approx(18.6206, abs=1e-3)
    assert frames[0].positions[0][0] == pytest.approx(2.3, abs=1e-3)


@pytest.mark.skipif(not HAS_LIB, reason="MDAnalysis / MDTraj / chemfiles がない")
def test_gromacs_analysis_uses_the_whole_trajectory():
    from adit.analysis.readers import load_run

    data = load_run(GROMACS)
    assert len(data.frames) == 11 and data.frame_source == "adit.xtc"
    assert data.frame_dt_fs == pytest.approx(50.0)
