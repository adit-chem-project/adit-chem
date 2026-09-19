"""Playing trajectories and vibrational modes in the analysis tab; reading normal modes from the code outputs."""

import gzip
import os
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from adit.analysis.modes import Modes, modes_from_hessian, read_g98_modes, read_modes, read_orca_modes
from adit.analysis.readers import frequencies_from_hessian
from adit.analysis.trajectory import Trajectory

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"
DATA = REPO / "tests" / "data"


def test_dftb_hessian_gives_the_same_frequencies_and_unit_scaled_vectors():
    water = read(EX / "water_generated" / "geom.out.gen")
    h = np.array((DATA / "water_hessian.out").read_text(encoding="utf-8").split(), dtype=float)
    freqs, vec = modes_from_hessian(h, water)
    assert freqs == pytest.approx(frequencies_from_hessian(h, water), abs=1e-3)   # cm-1; eigh vs eigvalsh differ only in the ~0 modes
    assert vec.shape == (9, 3, 3)
    assert np.linalg.norm(vec, axis=2).max(axis=1) == pytest.approx(np.ones(9))     # the largest atom moves 1
    # the mode vectors of the mass-weighted Hessian stay eigenvectors after the 1/sqrt(m) un-weighting:
    # they must satisfy H v = lambda M v (generalized eigenproblem) for every mode
    n = 9
    hm = 0.5 * (h[: n * n].reshape(n, n) + h[: n * n].reshape(n, n).T)
    m = np.repeat(water.get_masses() * 1822.888486209, 3)
    for k in range(n):
        v = vec[k].reshape(-1)
        lam = np.sign(freqs[k]) * (freqs[k] / 219474.6313705) ** 2
        assert hm @ v == pytest.approx(lam * m * v, abs=1e-8)


def test_mode_displacement_is_sinusoidal_with_the_given_amplitude():
    m = read_modes(EX / "xtb_vib_water_generated")
    assert m is not None and len(m) == 3 and m.source == "g98.out"
    assert m.frequencies_cm1 == pytest.approx([819.4274, 1100.232, 1656.7098])
    eq = m.atoms.get_positions()
    assert m.displaced(2, 0.5, 0.0) == pytest.approx(eq)
    d = m.displaced(2, 0.5, np.pi / 2) - eq
    assert np.linalg.norm(d, axis=1).max() == pytest.approx(0.5)
    assert np.linalg.norm(d, axis=1)[0] < np.linalg.norm(d, axis=1)[1]           # the H atoms move more than O in the stretch
    assert m.displaced(2, 0.5, 3 * np.pi / 2) - eq == pytest.approx(-d)
    # the g98 rows are read in the right column order: mode 3 (1656.7) on H2 is (-0.00, 0.10, -0.70) before scaling
    v = m.vectors[2][1] / np.linalg.norm(m.vectors[2][1])
    assert v == pytest.approx(np.array([0.0, 0.10, -0.70]) / np.hypot(0.10, 0.70), abs=1e-6)


def test_gaussian_and_orca_normal_modes(tmp_path):
    g = tmp_path / "gaussian.log"; g.write_bytes(gzip.decompress((DATA / "gaussian_dvb_raman.out.gz").read_bytes()))
    m = read_g98_modes(g)
    assert m is not None and len(m) == 54 and len(m.atoms) == 20
    assert m.frequencies_cm1[:2] == pytest.approx([53.1117, 84.6059])
    o = tmp_path / "output.log"; o.write_bytes(gzip.decompress((DATA / "orca_dvb_raman.out.gz").read_bytes()))
    m = read_orca_modes(o)
    assert m is not None and len(m) == 54 and len(m.atoms) == 20             # the six zero modes are dropped
    assert m.frequencies_cm1[0] == pytest.approx(36.86)
    assert np.linalg.norm(m.vectors, axis=2).max(axis=1) == pytest.approx(np.ones(54))


def test_unsupported_and_missing_modes_return_none():
    assert read_modes(EX / "dftb_md_water_generated") is None            # MD run without hessian.out
    assert read_modes(EX / "lammps_cu_nvt_generated") is None            # code without mode display


pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_playback_of_the_md_example_follows_the_frames(app):
    from adit.gui.playback import PlaybackPanel

    p = PlaybackPanel()
    assert p.load_run(EX / "dftb_md_water_generated")
    assert p.count() == 21 and p.slider.maximum() == 20 and not p.playing_modes()
    traj = Trajectory(EX / "dftb_md_water_generated" / "geo_end.xyz", "xyz")
    p.go(5)
    assert p.viewer.atoms().get_positions() == pytest.approx(traj[5].get_positions())
    assert p.position.text().startswith("フレーム 6 / 21") and "t = 25.0 fs" in p.position.text()   # 0.5 fs x every 10 steps
    assert p.series._markers and p.series._markers[0].get_xdata()[0] == pytest.approx(25.0)
    assert len(p.series._axes) == 2                                        # energy and temperature
    p.btn_next.click(); assert p.frame == 6
    p.btn_last.click(); assert p.frame == 20
    p.btn_next.click(); assert p.frame == 0                                # wraps around
    p.play(); assert p.is_playing(); p._tick(); assert p.frame == 1; p.pause(); assert not p.is_playing()
    p._on_series_click(101.0); assert p.frame == 20                        # clicking the plot jumps to the nearest frame
    p.slider.setValue(3); assert p.frame == 3
    assert p.grab().width() > 0


def test_playback_stride_and_skip_thin_the_frames(app):
    from adit.gui.playback import PlaybackPanel

    p = PlaybackPanel()
    assert p.load_run(EX / "dftb_md_water_generated", skip=1, stride=5)
    assert p.count() == 4                                                  # frames 1, 6, 11, 16
    p.go(1)
    assert "t = 30.0 fs" in p.position.text()
    p2 = PlaybackPanel()
    assert p2.load_run(EX / "dftb_md_water_generated", memory_mb=0.0005)   # 21 frames x 3 atoms need ~0.0014 MB: thinned automatically
    assert 1 < p2.count() < 21 and "間引いて" in p2.note.text()


def test_playback_of_vibrational_modes(app):
    from adit.gui.playback import PHASES_PER_PERIOD, PlaybackPanel

    p = PlaybackPanel()
    assert p.load_run(EX / "xtb_vib_water_generated")
    assert p.playing_modes() and p.count() == PHASES_PER_PERIOD and p.mode.count() == 3
    assert not p.what.isVisibleTo(p)                                        # only modes here, nothing to choose
    p.mode.setCurrentIndex(2); p.amplitude.setValue(0.3)
    eq = p._modes.atoms.get_positions()
    p.go(PHASES_PER_PERIOD // 4)
    d = p.viewer.atoms().get_positions() - eq
    assert np.linalg.norm(d, axis=1).max() == pytest.approx(0.3)
    assert "1656.7 cm" in p.position.text()
    p._on_series_click(820.0); assert p.mode.currentIndex() == 0
    p.go(0); assert p.viewer.atoms().get_positions() == pytest.approx(eq)


def test_analysis_panel_shows_the_player_after_a_run(app, tmp_path):
    import shutil

    from adit.gui.panels.analysis_panel import AnalysisPanel

    md = tmp_path / "md"; shutil.copytree(EX / "dftb_md_water_generated", md)      # the analysis writes next to the run
    opt = tmp_path / "opt"; shutil.copytree(EX / "water_generated", opt)
    panel = AnalysisPanel()
    panel.set_run_dir(md, ["O", "H"])
    panel.cb_rdf.setChecked(False); panel.cb_msd.setChecked(False)
    res = panel.run()
    assert res is not None
    assert panel.play_box.isVisibleTo(panel) and panel.playback.count() == 21
    panel.set_run_dir(opt, ["O", "H"])
    panel.run()
    assert not panel.play_box.isVisibleTo(panel)                           # one frame, no hessian: nothing to play, so no box
    assert panel.playback.count() == 0 and "フレームが 1 つ" in panel.playback.note.text()
