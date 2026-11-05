
import numpy as np
import pytest
from ase.build import bulk

from adit.analysis.xrd import XrdError, has_pymatgen, powder_pattern, read_measured, scale_to_100, write_csv

pytestmark = pytest.mark.skipif(not has_pymatgen(), reason="pymatgen is missing")


def test_silicon_111_peak_is_at_the_known_angle():
    pattern = powder_pattern(bulk("Si", "diamond", a=5.431))
    assert pattern.two_theta[int(np.argmax(pattern.intensity))] == pytest.approx(28.44, abs=0.1)
    assert pattern.intensity.max() == pytest.approx(100.0)
    assert pattern.radiation == "CuKa" and pattern.wavelength_ang == pytest.approx(1.54184, abs=1e-4)
    assert "判定していません" in pattern.note or "not assessed" in pattern.note


def test_range_and_radiation_are_checked():
    si = bulk("Si", "diamond", a=5.431)
    with pytest.raises(XrdError, match="線源|radiation"):
        powder_pattern(si, "CuKx")
    with pytest.raises(XrdError, match="範囲|range"):
        powder_pattern(si, two_theta_range=(90.0, 5.0))


def test_a_molecule_without_a_cell_is_refused():
    from ase import Atoms

    with pytest.raises(XrdError, match="周期系|periodic"):
        powder_pattern(Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]]))


def test_measured_file_skips_headers_and_scales_to_100(tmp_path):
    p = tmp_path / "meas.csv"
    p.write_text("two_theta,intensity\n10,5\n28.4,900\n47.3,600\n", encoding="utf-8")
    x, y = read_measured(p)
    assert list(x) == [10.0, 28.4, 47.3]
    assert scale_to_100(y).max() == pytest.approx(100.0)
    (tmp_path / "empty.csv").write_text("# 見出しだけ\n", encoding="utf-8")
    with pytest.raises(XrdError, match="2 点以上|at least two"):
        read_measured(tmp_path / "empty.csv")


def test_csv_has_one_row_per_peak(tmp_path):
    pattern = powder_pattern(bulk("Si", "diamond", a=5.431))
    text = write_csv(tmp_path / "xrd.csv", pattern).read_text(encoding="utf-8").strip().splitlines()
    assert text[0] == "two_theta_deg,intensity_rel,d_ang,hkl"
    assert len(text) == len(pattern.two_theta) + 1
