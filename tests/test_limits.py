
import sys
import time

import pytest
import numpy as np
from ase.build import bulk

from tests.conftest import cfg_for, water_spec
from adit.spec import AtomsData, BandSettings, DftbMethod, KPoints, Structure, Task
from adit.validate import MAX_BAND_POINTS, MAX_KPOINTS, _close_pairs, validate


def _bulk_spec(**kw):
    si = bulk("Si", "diamond", a=5.43)
    base = water_spec(method=DftbMethod(sk_set="fake-1-0"))
    return base.model_copy(update={"structure": Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(si)),
                                   "kpoints": KPoints(mode="mesh", mesh=(2, 2, 2)), **kw})


def _locs(spec):
    return [e.location for e in validate(spec, cfg_for(None))]


def test_band_points_over_limit_is_rejected_quickly():
    t0 = time.time()
    locs = _locs(_bulk_spec(task=Task(type="band_structure", bands=BandSettings(npoints=10_000_000))))
    assert "task.bands.npoints" in locs and time.time() - t0 < 5
    assert "task.bands.npoints" not in _locs(_bulk_spec(task=Task(type="band_structure", bands=BandSettings(npoints=MAX_BAND_POINTS))))


def test_kpoint_total_over_limit_is_rejected_quickly():
    t0 = time.time()
    assert "kpoints.density" in _locs(_bulk_spec(kpoints=KPoints(mode="density", density=1e6)))
    assert "kpoints.mesh" in _locs(_bulk_spec(kpoints=KPoints(mode="mesh", mesh=(1000, 1000, 1000))))
    assert time.time() - t0 < 5
    assert not any(l.startswith("kpoints") for l in _locs(_bulk_spec(kpoints=KPoints(mode="mesh", mesh=(100, 100, 100)))))
    assert MAX_KPOINTS == 1_000_000


def test_close_pairs_molecule_and_periodic_boundary():
    mol = Structure(source="preset", source_ref="x", atoms=AtomsData(symbols=["H", "H", "O"], positions=[(0, 0, 0), (0.3, 0, 0), (3, 0, 0)]))
    assert [(i, j) for i, j, _ in _close_pairs(mol)] == [(0, 1)]
    per = Structure(source="bulk", source_ref="x", atoms=AtomsData(symbols=["H", "H"], positions=[(0.1, 5, 5), (9.9, 5, 5)],
                                                                   cell=[(10, 0, 0), (0, 10, 0), (0, 0, 10)], pbc=(True, True, True)))
    (i, j, d), = _close_pairs(per)
    assert (i, j) == (0, 1) and abs(d - 0.2) < 1e-9
    si = bulk("Si", "diamond", a=5.43)
    assert _close_pairs(Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(si))) == []


def test_overlap_check_scales_for_thousands_of_atoms():
    resource = pytest.importorskip("resource", reason="Windows has no resource module")
    c = bulk("C", "diamond", a=3.567, cubic=True).repeat((8, 8, 8))
    st = Structure(source="bulk", source_ref="C", atoms=AtomsData.from_ase(c))
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.time()
    assert _close_pairs(st) == []
    unit = 1024 ** 2 if sys.platform == "darwin" else 1024
    grew_mb = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - before) / unit
    assert time.time() - t0 < 10 and grew_mb < 200, (time.time() - t0, grew_mb)
    assert len(c) == 4096 and np.isfinite(grew_mb)
