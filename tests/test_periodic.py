
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.codes.dftbplus import DftbPlusGenerator
from adit.codes.sk_sets import SKSet
from adit.spec import HARTREE_PER_BOHR_IN_EV_PER_ANG, AtomsData, CalculationSpec, DftbMethod, KPoints, Runtime, Structure, Task
from adit.structure import from_file
from adit.validate import validate
from tests.conftest import REAL_SK_ROOT, cfg_for, water_spec

REPO = Path(__file__).resolve().parent.parent
DFTB_EXE = os.environ.get("ADIT_DFTB_EXE") or shutil.which("dftb+") or ""
HAVE_TIO2 = (REAL_SK_ROOT / "mio-ext" / "README").is_file() and Path(DFTB_EXE).is_file()
gen = DftbPlusGenerator()


def tio2_spec(sk_set="mio-ext", **kw) -> CalculationSpec:
    atoms = from_file(REPO / "examples" / "dftb_tio2" / "geometry.gen")
    st = Structure(source="file", source_ref="examples/dftb_tio2/geometry.gen", atoms=AtomsData.from_ase(atoms))
    d = dict(structure=st, method=DftbMethod(sk_set=sk_set), kpoints=KPoints(mode="mesh", mesh=(4, 4, 4), shift=(0.5, 0.5, 0.5)),
             task=Task(type="single_point"), runtime=Runtime(profile="local", omp_threads=4, job_name="tio2"))
    d.update(kw)
    return CalculationSpec(**d)


def test_v1_spec_json_is_migrated():
    spec = CalculationSpec.load(REPO / "examples" / "water_generated" / "spec.json")
    assert spec.version == 2
    assert abs(spec.task.force_tolerance_ev_per_ang - 1e-4 * HARTREE_PER_BOHR_IN_EV_PER_ANG) < 1e-5 * 1e-4 * HARTREE_PER_BOHR_IN_EV_PER_ANG
    assert spec.kpoints is None and spec.structure.fixed_atoms == []


def test_periodic_requires_kpoints_and_molecule_forbids(sk_root, tmp_path):
    from tests.conftest import make_fake_skset
    make_fake_skset(sk_root, "ti-0-0", ["Ti", "O"])
    cfg = cfg_for(sk_root)
    spec = tio2_spec(sk_set="ti-0-0")
    assert validate(spec, cfg) == []
    no_k = spec.model_copy(update={"kpoints": None})
    assert [e.location for e in validate(no_k, cfg)] == ["kpoints"]
    mol = water_spec(kpoints=KPoints(mode="gamma"), task=Task(type="geometry_optimization", relax_cell="shape_and_volume"))
    assert sorted(e.location for e in validate(mol, cfg)) == ["kpoints", "task.relax_cell"]
    bad_fix = water_spec(structure=Structure(**{**water_spec().structure.model_dump(), "fixed_atoms": [0, 7]}))
    assert [e.location for e in validate(bad_fix, cfg)] == ["structure.fixed_atoms"]


def test_kpoints_density_and_shift(sk_root):
    cell = [(4.0, 0, 0), (0, 4.0, 0), (0, 0, 8.0)]
    assert KPoints(mode="density", density=2.0).resolved_mesh(cell) == (4, 4, 2)  # 2π/4 * 2 = 3.14 → 4, 2π/8*2 = 1.57 → 2
    assert KPoints(mode="gamma").resolved_mesh(cell) == (1, 1, 1)
    from tests.conftest import make_fake_skset
    make_fake_skset(sk_root, "ti-0-0", ["Ti", "O"])
    spec = tio2_spec(sk_set="ti-0-0", kpoints=KPoints(mode="mesh", mesh=(4, 4, 4), shift=(0.3, 0, 0)))
    assert [e.location for e in validate(spec, cfg_for(sk_root))] == ["kpoints.shift"]


def test_hsd_periodic_blocks(sk_root):
    from tests.conftest import make_fake_skset
    make_fake_skset(sk_root, "ti-0-0", ["Ti", "O"])
    skset = SKSet.from_dir(sk_root / "ti-0-0")
    spec = tio2_spec(sk_set="ti-0-0", task=Task(type="geometry_optimization", relax_cell="shape_and_volume"),
                     structure=Structure(**{**tio2_spec().structure.model_dump(), "fixed_atoms": [0, 1]}))
    hsd = gen.dftb_in_hsd(spec, skset)
    for key in ["KPointsAndWeights = SupercellFolding", "    4 0 0", "    0.5 0.5 0.5", "LatticeOpt = Yes",
                "MovedAtoms = 3 4 5 6", "GradElem [eV/AA]"]:
        assert key in hsd, key
    assert "Isotropic" not in hsd
    vol = gen.dftb_in_hsd(spec.model_copy(update={"task": Task(type="geometry_optimization", relax_cell="volume_only")}), skset)
    assert "Isotropic = Yes" in vol
    geo = gen.geometry_gen(spec)
    assert geo.splitlines()[0].split()[1] == "S"
    assert len(geo.strip().splitlines()) == 2 + 6 + 4


@pytest.mark.skipif(not HAVE_TIO2, reason="mio-ext (mio-1-1 + tiorg-0-1) か dftb+ が無い")
def test_tio2_matches_recipe_energy(tmp_path):
    from adit.project import write_project
    out = tmp_path / "tio2"
    write_project(tio2_spec(), cfg_for(REAL_SK_ROOT), out)
    env = {**os.environ, "PATH": f"{Path(DFTB_EXE).parent}:{os.environ.get('PATH', '')}"}
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=600, env=env)
    assert r.returncode == 0, r.stderr
    log = (out / "output.log").read_text(encoding="utf-8")
    e = float([l for l in log.splitlines() if "Total Energy" in l][-1].split()[2])
    ref = float([l for l in (REPO / "examples" / "dftb_tio2" / "output.log").read_text(encoding="utf-8").splitlines() if "Total Energy" in l][-1].split()[2])
    assert abs(e - ref) < 1e-7, (e, ref)
    assert "Periodic boundaries:         Yes" in log
