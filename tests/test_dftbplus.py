import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.codes.base import GenerationError
from adit.codes.dftbplus import DftbPlusGenerator
from adit.codes.sk_sets import SKSet
from adit.spec import DftbMethod, Structure, Task
from adit.validate import validate
from tests.conftest import REAL_SK_ROOT, cfg_for, make_fake_skset, water_spec

REPO = Path(__file__).resolve().parent.parent
DFTB_EXE = os.environ.get("ADIT_DFTB_EXE") or shutil.which("dftb+") or ""
HAVE_MIO = (REAL_SK_ROOT / "mio-1-1").is_dir()
HAVE_DFTB = Path(DFTB_EXE).is_file()

gen = DftbPlusGenerator()


def write_project(spec, skset, out: Path) -> Path:
    out.mkdir()
    for name, text in gen.generate(spec, skset).items():
        (out / name).write_text(text, encoding="utf-8")
    for rel, src in gen.files_to_copy(spec, skset).items():
        (out / rel).parent.mkdir(exist_ok=True)
        shutil.copy2(src, out / rel)
    return out


def test_copy_only_required_pairs_plus_docs(sk_root):
    skset = SKSet.from_dir(sk_root / "fake-1-0")
    files = gen.files_to_copy(water_spec(), skset)
    assert sorted(files) == sorted(["skf/H-H.skf", "skf/H-O.skf", "skf/O-H.skf", "skf/O-O.skf", "skf/LICENSE", "skf/README"])


def test_refuse_copy_without_docs(sk_root):
    skset = SKSet.from_dir(sk_root / "nodocs-0-1")
    with pytest.raises(GenerationError):
        gen.files_to_copy(water_spec(method=DftbMethod(sk_set="nodocs-0-1")), skset)


def test_hsd_matches_tutorial_keywords(sk_root):
    skset = SKSet.from_dir(sk_root / "fake-1-0")
    hsd = gen.dftb_in_hsd(water_spec(), skset)
    for key in ['<<< "geometry.gen"', "Driver = GeometryOptimization", "Optimizer = Rational {}", "MovedAtoms = 1:-1",
                "MaxSteps = 100", 'OutputPrefix = "geom.out"', "GradElem [eV/AA] = 0.00514221", "Scc = Yes",
                "SlaterKosterFiles = Type2FileNames", 'Prefix = "skf/"', 'O = "p"', 'H = "s"',
                "PrintForces = Yes", "ParserVersion = 14"]:
        assert key in hsd, key
    for absent in ["Charge", "SpinPolarisation", "Filling", "ThirdOrder", "Dispersion"]:
        assert absent not in hsd, absent


def test_hsd_single_point_and_options(sk_root):
    skset = SKSet.from_dir(sk_root / "fake-1-0")
    spec = water_spec(task=Task(type="single_point"),
                      method=DftbMethod(sk_set="fake-1-0", scc=False, filling_temperature=300, dispersion="lennard-jones"))
    s = water_spec().structure
    spec = spec.model_copy(update={"structure": Structure(**{**s.model_dump(), "charge": -1})})
    hsd = gen.dftb_in_hsd(spec, skset)
    assert "Driver {}" in hsd and "Scc = No" in hsd and "SccTolerance" not in hsd
    assert "Charge = -1" in hsd and "Temperature [K] = 300" in hsd
    assert "Dispersion = LennardJones" in hsd and "UFFParameters {}" in hsd


def test_missing_set_data_is_error(sk_root):
    skset = SKSet.from_dir(sk_root / "fake-1-0")
    with pytest.raises(GenerationError):
        gen.dftb_in_hsd(water_spec(method=DftbMethod(sk_set="fake-1-0", third_order=True)), skset)
    with pytest.raises(GenerationError):
        gen.dftb_in_hsd(water_spec(method=DftbMethod(sk_set="fake-1-0", dispersion="dftd3")), skset)
    s = water_spec().structure
    trip = Structure(**{**s.model_dump(), "multiplicity": 3})
    no_spinw = SKSet.from_dir(sk_root / "nodocs-0-1")
    with pytest.raises(GenerationError):
        gen.dftb_in_hsd(water_spec(structure=trip), no_spinw)
    assert "SpinConstants" in gen.dftb_in_hsd(water_spec(structure=trip), skset)
    locs = {e.location for e in validate(water_spec(structure=trip, method=DftbMethod(sk_set="fake-1-0", third_order=True, dispersion="dftd3")), cfg_for(sk_root))}
    assert locs == {"method.third_order", "method.d3_params"}
    locs = {e.location for e in validate(water_spec(structure=trip, method=DftbMethod(sk_set="nodocs-0-1")), cfg_for(sk_root))}
    assert "structure.multiplicity" in locs


def test_set_data_parsed_from_docs(tmp_path):
    d = make_fake_skset(tmp_path, "three-0-1", ["H", "O"])
    (d / "README").write_text("zeta = 4.00 (gamma^h function exponent; DampXHExponent in DFTB+)\n\n"
                              "List of all atomic Hubbard derivatives (atomic units):\n H = -0.1857\n O = -0.1575\n\nother\n", encoding="utf-8")
    (d / "spinw.txt").write_text("H:\n   -0.0717\n\nO:\n   -0.0352    -0.0296\n   -0.0296    -0.0278\n", encoding="utf-8")
    skset = SKSet.from_dir(d)
    assert skset.hubbard_derivs() == {"H": -0.1857, "O": -0.1575}
    assert skset.damping_exponent() == 4.0
    assert skset.spin_constants() == {"H": -0.0717, "O": -0.0278}
    s = water_spec().structure
    spec = water_spec(structure=Structure(**{**s.model_dump(), "charge": 1, "multiplicity": 2}),
                      method=DftbMethod(sk_set="three-0-1", third_order=True, dispersion="dftd3",
                                    d3_params={"s6": 1.0, "s8": 0.5883, "a1": 0.5719, "a2": 3.6017}))
    hsd = gen.dftb_in_hsd(spec, skset)
    for key in ["ThirdOrderFull = Yes", "O = -0.1575", "HCorrection = Damping", "Exponent = 4",
                "UnpairedElectrons = 1", "O = { -0.0278 }", "H = { -0.0717 }", "Dispersion = DftD3", "a2 = 3.6017"]:
        assert key in hsd, key


def test_geometry_gen_roundtrip(tmp_path):
    from adit.structure import from_file
    text = gen.geometry_gen(water_spec())
    p = tmp_path / "g.gen"; p.write_text(text, encoding="utf-8")
    a = from_file(p)
    assert a.get_chemical_symbols() == ["O", "H", "H"]
    assert text.splitlines()[0].split()[:2] == ["3", "C"]


@pytest.mark.skipif(not HAVE_MIO, reason="実物の mio-1-1 が無い")
def test_real_mio_generation(tmp_path):
    skset = SKSet.from_dir(REAL_SK_ROOT / "mio-1-1")
    out = write_project(water_spec(method=DftbMethod(sk_set="mio-1-1")), skset, tmp_path / "water")
    assert sorted(p.name for p in (out / "skf").iterdir()) == ["H-H.skf", "H-O.skf", "LICENSE", "O-H.skf", "O-O.skf", "README"]
    assert skset.spin_constants() is None or "H" in skset.spin_constants()


@pytest.mark.skipif(not (HAVE_MIO and HAVE_DFTB), reason="mio-1-1 か dftb+ が無い")
def test_real_run_matches_tutorial(tmp_path):
    skset = SKSet.from_dir(REAL_SK_ROOT / "mio-1-1")
    out = write_project(water_spec(method=DftbMethod(sk_set="mio-1-1")), skset, tmp_path / "water")
    r = subprocess.run([DFTB_EXE], cwd=out, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout[-2000:]
    assert "Geometry converged" in r.stdout
    from adit.structure import from_file
    a = from_file(out / "geom.out.gen")
    ref = from_file(REPO / "examples" / "water" / "geom.out.gen")
    assert abs(a.get_distance(0, 1) - ref.get_distance(0, 1)) < 1e-3
    assert abs(a.get_angle(1, 0, 2) - ref.get_angle(1, 0, 2)) < 0.05
    assert (out / "results.tag").is_file()


def c2_spec(**method):
    from adit.spec import AtomsData, KPoints

    atoms = AtomsData(symbols=["C", "C"], positions=[(0.0, 0.0, 0.0), (1.5, 1.5, 1.5)],
                      cell=[(3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0)], pbc=(True, True, True))
    return water_spec(structure=Structure(source="file", source_ref="c2", atoms=atoms), method=DftbMethod(sk_set="fake-1-0", **method),
                      task=Task(type="band_structure"), kpoints=KPoints(mode="mesh", mesh=(2, 2, 2)))


def test_non_scc_band_structure_does_not_hand_over_charges(sk_root):
    from adit.config import Profile

    skset = SKSet.from_dir(sk_root / "fake-1-0")
    profile = Profile(kind="direct")
    run = gen.run_command(c2_spec(scc=False), profile)
    assert "charges.bin" not in run and run.endswith("&& cd bands && dftb+ > output.log 2>&1")
    bands = gen.generate(c2_spec(scc=False), skset)["bands/dftb_in.hsd"]
    assert "ReadInitialCharges" not in bands and "Scc = No" in bands
    assert "電荷は引き継がず" in gen.readme_notes(c2_spec(scc=False), skset, {}).files[-1]
    run = gen.run_command(c2_spec(scc=True), profile)
    assert "&& cp charges.bin bands/ && cd bands && dftb+" in run
    assert "ReadInitialCharges = Yes" in gen.generate(c2_spec(scc=True), skset)["bands/dftb_in.hsd"]


def test_all_atoms_fixed_is_a_clear_error(sk_root):
    skset = SKSet.from_dir(sk_root / "fake-1-0")
    base = water_spec()
    spec = base.model_copy(update={"structure": base.structure.model_copy(update={"fixed_atoms": [0, 1, 2]})})
    with pytest.raises(GenerationError, match="MovedAtoms"):
        gen.generate(spec, skset)
    with pytest.raises(GenerationError, match="MovedAtoms"):
        gen.generate(spec.model_copy(update={"task": Task(type="molecular_dynamics")}), skset)
    partly = base.model_copy(update={"structure": base.structure.model_copy(update={"fixed_atoms": [0, 1]})})
    assert "MovedAtoms = 3" in gen.generate(partly, skset)["dftb_in.hsd"]
