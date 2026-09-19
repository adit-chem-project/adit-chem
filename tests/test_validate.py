from pathlib import Path

import pytest

from adit.codes.sk_sets import SKSet, discover_sets
from adit.spec import AtomsData, DftbMethod, Runtime, Structure, Task
from adit.validate import validate
from tests.conftest import REAL_SK_ROOT, cfg_for, water_spec


def locations(errs):
    return sorted({e.location for e in errs})


def test_clean_spec_has_no_errors(sk_root):
    assert validate(water_spec(), cfg_for(sk_root)) == []


def test_missing_element_in_set(sk_root):
    spec = water_spec(
        structure=Structure(
            source="preset", source_ref="x",
            atoms=AtomsData(symbols=["Ti", "O", "O"], positions=[(0, 0, 0), (0, 0, 1.6), (0, 0, -1.6)]),
        )
    )
    errs = validate(spec, cfg_for(sk_root))
    assert locations(errs) == ["method.sk_set"]
    assert "Ti" in errs[0].message


def test_unknown_set(sk_root):
    errs = validate(water_spec(method=DftbMethod(sk_set="nope-9-9")), cfg_for(sk_root))
    assert locations(errs) == ["method.sk_set"]


def test_missing_sk_root(tmp_path):
    errs = validate(water_spec(), cfg_for(tmp_path / "absent"))
    assert locations(errs) == ["method.sk_set"]
    assert validate(water_spec(), cfg_for(None))[0].location == "method.sk_set"


def test_set_without_license_or_readme(sk_root):
    errs = validate(water_spec(method=DftbMethod(sk_set="nodocs-0-1")), cfg_for(sk_root))
    msgs = " ".join(e.message for e in errs)
    assert "LICENSE" in msgs and "README" in msgs


def test_overlapping_atoms(sk_root):
    spec = water_spec(
        structure=Structure(
            source="preset", source_ref="x",
            atoms=AtomsData(symbols=["O", "H", "H"], positions=[(0, 0, 0), (0, 0, 0.3), (0, 0, -0.9)]),
        )
    )
    errs = validate(spec, cfg_for(sk_root))
    assert locations(errs) == ["structure.atoms"]
    assert "原子 1 (O) と 2 (H)" in errs[0].message


def test_empty_structure(sk_root):
    spec = water_spec(structure=Structure(source="file", source_ref="x", atoms=AtomsData(symbols=[], positions=[])))
    assert locations(validate(spec, cfg_for(sk_root))) == ["structure.atoms"]


def test_negative_steps_and_bad_numbers(sk_root):
    spec = water_spec(
        task=Task(type="geometry_optimization", max_steps=-1, force_tolerance_ev_per_ang=0.0),
        method=DftbMethod(sk_set="fake-1-0", scc_tolerance=0, max_scc_iterations=0, filling_temperature=-1),
        runtime=Runtime(profile="local", ncpus=0, omp_threads=0, walltime="1h", job_name="  "),
    )
    assert locations(validate(spec, cfg_for(sk_root))) == sorted([
        "task.max_steps", "task.force_tolerance_ev_per_ang",
        "method.scc_tolerance", "method.max_scc_iterations", "method.filling_temperature",
        "runtime.ncpus", "runtime.omp_threads", "runtime.walltime", "runtime.job_name",
    ])


def test_max_steps_ignored_for_single_point(sk_root):
    spec = water_spec(task=Task(type="single_point", max_steps=-5))
    assert validate(spec, cfg_for(sk_root)) == []


def test_multiplicity_parity(sk_root):
    s = water_spec().structure
    odd = Structure(**{**s.model_dump(), "multiplicity": 2})
    assert locations(validate(water_spec(structure=odd), cfg_for(sk_root))) == ["structure.multiplicity"]
    cation_doublet = Structure(**{**s.model_dump(), "charge": 1, "multiplicity": 2})
    assert validate(water_spec(structure=cation_doublet), cfg_for(sk_root)) == []


def test_output_dir_parent_must_exist(sk_root, tmp_path):
    assert validate(water_spec(), cfg_for(sk_root), output_dir=tmp_path / "new") == []
    errs = validate(water_spec(), cfg_for(sk_root), output_dir=tmp_path / "no" / "such" / "dir")
    assert locations(errs) == ["output_dir"]
    (tmp_path / "file").write_text("x", encoding="utf-8")
    errs = validate(water_spec(), cfg_for(sk_root), output_dir=tmp_path / "file")
    assert locations(errs) == ["output_dir"] and "file" in errs[0].message


def test_skset_reads_shells(sk_root):
    s = SKSet.from_dir(sk_root / "fake-1-0")
    assert s.elements == ["C", "H", "N", "O"]
    assert s.max_angular_momentum("H") == "s" and s.max_angular_momentum("O") == "p"
    assert sorted(s.required_files(["O", "H"])) == [("H", "H"), ("H", "O"), ("O", "H"), ("O", "O")]
    assert set(discover_sets(sk_root)) == {"fake-1-0", "nodocs-0-1"}


@pytest.mark.skipif(not (REAL_SK_ROOT / "mio-1-1").is_dir(), reason="実物の mio-1-1 が無い")
def test_real_mio_set():
    s = SKSet.from_dir(REAL_SK_ROOT / "mio-1-1")
    assert s.elements == ["C", "H", "N", "O", "P", "S"]
    assert s.max_angular_momentum("H") == "s"
    assert s.max_angular_momentum("O") == "p"
    assert s.max_angular_momentum("S") == "d"
    assert set(s.doc_files()) == {"LICENSE", "README"}
    assert validate(water_spec(method=DftbMethod(sk_set="mio-1-1")), cfg_for(REAL_SK_ROOT)) == []


def test_overlap_across_periodic_boundary(sk_root):
    from tests.conftest import make_fake_skset
    make_fake_skset(sk_root, "si-0-0", ["Si"])
    from adit.spec import KPoints
    atoms = AtomsData(symbols=["Si", "Si"], positions=[(0.1, 0, 0), (4.9, 0, 0)], cell=[(5, 0, 0), (0, 5, 0), (0, 0, 5)], pbc=(True, True, True))
    spec = water_spec(structure=Structure(source="file", source_ref="x", atoms=atoms), method=DftbMethod(sk_set="si-0-0"),
                      kpoints=KPoints(mode="gamma"), task=Task(type="single_point"))
    errs = validate(spec, cfg_for(sk_root))
    assert [e.location for e in errs] == ["structure.atoms"] and "0.200" in errs[0].message
