

import numpy as np
import pytest
from ase.build import bulk

from tests.conftest import cfg_for, water_spec
from adit.analysis import AnalysisOptions, run_analysis
from adit.spec import AtomsData, EspressoMethod, KPoints, MDSettings, Structure, Task, XtbMethod


def si_spec(**kw):
    s = bulk("Si", "diamond", a=5.43)
    return water_spec(structure=Structure(source="preset", source_ref="Si", atoms=AtomsData.from_ase(s)),
                      method=EspressoMethod(ecutwfc=30), kpoints=KPoints(), task=Task(type="single_point"), **kw)


def _md_run(tmp_path):
    from adit.project import write_project

    cfg = cfg_for(None)
    d = tmp_path / "md"
    spec = water_spec(task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NVE", steps=4, dump_interval=1)))
    cfg = cfg_for(tmp_path / "slakos")
    from tests.conftest import make_fake_skset
    make_fake_skset(tmp_path / "slakos", "fake-1-0", ["H", "C", "N", "O"])
    write_project(spec, cfg, d)
    (d / "geo_end.xyz").write_text("".join(
        f"3\nMD iter: {k}, time: {k:.2f}\nO 0.0 0.0 {0.01 * k:.4f}\nH 0.0 0.0 0.9\nH 0.0 0.9 0.0\n" for k in range(5)), encoding="utf-8")
    (d / "md.out").write_text("".join(f"{k:>8d}{-4.0 - 0.001 * k:20.10f}{-4.0:20.10f}{-4.0:20.10f}{300.0 + k:20.10f}\n" for k in range(5)), encoding="utf-8")
    (d / "output.log").write_text("Total Energy:  -4.0 H  -108.8 eV\n", encoding="utf-8")
    return d


@pytest.mark.parametrize("kw, word", [
    ({"skip_frames": -3}, "skip"),
    ({"stride": 0}, "stride"),
    ({"rdf_rmax": 0.0}, "動径分布関数"),
    ({"rdf_rmax": -1.0}, "動径分布関数"),
    ({"dos_sigma": 0.0}, "状態密度"),
    ({"zdens_bin_ang": 0.0}, "z 方向"),
    ({"zdens_bin_ang": -1.0}, "z 方向"),
    ({"memory_budget_mb": 0}, "MSD"),
    ({"memory_budget_mb": -5}, "MSD"),
    ({"bands_window_ev": 0.0}, "バンド"),
    ({"rdf_rmax": float("nan")}, "動径分布関数"),
])
def test_analysis_options_are_checked(tmp_path, kw, word):
    d = _md_run(tmp_path)
    with pytest.raises(ValueError) as ex:
        run_analysis(d, AnalysisOptions(**kw))
    assert word in str(ex.value)


def test_negative_skip_used_to_take_the_last_frames(tmp_path):
    d = _md_run(tmp_path)
    assert run_analysis(d, AnalysisOptions()).tables["trajectory"]["n_frames_used"] == 5
    with pytest.raises(ValueError):
        run_analysis(d, AnalysisOptions(skip_frames=-3))


def test_compare_conditions_csv_has_values(tmp_path):
    from tests.conftest import make_fake_skset
    from adit.analysis.compare import analyze_compare
    from adit.compare_sets import write_compare_set

    make_fake_skset(tmp_path / "slakos", "fake-1-0", ["H", "C", "N", "O"])
    cfg = cfg_for(tmp_path / "slakos")
    spec = water_spec(method=XtbMethod(), task=Task(type="single_point"))
    out = tmp_path / "cmp"
    dirs = write_compare_set(spec, cfg, out, {"kind": "solvation", "structure": "base",
                                              "solvent": {"method": {"solvation": "alpb", "solvent": "water"}}}, tmp_path)
    for p, e in zip(dirs, (-5.0, -5.1)):
        (p / "output.log").write_text(f"          | TOTAL ENERGY  {e:.12f} Eh   |\n", encoding="utf-8")
    analyze_compare(out)
    rows = [l.split(",") for l in (out / "compare_conditions.csv").read_text(encoding="utf-8").splitlines()]
    head = rows[0]
    assert head[:3] == ["setting", "solvent", "gas"]
    table = {r[0]: r[1:3] for r in rows[1:]}
    assert table["task.type"] == ["single_point", "single_point"]
    assert table["method.solvation"] == ["alpb", "none"]
    assert all(v != ["", ""] for k, v in table.items() if k.startswith("method.code"))


def test_elastic_rejects_non_finite_and_indistinguishable_strains():
    from adit.elastic_setup import ElasticError, plan

    s = si_spec()
    for bad in ([float("nan")], [float("inf")], [0.01, float("nan")]):
        with pytest.raises(ElasticError):
            plan(s, bad, [1])
    with pytest.raises(ElasticError) as ex:
        plan(s, [1e-9, 1e-7], [1])
    assert "名前" in str(ex.value) or "name" in str(ex.value)
    assert len(plan(s, [-0.01, 0.01], [1])) == 3


def test_phonon_rejects_non_finite_displacement_and_width():
    from adit.phonon_setup import PhononError, _check_dos, plan

    s = si_spec()
    for bad in (float("nan"), float("inf"), 0.0, -0.1):
        with pytest.raises(PhononError):
            plan(s, (1, 1, 1), distance=bad, backend="ase")
    for bad in (float("nan"), 0.0, -1.0):
        with pytest.raises(PhononError):
            _check_dos(s.atoms, "ase", (2, 2, 2), bad)


def test_compare_set_rejects_duplicate_names_and_non_finite_nu(tmp_path):
    from adit.compare_sets import CompareSetError, plan_set

    s = si_spec()
    dup = {"kind": "reaction", "members": [{"name": "a", "structure": "base", "nu": 1},
                                           {"name": "a", "structure": "base", "nu": -1}]}
    with pytest.raises(CompareSetError) as ex:
        plan_set(s, dup, tmp_path)
    assert "name" in str(ex.value)
    nan = {"kind": "reaction", "members": [{"name": "a", "structure": "base", "nu": float("nan")},
                                           {"name": "b", "structure": "base", "nu": -1}]}
    with pytest.raises(CompareSetError):
        plan_set(s, nan, tmp_path)


def test_compare_set_error_has_no_raw_pydantic_text(tmp_path):
    from adit.compare_sets import CompareSetError, plan_set

    bad = {"kind": "reaction", "members": [{"name": "a", "structure": "base", "nu": 1, "charge": "x"},
                                           {"name": "b", "structure": "base", "nu": -1}]}
    with pytest.raises(CompareSetError) as ex:
        plan_set(si_spec(), bad, tmp_path)
    t = str(ex.value)
    assert "全電荷" in t
    assert not any(w in t for w in ("pydantic.dev", "input_value", "validation error", "[type="))


@pytest.mark.parametrize("name, text", [
    ("phonons.json", '{"backend": "ase", "dim": null, "distance_ang": null, "dirs": []}'),
    ("phonons.json", "{ broken"),
    ("elastic.json", '{"components": [1], "strains": [0.01], "runs": [], "c_gpa": [[null,null,null,null,null,null]]}'),
    ("compare.json", '{"reactions": [{"terms": [{"dir": "a"}]}]}'),
    ("compare.json", '{"reactions": [{"name": "r", "terms": [{"dir": "a", "nu": NaN}]}]}'),
    ("compare.json", "{ broken"),
])
def test_broken_collection_files_give_a_reason(tmp_path, name, text):
    d = tmp_path / name.split(".")[0]
    d.mkdir()
    (d / name).write_text(text, encoding="utf-8")
    res = run_analysis(d, AnalysisOptions())
    assert res.summary_text()


def test_load_compare_file_reports_missing_nu(tmp_path):
    from adit.analysis.compare import CompareError, load_compare_file

    p = tmp_path / "compare.json"
    p.write_text('{"reactions": [{"terms": [{"dir": "a"}]}]}', encoding="utf-8")
    with pytest.raises(CompareError):
        load_compare_file(p)
    p.write_text("{ broken", encoding="utf-8")
    with pytest.raises(CompareError):
        load_compare_file(p)


def test_stages_keep_the_continuation_handoff(tmp_path):
    from tests.conftest import make_fake_skset
    from adit.continuation import continue_from
    from adit.project import write_project
    from adit.stages import Stage, write_stages

    make_fake_skset(tmp_path / "slakos", "fake-1-0", ["H", "C", "N", "O"])
    cfg = cfg_for(tmp_path / "slakos")
    md = Task(type="molecular_dynamics", md=MDSettings(ensemble="NVT", steps=10, dump_interval=1))
    prev = tmp_path / "prev"
    write_project(water_spec(method=XtbMethod(), task=md), cfg, prev)
    (prev / "xtb.trj").write_text("3\n energy: -5.0\nO 0.0 -1.0 0.0\nH 0.0 0.0 0.78\nH 0.0 0.0 -0.78\n", encoding="utf-8")
    (prev / "mdrestart").write_text("fake restart\n", encoding="utf-8")
    spec = continue_from(prev)
    assert spec.handoff is not None and spec.handoff.files == {"mdrestart": "mdrestart"}
    dirs = write_stages(spec, cfg, tmp_path / "st",
                        [Stage("nvt", {"task": {"type": "molecular_dynamics"}}), Stage("prod", {"task": {"type": "molecular_dynamics"}})])
    assert (dirs[0] / "mdrestart").is_file()
    assert "restart=true" in (dirs[0] / "xtb.inp").read_text(encoding="utf-8")
    from adit.spec import CalculationSpec
    s2 = CalculationSpec.load(dirs[1] / "spec.json")
    assert s2.handoff is not None and s2.handoff.at_run and s2.handoff.previous_dir == "../stage_01_nvt"


def test_stress_and_force_units(tmp_path):
    from adit.outputs import HARTREE_BOHR3_TO_GPA, HARTREE_BOHR_TO_EV_ANG, RY_BOHR_TO_EV_ANG, OutputError, read_forces, read_stress

    c = bulk("C", "diamond", a=3.567)
    st = Structure(source="preset", source_ref="C", atoms=AtomsData.from_ase(c))
    d = tmp_path / "dftb"
    d.mkdir()
    water_spec(structure=st, kpoints=KPoints(), task=Task(type="single_point")).save(d / "spec.json")
    s_au = 1.0e-4
    body = (f"\n Total Forces\n    1     0.0100000000    0.0000000000    0.0000000000\n"
            f"    2    -0.0100000000    0.0000000000    0.0000000000\n\n Total stress tensor\n"
            f"   {s_au:.12E}   0.0  0.0\n   0.0   {s_au:.12E}  0.0\n   0.0   0.0   {s_au:.12E}\n\n")
    (d / "detailed.out").write_text(body + f"Pressure:                     {s_au:.6E} au    1.0 Pa\n", encoding="utf-8")
    assert np.allclose(np.diag(read_stress(d)), -s_au * HARTREE_BOHR3_TO_GPA)
    assert np.allclose(read_forces(d)[0], [0.01 * HARTREE_BOHR_TO_EV_ANG, 0, 0])
    (d / "detailed.out").write_text(body + f"Pressure:                     {-s_au:.6E} au    1.0 Pa\n", encoding="utf-8")
    assert np.allclose(np.diag(read_stress(d)), s_au * HARTREE_BOHR3_TO_GPA)
    (d / "detailed.out").write_text(body + "Pressure:                     5.000000E-04 au    1.0 Pa\n", encoding="utf-8")
    with pytest.raises(OutputError):
        read_stress(d)
    q = tmp_path / "qe"
    q.mkdir()
    si_spec().model_copy(update={"structure": st}).save(q / "spec.json")
    (q / "output.log").write_text(
        "\n     total   stress  (Ry/bohr**3)                   (kbar)     P=       20.11\n"
        "   0.00013671   0.00000000   0.00000000           20.11      0.00      0.00\n"
        "   0.00000000   0.00013671   0.00000000            0.00     20.11      0.00\n"
        "   0.00000000   0.00000000   0.00013671            0.00      0.00     20.11\n\n"
        "     Forces acting on atoms (cartesian axes, Ry/au):\n\n"
        "     atom    1 type  1   force =     0.01000000    0.00000000    0.00000000\n"
        "     atom    2 type  1   force =    -0.01000000    0.00000000    0.00000000\n", encoding="utf-8")
    assert np.allclose(np.diag(read_stress(q)), -2.011)
    assert np.allclose(read_forces(q)[0], [0.01 * RY_BOHR_TO_EV_ANG, 0, 0])


def test_qe_neb_readme_matches_the_written_nstep(tmp_path):
    from adit.neb_setup import neb_in, NebOptions

    text = neb_in(["&CONTROL\n   calculation = 'scf'\n/\nATOMIC_SPECIES\n Si 28.0 Si.upf\n\nATOMIC_POSITIONS crystal\n Si 0.0 0.0 0.0\n"] * 3,
                  NebOptions(images=1), 1)
    assert "nstep_path = 1" in text and "num_of_images = 3" in text


def test_unit_constants_match_hand_calculation():
    from adit.analysis.compare import EV_KCAL_MOL, EV_KJ_MOL
    from adit.analysis.readers import EV_ANG3_GPA
    from adit.analysis.thermo import EV_J_MOL
    from adit.outputs import EV_ANG3_TO_GPA, HARTREE_BOHR3_TO_GPA, HARTREE_BOHR_TO_EV_ANG, RY_BOHR_TO_EV_ANG

    e, a0 = 1.602176634e-19, 0.529177210903e-10
    assert HARTREE_BOHR_TO_EV_ANG == pytest.approx(27.211386245988 / 0.529177210903, rel=1e-6)
    assert RY_BOHR_TO_EV_ANG == pytest.approx(HARTREE_BOHR_TO_EV_ANG / 2, rel=1e-12)
    assert HARTREE_BOHR3_TO_GPA == pytest.approx(27.211386245988 * e / a0 ** 3 / 1e9, rel=1e-6)
    from adit.spec import HARTREE_PER_BOHR_IN_EV_PER_ANG
    assert HARTREE_BOHR_TO_EV_ANG == HARTREE_PER_BOHR_IN_EV_PER_ANG
    assert EV_ANG3_TO_GPA == pytest.approx(e / 1e-30 / 1e9, rel=1e-8) and EV_ANG3_GPA == EV_ANG3_TO_GPA
    assert EV_J_MOL == pytest.approx(96485.33212, rel=1e-9)
    assert EV_KJ_MOL == pytest.approx(EV_J_MOL / 1000, rel=1e-8)
    assert EV_KCAL_MOL == pytest.approx(EV_KJ_MOL / 4.184, rel=1e-8)


def test_conformers_check_the_atom_order_of_fixed_atoms(tmp_path):
    pytest.importorskip("rdkit")
    from adit.conformers import ConformerError, write_conformers

    from tests.conftest import make_fake_skset
    make_fake_skset(tmp_path / "slakos", "fake-1-0", ["H", "C", "N", "O"])
    cfg = cfg_for(tmp_path / "slakos")
    ok = water_spec(method=XtbMethod(), structure=Structure(source="smiles", source_ref="O",
                                                            atoms=water_spec().structure.atoms, fixed_atoms=[0]))
    write_conformers(ok, cfg, tmp_path / "c1", n=2, rmsd=0.3)
    st = water_spec().structure.atoms.model_copy(update={"symbols": ["H", "O", "H"]})
    bad = water_spec(method=XtbMethod(), structure=Structure(source="smiles", source_ref="O", atoms=st, fixed_atoms=[0]))
    with pytest.raises(ConformerError) as ex:
        write_conformers(bad, cfg, tmp_path / "c2", n=2, rmsd=0.3)
    assert "固定原子" in str(ex.value) or "fixed-atom" in str(ex.value)


def test_vasp_neb_stops_above_two_digit_image_dirs(tmp_path):
    from ase import Atoms

    from tests.conftest import make_fake_skset
    from adit.neb_setup import NebError, NebOptions, write_neb
    from adit.spec import VaspMethod

    make_fake_skset(tmp_path / "slakos", "fake-1-0", ["H", "C", "N", "O"])
    cfg = cfg_for(tmp_path / "slakos")
    a = Atoms("OH", positions=[(4, 4, 4), (5, 4, 4)], cell=[8, 8, 8], pbc=True)
    b = Atoms("OH", positions=[(4, 4, 4), (4, 5, 4)], cell=[8, 8, 8], pbc=True)
    from ase.io import write as ase_write
    ase_write(tmp_path / "a.extxyz", a)
    ase_write(tmp_path / "b.extxyz", b)
    spec = water_spec(structure=Structure(source="file", source_ref="a", atoms=AtomsData.from_ase(a), multiplicity=2),
                      method=VaspMethod(ispin=2), kpoints=KPoints(), task=Task(type="geometry_optimization", max_steps=50))
    with pytest.raises(NebError) as ex:
        write_neb(spec, cfg, tmp_path / "neb", "a.extxyz", "b.extxyz", NebOptions(images=99), where=tmp_path)
    assert "00" in str(ex.value)
    assert not (tmp_path / "neb").exists()
