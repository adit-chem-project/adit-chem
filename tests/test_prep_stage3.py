
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, molecule
from ase.io import write

from tests.conftest import cfg_for, water_spec
from adit.spec import (AtomsData, CalculationSpec, Cp2kMethod, KPoints, MlipMethod, OrcaMethod, Structure, Task,
                        VaspMethod, XtbMethod)

from adit.spec import HARTREE_PER_BOHR3_IN_GPA as HA_BOHR3_GPA, HARTREE_PER_BOHR_IN_EV_PER_ANG as HA_BOHR


def _box_oh():
    a = Atoms("OH", positions=[(4, 4, 4), (5, 4, 4)], cell=[8, 8, 8], pbc=True)
    b = Atoms("OH", positions=[(4, 4, 4), (4, 5, 4)], cell=[8, 8, 8], pbc=True)
    return a, b


# ---- NEB ----
def test_endpoint_checks():
    from adit.neb_setup import check_endpoints

    a, b = _box_oh()
    assert check_endpoints(a, b) == []
    assert "原子の数" in check_endpoints(a, a + Atoms("H", [(1, 1, 1)]))[0] or "numbers" in check_endpoints(a, a + Atoms("H", [(1, 1, 1)]))[0]
    swapped = Atoms("HO", positions=b.positions, cell=b.cell, pbc=True)
    assert any("1" in e for e in check_endpoints(a, swapped))
    c = b.copy()
    c.cell = [8, 8, 9]
    assert check_endpoints(a, c)


def test_neb_vasp_native(sk_root, tmp_path):
    from adit.neb_setup import NebError, NebOptions, write_neb

    a, b = _box_oh()
    write(tmp_path / "a.extxyz", a)
    write(tmp_path / "b.extxyz", b)
    spec = water_spec(structure=Structure(source="file", source_ref="a", atoms=AtomsData.from_ase(a), multiplicity=2), method=VaspMethod(ispin=2),
                      kpoints=KPoints(), task=Task(type="geometry_optimization", max_steps=50))
    out = tmp_path / "neb"
    dirs = write_neb(spec, cfg_for(sk_root), out, "a.extxyz", "b.extxyz", NebOptions(images=3, vasp_spring=-5.0), where=tmp_path)
    assert [p.name for p in sorted(out.glob("0?"))] == ["00", "01", "02", "03", "04"]
    incar = (out / "INCAR").read_text(encoding="utf-8")
    assert "IMAGES = 3" in incar and "SPRING = -5" in incar and "LCLIMB" not in incar and "IBRION = 2" in incar
    assert not (out / "POSCAR").exists() and (out / "04" / "POSCAR").is_file()
    assert "NSW = 0" in (out / "endpoints" / "initial" / "INCAR").read_text(encoding="utf-8") and len(dirs) == 3
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "最高エネルギーの像は鞍点そのものとは限らない" in readme
    assert "  POSCAR        構造" not in readme and "00 … 04/POSCAR" in readme
    js = json.loads((out / "neb.json").read_text(encoding="utf-8"))
    assert js["n_images_total"] == 5 and len(js["path_length_ang"]) == 5
    from ase.io import read
    mid = read(out / "02" / "POSCAR")
    assert np.allclose(sorted(mid.positions[:, 0]), sorted(read(out / "00" / "POSCAR").positions[:, 0]), atol=1.0)
    with pytest.raises(NebError, match="VTST"):
        write_neb(spec, cfg_for(sk_root), tmp_path / "n2", "a.extxyz", "b.extxyz", NebOptions(images=3, climb=True), where=tmp_path)
    with pytest.raises(NebError):
        write_neb(spec, cfg_for(sk_root), tmp_path / "n3", "spec", "b.extxyz", NebOptions(images=0), where=tmp_path)


def test_neb_qe_native(sk_root, tmp_path):
    from tests.test_espresso import make_fake_upf, si_spec
    from adit.neb_setup import NebOptions, write_neb

    spec = si_spec().model_copy(update={"task": Task(type="single_point", max_steps=30)})
    end = spec.atoms
    end.positions[1] += [0.3, 0.0, 0.0]
    write(tmp_path / "end.extxyz", end)
    cfg = cfg_for(sk_root)
    cfg.pseudo_root = str(make_fake_upf(tmp_path / "pseudo"))
    out = tmp_path / "qe"
    write_neb(spec, cfg, out, "spec", "end.extxyz", NebOptions(images=3, climb=True, qe_opt_scheme="broyden"), where=tmp_path)
    t = (out / "neb.in").read_text(encoding="utf-8")
    assert t.startswith("BEGIN\nBEGIN_PATH_INPUT\n&PATH\n")
    for w in ("num_of_images = 5", "nstep_path = 30", "CI_scheme = 'auto'", "opt_scheme = 'broyden'", "END_PATH_INPUT", "BEGIN_ENGINE_INPUT", "END_ENGINE_INPUT"):
        assert w in t, w
    assert t.count("INTERMEDIATE_IMAGE") == 3 and t.count("FIRST_IMAGE") == 1 and t.count("LAST_IMAGE") == 1
    assert "calculation" not in t and "&CELL" not in t and t.index("END_POSITIONS") < t.index("K_POINTS")
    assert "neb.x -inp neb.in" in (out / "submit.sh").read_text(encoding="utf-8") and not (out / "pw.in").exists()
    assert "No climbing image is used" not in (out / "README.txt").read_text(encoding="utf-8")


def test_neb_cp2k_native(sk_root, tmp_path, monkeypatch):
    from adit.neb_setup import NebError, NebOptions, write_neb

    monkeypatch.delenv("CP2K_DATA_DIR", raising=False)
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    a = Atoms("H2", positions=[(4, 4, 4), (4.74, 4, 4)], cell=[8, 8, 8], pbc=True)
    b = Atoms("H2", positions=[(4, 4, 4), (4, 4.74, 4)], cell=[8, 8, 8], pbc=True)
    write(tmp_path / "b.extxyz", b)
    m = Cp2kMethod(xc="PBE", cutoff_ry=200, rel_cutoff_ry=40, basis={"H": "DZVP-MOLOPT-GTH"}, potential={"H": "GTH-PBE-q1"})
    spec = water_spec(structure=Structure(source="file", source_ref="a", atoms=AtomsData.from_ase(a)), method=m, kpoints=KPoints(), task=Task(max_steps=40))
    out = tmp_path / "cp"
    write_neb(spec, cfg_for(sk_root), out, "spec", "b.extxyz", NebOptions(images=3, climb=True, cp2k_k_spring=0.05), where=tmp_path)
    t = (out / "cp2k.inp").read_text(encoding="utf-8")
    assert "RUN_TYPE BAND" in t and "NUMBER_OF_REPLICA 5" in t and "BAND_TYPE CI-NEB" in t and "K_SPRING 0.05" in t and "MAX_STEPS 40" in t
    assert t.count("COORD_FILE_NAME replica_") == 5 and (out / "replica_04.xyz").is_file()
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "像ごとのエネルギーと像間の距離を読みます" in readme and "まだ読みません" not in readme
    mol = spec.model_copy(update={"structure": Structure(source="preset", source_ref="H2", atoms=AtomsData.from_ase(molecule("H2"))), "kpoints": None})
    write(tmp_path / "m.xyz", molecule("H2"))
    with pytest.raises(NebError):
        write_neb(mol.model_copy(update={"method": m.model_copy(update={"poisson_solver": "MT", "isolated_box_ang": 8})}), cfg_for(sk_root),
                  tmp_path / "x", "spec", "m.xyz", NebOptions(images=2), where=tmp_path)


def test_neb_images_mode_any_code(sk_root, tmp_path):
    from adit.neb_setup import NebOptions, write_neb

    spec = water_spec(task=Task(type="geometry_optimization"))
    end = spec.atoms
    end.positions[1] += [0.0, 0.3, 0.0]
    write(tmp_path / "e.xyz", end)
    dirs = write_neb(spec, cfg_for(sk_root), tmp_path / "img", "spec", "e.xyz", NebOptions(images=2, mode="images"), where=tmp_path)
    assert [d.name for d in dirs] == ["image_00", "image_01", "image_02", "image_03"]
    assert "Driver {}" in (dirs[1] / "dftb_in.hsd").read_text(encoding="utf-8")
    assert json.loads((tmp_path / "img" / "scan.json").read_text(encoding="utf-8"))["values"] == ["00", "01", "02", "03"]


def _diamond_spec(**kw):
    kw.setdefault("task", Task(type="single_point"))
    return water_spec(structure=Structure(source="bulk", source_ref="C", atoms=AtomsData.from_ase(bulk("C", "diamond", a=3.57))),
                      kpoints=KPoints(mode="mesh", mesh=(2, 2, 2)), **kw)


def _fake_dftb_forces(d: Path):
    from ase.calculators.emt import EMT

    a = CalculationSpec.load(d / "spec.json").atoms
    a.calc = EMT()
    f = a.get_forces() / HA_BOHR
    rows = "\n".join(f"{i + 1:6d} {x:22.14f} {y:22.14f} {z:22.14f}" for i, (x, y, z) in enumerate(f))
    (d / "detailed.out").write_text(f"Total Forces\n{rows}\n\nMaximal derivative component: 0.0 au\n", encoding="utf-8")


def test_phonons_ase_backend_collect(sk_root, tmp_path):
    from adit.analysis.phonons import read_band_yaml, read_total_dos
    from adit.phonon_setup import PhononError, collect, write_phonons

    out = tmp_path / "ph"
    with pytest.raises(PhononError):
        write_phonons(_diamond_spec(), cfg_for(sk_root), tmp_path / "nw", (2, 2, 2), backend="ase", dos_mesh=(8, 8, 8))
    with pytest.raises(PhononError):
        write_phonons(_diamond_spec(), cfg_for(sk_root), tmp_path / "big", (2, 2, 2), backend="ase", dos_mesh=(60, 60, 60), dos_width_thz=0.5)
    dirs = write_phonons(_diamond_spec(), cfg_for(sk_root), out, (2, 2, 2), backend="ase", dos_mesh=(8, 8, 8), dos_width_thz=0.5)
    rec = json.loads((out / "phonons.json").read_text(encoding="utf-8"))
    assert len(dirs) == 12 and rec["ase_keys"][:2] == ["0x-", "0x+"] and rec["distance_ang"] == 0.01
    assert "KPointsAndWeights" in (dirs[0] / "dftb_in.hsd").read_text(encoding="utf-8") and (out / "phonon_collect.py").is_file()
    with pytest.raises(PhononError, match="12"):
        collect(out)
    for d in dirs:
        _fake_dftb_forces(d)
    res = collect(out)
    b = read_band_yaml(out / "band.yaml")
    assert b["frequency"].shape[1] == 6 and b["labels"][0][0] == "G"
    assert np.sort(np.abs(b["frequency"][0]))[:3].max() < 0.05
    e, g, _ = read_total_dos(out / "total_dos.dat")
    assert len(e) > 10 and g.max() > 0 and "THz" in res["summary"]


def test_phonons_phonopy_backend(sk_root, tmp_path):
    pytest.importorskip("phonopy")
    from adit.analysis.phonons import read_band_yaml
    from adit.phonon_setup import collect, write_phonons

    out = tmp_path / "pp"
    dirs = write_phonons(_diamond_spec(), cfg_for(sk_root), out, (2, 2, 2), backend="phonopy", dos_mesh=(10, 10, 10))
    assert (out / "phonopy_disp.yaml").is_file() and len(dirs) == json.loads((out / "phonons.json").read_text(encoding="utf-8"))["n_displacements"]
    for d in dirs:
        _fake_dftb_forces(d)
    collect(out)
    assert read_band_yaml(out / "band.yaml")["frequency"].shape[1] == 6 and (out / "total_dos.dat").is_file()


def test_ase_backend_matches_ase_phonons_run(sk_root, tmp_path):
    from ase.calculators.emt import EMT
    from ase.phonons import Phonons

    from adit.analysis.phonons import read_band_yaml
    from adit.phonon_setup import BAND_POINTS, _thz_per_ev, band_paths, collect, write_phonons

    spec = _diamond_spec()
    out = tmp_path / "mine"
    for d in write_phonons(spec, cfg_for(sk_root), out, (2, 2, 2), backend="ase"):
        _fake_dftb_forces(d)
    collect(out)
    mine = read_band_yaml(out / "band.yaml")["frequency"]
    ph = Phonons(spec.atoms, EMT(), supercell=(2, 2, 2), delta=0.01, name=str(tmp_path / "direct"))
    ph.run()
    ph.read()
    paths, _, _, _ = band_paths(spec.atoms.cell, BAND_POINTS)
    ref = np.asarray(ph.band_structure(np.vstack(paths), verbose=False)) * _thz_per_ev()
    assert mine.shape == ref.shape and np.abs(mine - ref).max() < 1e-5


def test_phonon_checks(sk_root, tmp_path):
    from adit.phonon_setup import PhononError, parse_dim, write_phonons

    assert parse_dim("2x2x1") == (2, 2, 1)
    for bad in ("2x2", "0x1x1", "a"):
        with pytest.raises(PhononError):
            parse_dim(bad)
    with pytest.raises(PhononError):
        write_phonons(water_spec(), cfg_for(sk_root), tmp_path / "w", (2, 2, 2), backend="ase")
    with pytest.raises(PhononError):
        write_phonons(_diamond_spec().model_copy(update={"method": XtbMethod()}), cfg_for(sk_root), tmp_path / "x", (2, 2, 2), backend="ase")


C_TRUE = np.array([[500, 100, 100, 0, 0, 0], [100, 500, 100, 0, 0, 0], [100, 100, 500, 0, 0, 0],
                   [0, 0, 0, 300, 0, 0], [0, 0, 0, 0, 300, 0], [0, 0, 0, 0, 0, 300]], dtype=float)


def _fake_stress(d: Path, strain_voigt: np.ndarray, flip: bool):
    from adit.elastic_setup import VOIGT

    s_v = C_TRUE @ strain_voigt + np.array([0.3, 0.3, 0.3, 0, 0, 0])
    s = np.zeros((3, 3))
    for k, (i, j) in VOIGT.items():
        s[i, j] = s[j, i] = s_v[k - 1]
    s_au = s / HA_BOHR3_GPA * (-1 if flip else 1)
    p = -np.trace(s) / 3 / HA_BOHR3_GPA
    rows = "\n".join(" ".join(f"{x:18.12f}" for x in r) for r in s_au)
    (d / "detailed.out").write_text(f"Total stress tensor\n{rows}\n\nPressure:   {p:.10E} au  {p * 2.9421e13:.6E} Pa\n", encoding="utf-8")


@pytest.mark.parametrize("flip", [False, True])
def test_elastic_generate_and_fit(sk_root, tmp_path, flip):
    from adit.elastic_setup import collect, write_elastic

    spec = _diamond_spec()
    out = tmp_path / "el"
    dirs = write_elastic(spec, cfg_for(sk_root), out, [-0.01, 0.01], [1, 4])
    assert [d.name for d in dirs] == ["e0", "e1_-0.01", "e1_+0.01", "e4_-0.01", "e4_+0.01"]
    c0 = spec.atoms.cell.array
    c1 = CalculationSpec.load(out / "e1_+0.01" / "spec.json").atoms.cell.array
    assert np.allclose(c1, c0 @ (np.eye(3) + np.diag([0.01, 0, 0])))
    rec = json.loads((out / "elastic.json").read_text(encoding="utf-8"))
    for r in rec["runs"]:
        v = np.zeros(6)
        if r["component"]:
            v[r["component"] - 1] = r["strain"]
        _fake_stress(out / r["dir"], v, flip)
    res = collect(out)
    c = res["c_gpa"]
    assert abs(c[0][0] - 500) < 1e-6 and abs(c[1][0] - 100) < 1e-6 and abs(c[3][3] - 300) < 1e-6 and c[0][1] is None
    assert (out / "elastic_constants.csv").is_file()


def test_elastic_checks_and_qe_stress(sk_root, tmp_path):
    from adit.elastic_setup import ElasticError, write_elastic
    from adit.outputs import read_stress

    for strains, comps in (([0.0, 0.01], [1]), ([0.01, 0.01], [1]), ([0.01], [7]), ([1.5], [1])):
        with pytest.raises(ElasticError):
            write_elastic(_diamond_spec(), cfg_for(sk_root), tmp_path / "e", strains, comps)
    with pytest.raises(ElasticError):
        write_elastic(_diamond_spec(task=Task(type="geometry_optimization", relax_cell="shape_and_volume")), cfg_for(sk_root), tmp_path / "e", [0.01], [1])
    from tests.test_espresso import si_spec
    d = tmp_path / "qe"
    d.mkdir()
    si_spec().save(d / "spec.json")
    (d / "output.log").write_text("          total   stress  (Ry/bohr**3)                   (kbar)     P=       20.11\n"
                                  "   0.00013671  -0.00000000  -0.00000000           20.11       -0.00       -0.00\n"
                                  "  -0.00000000   0.00013671  -0.00000000           -0.00       20.11       -0.00\n"
                                  "  -0.00000000  -0.00000000   0.00013671           -0.00       -0.00       20.11\n", encoding="utf-8")
    assert np.allclose(np.diag(read_stress(d)), -2.011)


def test_orca_ts_and_irc_inputs(sk_root, tmp_path):
    from adit.ts_setup import TsError, ts_spec, write_ts
    from adit.validate import validate

    spec = water_spec(method=OrcaMethod(method="HF", basis="def2-SVP"))
    ts = ts_spec(spec, "ts").model_copy(update={"method": ts_spec(spec, "ts").method.model_copy(update={"ts_recalc_hess": 5})})
    write_ts(ts, cfg_for(sk_root), tmp_path / "ts", "ts")
    t = (tmp_path / "ts" / "orca.inp").read_text(encoding="utf-8")
    assert "! HF def2-SVP OptTS Freq" in t and "Calc_Hess true" in t and "Recalc_Hess 5" in t
    irc = ts_spec(spec, "irc")
    irc = irc.model_copy(update={"method": irc.method.model_copy(update={"irc_max_iter": 20, "irc_direction": "forward"})})
    write_ts(irc, cfg_for(sk_root), tmp_path / "irc", "irc")
    t = (tmp_path / "irc" / "orca.inp").read_text(encoding="utf-8")
    assert "! HF def2-SVP Freq IRC" in t and "%irc\n   MaxIter 20\n   Direction forward\nend" in t
    bad = spec.model_copy(update={"method": spec.method.model_copy(update={"ts_search": True}), "task": Task(type="single_point")})
    assert any(e.location == "method.ts_search" for e in validate(bad, cfg_for(sk_root)))
    both = spec.model_copy(update={"method": spec.method.model_copy(update={"ts_search": True, "irc": True})})
    assert any(e.location == "method.irc" for e in validate(both, cfg_for(sk_root)))
    dirs = write_ts(spec, cfg_for(sk_root), tmp_path / "both", "ts+irc")
    assert [d.name for d in dirs] == ["stage_01_ts", "stage_02_irc"]
    assert "OptTS" in (dirs[0] / "orca.inp").read_text(encoding="utf-8") and "Freq IRC" in (dirs[1] / "orca.inp").read_text(encoding="utf-8")
    assert "handoff.py orca ../stage_01_ts geometry_optimization" in (dirs[1] / "submit.sh").read_text(encoding="utf-8")
    with pytest.raises(TsError):
        write_ts(water_spec(method=XtbMethod()), cfg_for(sk_root), tmp_path / "x", "ts")


def _stubs(root: Path) -> Path:
    (root / "mace" / "calculators").mkdir(parents=True)
    (root / "mace" / "__init__.py").write_text("", encoding="utf-8")
    (root / "mace" / "calculators" / "__init__.py").write_text(
        "from ase.calculators.emt import EMT\nCALLS = []\n\ndef mace_mp(**kw):\n    return EMT()\n\ndef mace_off(**kw):\n    return EMT()\n", encoding="utf-8")
    (root / "tblite").mkdir()
    (root / "tblite" / "__init__.py").write_text("", encoding="utf-8")
    (root / "tblite" / "ase.py").write_text("from ase.calculators.emt import EMT\n\ndef TBLite(**kw):\n    assert kw['method'] == 'GFN2-xTB'\n    return EMT()\n", encoding="utf-8")
    (root / "sella").mkdir()
    (root / "sella" / "__init__.py").write_text(
        "from ase.optimize import BFGS\n\nclass Sella(BFGS):\n    pass\n\nclass IRC(BFGS):\n    def run(self, fmax=0.05, steps=10, direction='forward'):\n        return True\n", encoding="utf-8")
    return root


def _run(d: Path, script: str, stubs: Path) -> dict:
    import os
    env = dict(os.environ, PYTHONPATH=str(stubs), OMP_NUM_THREADS="1")
    r = subprocess.run([sys.executable, script], cwd=d, env=env, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    return json.loads((d / "results.json").read_text(encoding="utf-8"))


def _cu_spec(task: Task, **mkw):
    return water_spec(structure=Structure(source="bulk", source_ref="Cu", atoms=AtomsData.from_ase(bulk("Cu", "fcc", a=3.6, cubic=True))),
                      method=MlipMethod(model_family="mace_mp", model="small", **mkw), task=task)


def test_mlip_generator_and_script_runs(sk_root, tmp_path):
    from adit.project import write_project

    stubs = _stubs(tmp_path / "stubs")
    cases = {
        "sp": _cu_spec(Task(type="single_point")),
        "opt": _cu_spec(Task(type="geometry_optimization", optimizer="LBFGS", max_steps=5, relax_cell="volume_only")),
        "npt": _cu_spec(Task(type="molecular_dynamics", md={"ensemble": "NPT", "thermostat": "nose_hoover", "steps": 3, "dump_interval": 1})),
        "nvt": water_spec(method=MlipMethod(model_family="mace_off"), task=Task(type="molecular_dynamics", md={"ensemble": "NVT", "thermostat": "langevin", "steps": 3, "dump_interval": 1})),
        "vib": water_spec(method=MlipMethod(model_family="mace_mp", dtype="float64"), task=Task(type="vibrations")),
    }
    for name, s in cases.items():
        d = tmp_path / name
        write_project(s, cfg_for(sk_root), d)
        assert "python3 run_mlip.py > output.log" in (d / "submit.sh").read_text(encoding="utf-8")
        res = _run(d, "run_mlip.py", stubs)
        assert "energy_ev" in res and res["versions"]["ase"] != "not installed", name
        if name == "sp":
            assert len(res["stress_voigt_ev_per_ang3"]) == 6
        if name == "vib":
            assert len(res["frequencies_cm1"]) == 9
        if name in ("opt", "npt", "nvt"):
            assert (d / "trajectory.extxyz").is_file()
    st = json.loads((tmp_path / "opt" / "mlip_settings.json").read_text(encoding="utf-8"))
    assert st["model"] == "small" and st["relax_cell"] == "volume_only" and st["optimizer"] == "LBFGS"
    readme = (tmp_path / "sp" / "README.txt").read_text(encoding="utf-8")
    assert "pip install ase mace-torch" in readme and "MIT" in readme and "pair_style = mace no_domain_decomposition" in readme


def test_mlip_checks_and_model_file(sk_root, tmp_path):
    from adit.project import write_project
    from adit.validate import validate

    cfg = cfg_for(sk_root)
    locs = lambda s: {e.location for e in validate(s, cfg)}  # noqa: E731
    assert "method.model_family" in locs(water_spec(method=MlipMethod()))
    assert "structure.charge" in locs(water_spec(method=MlipMethod(model_family="chgnet"), structure=water_spec().structure.model_copy(update={"charge": 1})))
    assert "method.model" in locs(water_spec(method=MlipMethod(model_family="chgnet", model="x")))
    assert "method.dispersion" in locs(water_spec(method=MlipMethod(model_family="mace_off", dispersion=True)))
    assert "task.optimizer" in locs(water_spec(method=MlipMethod(model_family="mace_mp")))
    assert "task.type" in locs(_cu_spec(Task(type="band_structure")))
    assert "task.md.thermostat" in locs(_cu_spec(Task(type="molecular_dynamics", md={"ensemble": "NPT", "thermostat": "berendsen"})))
    model = tmp_path / "my.model"
    model.write_bytes(b"x")
    s = _cu_spec(Task(type="single_point")).model_copy(update={"method": MlipMethod(model_family="mace_mp", model=str(model))})
    write_project(s, cfg, tmp_path / "m")
    assert (tmp_path / "m" / "my.model").is_file() and json.loads((tmp_path / "m" / "mlip_settings.json").read_text(encoding="utf-8"))["model"] == "my.model"
    again = CalculationSpec.from_json(s.to_json())
    assert again == s and again.method.code == "mlip"


def test_sella_script(sk_root, tmp_path):
    from adit.ts_setup import TsError, write_sella

    stubs = _stubs(tmp_path / "stubs")
    spec = water_spec(method=XtbMethod(), task=Task(type="geometry_optimization", max_steps=5))
    out = tmp_path / "sella"
    write_sella(spec, cfg_for(sk_root), out)
    st = json.loads((out / "sella_settings.json").read_text(encoding="utf-8"))
    assert st["calculator"] == {"kind": "tblite", "method": "GFN2-xTB", "charge": 0, "multiplicity": 1, "accuracy": 1.0} and st["irc"] is True
    res = _run(out, "run_sella.py", stubs)
    assert "ts" in res and "irc_forward" in res and "irc_reverse" in res and (out / "ts.extxyz").is_file()
    with pytest.raises(TsError):
        write_sella(water_spec(method=OrcaMethod()), cfg_for(sk_root), tmp_path / "o")
    with pytest.raises(TsError):
        write_sella(water_spec(method=XtbMethod(solvation="alpb", solvent="water")), cfg_for(sk_root), tmp_path / "p")


def test_old_spec_json_still_loads():
    repo = Path(__file__).resolve().parent.parent
    for p in sorted(repo.glob("examples/*/spec.json"))[:40]:
        CalculationSpec.load(p)
    o = OrcaMethod.model_validate({"code": "orca", "method": "HF"})
    assert o.ts_search is False and o.irc is False and o.ts_calc_hess is True and o.irc_direction == "both"
