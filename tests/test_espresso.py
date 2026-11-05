
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.codes.upf import UpfLibrary
from adit.project import ProjectError, build_project, write_project
from adit.spec import AtomsData, CalculationSpec, EspressoMethod, KPoints, MDSettings, Runtime, Structure, Task
from adit.structure import from_bulk
from tests.conftest import cfg_for, water_spec, run_log_tail

PWX = os.environ.get("ADIT_PWX_EXE") or shutil.which("pw.x") or ""
REPO = Path(__file__).resolve().parent.parent
REAL_PSEUDO = REPO / "pseudo"


def make_fake_upf(root: Path, pseudo_set="fake-pbe") -> Path:
    d = root / pseudo_set
    d.mkdir(parents=True)
    for el, z, wfc in (("Si", 4.0, 44.0), ("H", 1.0, 46.0), ("O", 6.0, 47.0)):
        (d / f"{el}.pbe-fake.UPF").write_text(
            f'<UPF version="2.0.1">\n<PP_INFO>fake</PP_INFO>\n<PP_HEADER\n   element="{el}"\n   pseudo_type="USPP"\n   functional="PBE"\n'
            f'   z_valence="{z:.12e}"\n   wfc_cutoff="{wfc:.12e}"\n   rho_cutoff="{wfc * 4:.12e}"/>\n</UPF>\n', encoding="utf-8")
    (d / "Si.pbe-second.UPF").write_text('<UPF version="2.0.1">\n<PP_HEADER element="Si" z_valence="4.0" functional="PBE" pseudo_type="PAW"/>\n</UPF>\n', encoding="utf-8")
    (d / "LICENSE").write_text("fake license\n", encoding="utf-8")
    return root


def si_spec(pseudo_set="fake-pbe", **kw) -> CalculationSpec:
    a = from_bulk("Si")
    d = dict(structure=Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(a)),
             method=EspressoMethod(pseudo_set=pseudo_set, pseudo={"Si": "Si.pbe-fake.UPF"}, ecutwfc=40, conv_thr=1e-8),
             kpoints=KPoints(mode="mesh", mesh=(4, 4, 4), shift=(0.5, 0.5, 0.5)),
             task=Task(type="geometry_optimization", max_steps=20), runtime=Runtime(profile="local", mpiprocs=2, omp_threads=1, job_name="si"))
    d.update(kw)
    return CalculationSpec(**d)


@pytest.fixture
def cfg_qe(sk_root, tmp_path):
    cfg = cfg_for(sk_root)
    cfg.pseudo_root = str(make_fake_upf(tmp_path / "pseudo"))
    return cfg


def test_upf_library(tmp_path):
    lib = UpfLibrary(make_fake_upf(tmp_path), "fake-pbe")
    assert lib.files_for("Si") == ["Si.pbe-fake.UPF", "Si.pbe-second.UPF"] and lib.files_for("H") == ["H.pbe-fake.UPF"]
    h = lib.header("Si.pbe-fake.UPF")
    assert (h.element, h.z_valence, h.functional, h.wfc_cutoff) == ("Si", 4.0, "PBE", 44.0)
    assert "LICENSE" in lib.doc_files()


def test_pw_in_and_copies(cfg_qe):
    files = build_project(si_spec(), cfg_qe)
    inp = files.texts["pw.in"]
    for key in ["calculation      = 'relax'", "nstep            = 20", "ecutwfc          = 40", "conv_thr         = 1e-08",
                "K_POINTS automatic", "4 4 4  1 1 1", "Si.pbe-fake.UPF", "ATOMIC_POSITIONS crystal", "forc_conv_thr"]:
        assert key in inp, key
    assert sorted(files.copies) == ["pseudo/LICENSE", "pseudo/Si.pbe-fake.UPF"]
    assert files.texts["submit.sh"].rstrip().endswith("mpirun -np 2 pw.x -in pw.in > output.log 2>&1")


def test_ambiguous_or_missing_pseudo(cfg_qe):
    spec = si_spec(method=EspressoMethod(pseudo_set="fake-pbe", ecutwfc=40))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_qe)
    assert "複数" in ex.value.errors[0].message
    spec = si_spec(method=EspressoMethod(pseudo_set="fake-pbe", pseudo={"Si": "H.pbe-fake.UPF"}, ecutwfc=40))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_qe)
    assert "使えません" in ex.value.errors[0].message
    spec = si_spec(method=EspressoMethod(pseudo_set="nope", ecutwfc=40))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_qe)
    assert ex.value.errors[0].location == "method.pseudo_set"
    with pytest.raises(ProjectError) as ex:
        build_project(si_spec(method=EspressoMethod(pseudo_set="fake-pbe", pseudo={"Si": "Si.pbe-fake.UPF"})), cfg_qe)
    assert ex.value.errors[0].location == "method.ecutwfc"


def test_vc_relax_smearing_fix_and_extra(cfg_qe):
    st = si_spec().structure
    spec = si_spec(structure=Structure(**{**st.model_dump(), "fixed_atoms": [0], "fixed_axes": {"1": (True, True, False)}}),
                   method=EspressoMethod(pseudo_set="fake-pbe", pseudo={"Si": "Si.pbe-fake.UPF"}, ecutwfc=30, ecutrho=240,
                                         occupations="smearing", smearing="marzari-vanderbilt", degauss=0.02,
                                         extra={"system": {"nbnd": 12}, "electrons": {"diagonalization": "cg"}}),
                   task=Task(type="geometry_optimization", relax_cell="volume_only", max_steps=5))
    inp = build_project(spec, cfg_qe).texts["pw.in"]
    for key in ["calculation      = 'vc-relax'", "cell_dofree      = 'volume'", "smearing         = 'marzari-vanderbilt'", "degauss          = 0.02",
                "ecutrho          = 240", "nbnd             = 12", "diagonalization  = 'cg'"]:
        assert key in inp, key
    lines = inp.splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith("ATOMIC_POSITIONS"))
    assert lines[k + 1].split()[4:] == ["0", "0", "0"] and lines[k + 2].split()[4:] == ["1", "1", "0"]  # if_pos


def test_assume_isolated_is_optional(cfg_qe):
    method = EspressoMethod(pseudo_set="fake-pbe", pseudo={"Si": "Si.pbe-fake.UPF"}, ecutwfc=40, assume_isolated="2D")
    inp = build_project(si_spec(method=method), cfg_qe).texts["pw.in"].lower()
    assert "assume_isolated" in inp and "'2d'" in inp


def test_molecule_needs_box_and_parity(cfg_qe):
    spec = water_spec(method=EspressoMethod(pseudo_set="fake-pbe", ecutwfc=30))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_qe)
    assert any("周期セル" in e.message for e in ex.value.errors)


@pytest.mark.parametrize("ensemble", ["NVE", "NVT", "NPT"])
def test_md_trajectory_interval_is_written_to_control(cfg_qe, ensemble):
    spec = si_spec(task=Task(type="molecular_dynamics", md=MDSettings(ensemble=ensemble, dump_interval=3)))
    inp = build_project(spec, cfg_qe).texts["pw.in"]
    control = inp.split("&CONTROL", 1)[1].split("\n/", 1)[0]
    assert "iprint           = 3" in control


@pytest.mark.parametrize("task_type", ["single_point", "geometry_optimization", "vibrations"])
def test_non_md_does_not_receive_md_trajectory_interval(cfg_qe, task_type):
    spec = si_spec(task=Task(type=task_type, md=MDSettings(dump_interval=3)))
    inp = build_project(spec, cfg_qe).texts["pw.in"]
    assert "iprint" not in inp


@pytest.mark.skipif(not (PWX and (REAL_PSEUDO / "pslibrary").is_dir()), reason="pw.x or pseudo/pslibrary is missing")
def test_pwx_real_run(sk_root, tmp_path):
    from adit.results import summarize_run
    cfg = cfg_for(sk_root); cfg.pseudo_root = str(REAL_PSEUDO)
    spec = si_spec(pseudo_set="pslibrary", method=EspressoMethod(pseudo_set="pslibrary", ecutwfc=44, ecutrho=175, conv_thr=1e-8))
    out = tmp_path / "si"
    write_project(spec, cfg, out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900,
                       env={**os.environ, "PATH": f"{Path(PWX).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0, f"{r.stderr}\n{run_log_tail(out)}"
    log = (out / "output.log").read_text(encoding="utf-8")
    assert "JOB DONE" in log
    s = summarize_run(out)
    assert s.converged and s.mermin_energy_hartree is not None and abs(s.mermin_energy_hartree - (-22.8395669929 / 2)) < 1e-4
