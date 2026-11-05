
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from ase.io import read

from adit.codes.gromacs import gmx_top_dir, scan_includes
from adit.config import Profile, save_config
from adit.project import ProjectError, build_project, write_project
from adit.spec import AtomsData, CalculationSpec, GromacsMethod, MDSettings, Runtime, Structure, Task
from tests.conftest import cfg_for

GRO = """water
    3
    1SOL     OW    1   0.126   0.639   0.322
    1SOL    HW1    2   0.187   0.713   0.290
    1SOL    HW2    3   0.108   0.577   0.248
   3.00000   3.00000   3.00000
"""


@pytest.fixture
def src(tmp_path) -> Path:
    d = tmp_path / "src"
    (d / "myff.ff").mkdir(parents=True)
    (d / "myff.ff" / "forcefield.itp").write_text('[ defaults ]\n1 2 yes 0.5 0.5\n#include "ffnonbonded.itp"\n', encoding="utf-8")
    (d / "myff.ff" / "ffnonbonded.itp").write_text("[ atomtypes ]\n", encoding="utf-8")
    (d / "mol.itp").write_text("[ moleculetype ]\nSOL 2\n", encoding="utf-8")
    (d / "posre.itp").write_text("[ position_restraints ]\n", encoding="utf-8")
    (d / "topol.top").write_text('#include "myff.ff/forcefield.itp"\n#include "mol.itp"\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'
                                 "[ system ]\nx\n[ molecules ]\nSOL 1\n", encoding="utf-8")
    (d / "conf.gro").write_text(GRO, encoding="utf-8")
    return d


def spec_for(src, task=None, **m) -> CalculationSpec:
    atoms = read(src / "conf.gro")
    m = {"topology_file": str(src / "topol.top"), "structure_file": str(src / "conf.gro"), **m}
    return CalculationSpec(structure=Structure(source="file", source_ref=str(src / "conf.gro"), atoms=AtomsData.from_ase(atoms)),
                           method=GromacsMethod(**m),
                           task=task or Task(type="geometry_optimization", optimizer="SteepestDescent", max_steps=100, force_tolerance_ev_per_ang=1.0),
                           runtime=Runtime(profile="local", mpiprocs=1, omp_threads=2, job_name="w"))


def mdp(text: str) -> dict[str, str]:
    return {k.strip(): v.strip() for k, v in (l.split("=", 1) for l in text.splitlines() if "=" in l and not l.startswith(";"))}


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def md(ensemble="NVT", thermostat="csvr", **kw) -> Task:
    return Task(type="molecular_dynamics", md=MDSettings(ensemble=ensemble, thermostat=thermostat, temperature_k=300, timestep_fs=2.0,
                                                         steps=500, dump_interval=50, coupling_time_fs=100, **kw))


def test_scan_includes(src, tmp_path):
    copies, problems, lib = scan_includes(src / "topol.top", "", None)
    assert sorted(copies) == ["mol.itp", "myff.ff/ffnonbonded.itp", "myff.ff/forcefield.itp", "posre.itp"] and problems == [] and lib == []
    lib_dir = tmp_path / "gmxtop"
    (lib_dir / "oplsaa.ff").mkdir(parents=True)
    (lib_dir / "oplsaa.ff" / "spce.itp").write_text("x\n", encoding="utf-8")
    (src / "posre.itp").unlink()
    _, problems, _ = scan_includes(src / "topol.top", "-DPOSRES", lib_dir)
    assert len(problems) == 1 and "posre.itp" in problems[0]
    _, problems, lib = scan_includes(src / "topol.top", "-DPOSRES", None)
    assert problems == [] and lib == ["posre.itp (未確認)"]
    top2 = src / "t2.top"
    top2.write_text('#include "oplsaa.ff/spce.itp"\n#include "../outside.itp"\n', encoding="utf-8")
    (tmp_path / "outside.itp").write_text("x\n", encoding="utf-8")
    copies, problems, lib = scan_includes(top2, "", lib_dir)
    assert lib == ["oplsaa.ff/spce.itp"] and len(problems) == 1 and "outside.itp" in problems[0]


def test_em(sk_root, src):
    files = build_project(spec_for(src), cfg_for(sk_root))
    d = mdp(files.texts["grompp.mdp"])
    assert d["integrator"] == "steep" and d["nsteps"] == "100" and d["emtol"] == "964.8533212"
    assert (d["coulombtype"], d["rcoulomb"], d["rvdw"], d["constraints"], d["cutoff-scheme"]) == ("Cut-off", "1", "1", "none", "Verlet")
    assert "rlist" not in d and "tcoupl" not in d and "gen-vel" not in d
    assert sorted(files.copies) == ["conf.gro", "mol.itp", "myff.ff/ffnonbonded.itp", "myff.ff/forcefield.itp", "posre.itp", "topol.top"]
    run = files.texts["submit.sh"].rstrip().splitlines()[-1]
    assert run == ("gmx grompp -f grompp.mdp -c conf.gro -p topol.top -o adit.tpr > grompp.log 2>&1 && "
                   "gmx mdrun -deffnm adit -ntmpi 1 -ntomp 2 > output.log 2>&1")
    assert "GROMACS" in files.texts["README.txt"] and "adit.edr" in files.texts["README.txt"]


def test_md_stages(sk_root, src, tmp_path):
    d = mdp(build_project(spec_for(src, md()), cfg_for(sk_root)).texts["grompp.mdp"])
    assert (d["integrator"], d["dt"], d["nsteps"], d["nstxout-compressed"], d["tcoupl"], d["tc-grps"], d["tau-t"], d["ref-t"]) == \
        ("md", "0.002", "500", "50", "V-rescale", "System", "0.1", "300")
    assert (d["pcoupl"], d["gen-vel"], d["gen-temp"], d["gen-seed"]) == ("no", "yes", "300", "12345")
    readme = build_project(spec_for(src, md()), cfg_for(sk_root)).texts["README.txt"]
    assert "dt = 2 fs、constraints = none" in readme and "gen-seed = 12345" in readme
    d = mdp(build_project(spec_for(src, md(thermostat="langevin")), cfg_for(sk_root)).texts["grompp.mdp"])
    assert d["integrator"] == "sd" and "tcoupl" not in d and d["tau-t"] == "0.1"
    cpt = tmp_path / "state.cpt"
    cpt.write_bytes(b"x")
    files = build_project(spec_for(src, md("NPT", pressure_bar=1.0, barostat_time_fs=2000), pcoupl="Parrinello-Rahman",
                                   compressibility_per_bar=4.5e-5, checkpoint_file=str(cpt), define="-DPOSRES"), cfg_for(sk_root))
    d = mdp(files.texts["grompp.mdp"])
    assert (d["pcoupl"], d["pcoupltype"], d["tau-p"], d["ref-p"], d["compressibility"]) == ("Parrinello-Rahman", "isotropic", "2", "1", "4.5e-05")
    assert (d["continuation"], d["gen-vel"], d["define"]) == ("yes", "no", "-DPOSRES") and "gen-temp" not in d
    assert files.copies["prev.cpt"] == cpt
    assert "-r conf.gro -t prev.cpt -o" not in files.texts["submit.sh"] and "-o adit.tpr -r conf.gro -t prev.cpt" in files.texts["submit.sh"]


def test_extra_mdp_and_mpi_command(sk_root, src):
    d = mdp(build_project(spec_for(src, md(), extra_mdp={"constraints": "h-bonds", "nstcalcenergy": 50, "comm-mode": "Linear"}),
                          cfg_for(sk_root)).texts["grompp.mdp"])
    assert (d["constraints"], d["nstcalcenergy"], d["comm-mode"]) == ("h-bonds", "50", "Linear")
    cfg = cfg_for(sk_root)
    cfg.profiles["local"] = Profile(kind="direct", commands={"gromacs": "mpirun -np {mpiprocs} gmx_mpi"})
    s = spec_for(src)
    s = s.model_copy(update={"runtime": s.runtime.model_copy(update={"mpiprocs": 4})})
    run = build_project(s, cfg).texts["submit.sh"].rstrip().splitlines()[-1]
    assert run.startswith("gmx_mpi grompp ") and "mpirun -np 4 gmx_mpi mdrun -deffnm adit -ntomp 2 >" in run and "-ntmpi" not in run


def test_things_grompp_would_reject(sk_root, src):
    cfg = cfg_for(sk_root)
    assert any(e.location == "task.md.thermostat" for e in errors(spec_for(src, md(thermostat="berendsen")), cfg))
    assert any(e.location == "task.md.thermostat" for e in errors(spec_for(src, md(thermostat="andersen")), cfg))
    locs = {e.location for e in errors(spec_for(src, md("NPT")), cfg)}
    assert {"method.pcoupl", "method.compressibility_per_bar"} <= locs
    pr = errors(spec_for(src, md("NPT"), pcoupl="Parrinello-Rahman", compressibility_per_bar=4.5e-5), cfg)
    assert [e.location for e in pr] == ["method.pcoupl"] and ".cpt" in pr[0].message
    build_project(spec_for(src, md("NPT"), pcoupl="C-rescale", compressibility_per_bar=4.5e-5), cfg)
    assert any(e.location == "method.rcoulomb_nm" for e in errors(spec_for(src, rcoulomb_nm=1.45, rvdw_nm=1.45), cfg))
    build_project(spec_for(src, rcoulomb_nm=1.4, rvdw_nm=1.4), cfg)
    build_project(spec_for(src, md(), rcoulomb_nm=1.45, rvdw_nm=1.45), cfg)


def test_inputs_and_unsupported(sk_root, src):
    cfg = cfg_for(sk_root)
    assert {e.location for e in errors(spec_for(src, topology_file="", structure_file=""), cfg)} >= {"method.topology_file", "method.structure_file"}
    s = spec_for(src)
    four = s.structure.model_copy(update={"atoms": AtomsData(symbols=["O", "H", "H", "H"], positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
                                                             cell=s.structure.atoms.cell, pbc=(True, True, True))})
    assert any(e.location == "structure.atoms" for e in errors(s.model_copy(update={"structure": four}), cfg))
    assert any(e.location == "structure.charge" for e in errors(s.model_copy(update={"structure": s.structure.model_copy(update={"charge": -1})}), cfg))
    assert any(e.location == "structure.fixed_atoms" for e in errors(s.model_copy(update={"structure": s.structure.model_copy(update={"fixed_atoms": [0]})}), cfg))
    assert any(e.location == "task.type" for e in errors(spec_for(src, Task(type="vibrations")), cfg))
    assert any(e.location == "task.relax_cell" for e in errors(spec_for(src, Task(type="geometry_optimization", relax_cell="shape_and_volume")), cfg))
    assert any(e.location == "method.extra_mdp" for e in errors(spec_for(src, extra_mdp={"bad key": 1}), cfg))
    d = mdp(build_project(spec_for(src, Task(type="single_point")), cfg).texts["grompp.mdp"])
    assert (d["integrator"], d["nsteps"]) == ("md", "0") and "gen-vel" not in d


def test_spec_roundtrip_and_cli(sk_root, src, tmp_path, capsys):
    from adit.cli import main
    spec = spec_for(src, md(), coulombtype="PME", rcoulomb_nm=1.2, rvdw_nm=1.2, extra_mdp={"nstcalcenergy": 50})
    assert CalculationSpec.from_json(spec.to_json()) == spec
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec.save(tmp_path / "spec.json")
    assert main([str(tmp_path / "spec.json"), str(tmp_path / "calc"), "--config", str(cfg_path)]) == 0
    assert (tmp_path / "calc" / "myff.ff" / "ffnonbonded.itp").is_file() and "coulombtype        = PME" in (tmp_path / "calc" / "grompp.mdp").read_text(encoding="utf-8")


GMX = shutil.which("gmx") or ""


@pytest.mark.skipif(not (GMX and gmx_top_dir() and (gmx_top_dir() / "spc216.gro").is_file()), reason="gmx or the spc216.gro shipped with GROMACS is missing")
def test_gromacs_real_run(sk_root, tmp_path):
    top_dir = gmx_top_dir()
    s = tmp_path / "src"
    s.mkdir()
    shutil.copy(top_dir / "spc216.gro", s / "conf.gro")
    (s / "topol.top").write_text('#include "oplsaa.ff/forcefield.itp"\n#include "oplsaa.ff/spce.itp"\n[ system ]\nSPC/E\n[ molecules ]\nSOL 216\n', encoding="utf-8")
    spec = spec_for(s, Task(type="geometry_optimization", optimizer="SteepestDescent", max_steps=20, force_tolerance_ev_per_ang=1.0),
                    coulombtype="PME", rcoulomb_nm=0.85, rvdw_nm=0.85)
    spec = spec.model_copy(update={"runtime": spec.runtime.model_copy(update={"omp_threads": 1})})
    out = tmp_path / "em"
    write_project(spec, cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=300,
                       env={**os.environ, "PATH": f"{Path(GMX).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0, (out / "grompp.log").read_text(encoding="utf-8")[-2000:]
    assert "Steepest Descents" in (out / "output.log").read_text(encoding="utf-8") and (out / "adit.gro").is_file()
