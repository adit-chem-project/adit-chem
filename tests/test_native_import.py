import json

import pytest

from adit.native_import import import_native, NativeImportError, MAX_BYTES


POSCAR = "Si\n1\n3 0 0\n0 3 0\n0 0 3\nSi\n1\nDirect\n0 0 0\n"
QE = """&CONTROL
 calculation = 'scf'
/
&SYSTEM
 ibrav = 0, nat = 1, ntyp = 1, ecutwfc = 30
/
&ELECTRONS
 conv_thr = 1d-8
/
ATOMIC_SPECIES
Si 28.085 Si.UPF
CELL_PARAMETERS angstrom
3 0 0
0 3 0
0 0 3
ATOMIC_POSITIONS crystal
Si 0 0 0
K_POINTS automatic
2 2 2 1 0 1
"""
DATA = """Small explicit type example

1 atoms
1 atom types

0 3 xlo xhi
0 3 ylo yhi
0 3 zlo zhi

Masses

1 28.085 # Si

Atoms # atomic

1 1 0 0 0
"""


def vasp(tmp_path, incar="ENCUT=400; EDIFF=1D-6\nNSW=0\nIBRION=-1\n"):
    (tmp_path / "INCAR").write_text(incar, encoding="utf-8")
    (tmp_path / "POSCAR").write_text(POSCAR, encoding="utf-8")
    (tmp_path / "KPOINTS").write_text("mesh\n0\nGamma\n2 2 2\n0 0 0\n", encoding="utf-8")
    return import_native(tmp_path)


def test_vasp_import_provenance_and_defaults(tmp_path):
    result = vasp(tmp_path)
    assert result.spec is not None, result.report()
    assert result.spec.method.encut == 400
    assert result.spec.method.ediff == 1e-6
    assert result.spec.kpoints.mesh == (2, 2, 2)
    assert result.provenance["method.encut"]["line"] == 1
    assert len(result.files[str(tmp_path / "INCAR")]["sha256"]) == 64
    assert "method.potcar_set" in result.inferred_defaults
    assert result.report()["ready_to_regenerate"] is False
    assert result.unresolved_dependencies[0]["kind"] == "POTCAR library"
    json.dumps(result.report())


@pytest.mark.parametrize("extra", ["NELECT=3", "NSW=10", "ENCUT=300\nENCUT=400", "MAGMOM=1000000000*1"])
def test_vasp_unsupported_never_returns_spec(tmp_path, extra):
    result = vasp(tmp_path, extra)
    assert result.spec is None
    assert result.unknown or result.unsupported


def test_poscar_huge_count_rejected_before_ase(tmp_path, monkeypatch):
    vasp(tmp_path)
    (tmp_path / "POSCAR").write_text(POSCAR.replace("\n1\nDirect", "\n999999999\nDirect"), encoding="utf-8")
    def fail(*args, **kwargs):
        pytest.fail("ASE called before count preflight")
    monkeypatch.setattr("ase.io.read", fail)
    assert import_native(tmp_path).spec is None


def test_qe_import(tmp_path):
    path = tmp_path / "pw.in"
    path.write_text(QE, encoding="utf-8")
    result = import_native(path)
    assert result.spec is not None, result.report()
    assert result.spec.method.conv_thr == 1e-8
    assert result.spec.method.pseudo == {"Si": "Si.UPF"}
    assert result.spec.kpoints.shift == (0.5, 0, 0.5)
    assert result.provenance["method.ecutwfc"]["line"] == 5
    assert result.unresolved_dependencies[0]["file"] == "Si.UPF"
    assert result.unresolved_dependencies[0]["line"] == 11
    assert result.report()["ready_to_regenerate"] is False


@pytest.mark.parametrize("text", [QE.replace("'scf'", "'md'"), QE.replace("ecutwfc = 30", "ecutwfc = 30, noncolin = .true."), QE + "HUBBARD atomic\nU Si-3p 1\n"])
def test_qe_unsupported_is_explicit(tmp_path, text):
    path = tmp_path / "pw.in"
    path.write_text(text, encoding="utf-8")
    result = import_native(path)
    assert result.spec is None
    assert result.unknown or result.unsupported


def lammps(tmp_path, data=DATA, extra=""):
    (tmp_path / "data.lammps").write_text(data, encoding="utf-8")
    (tmp_path / "in.lammps").write_text("units metal\natom_style atomic\nboundary p p p\nread_data data.lammps\npair_style lj/cut 2\npair_coeff * * 1 1\nrun 0\n" + extra, encoding="utf-8")
    return import_native(tmp_path)


def test_lammps_explicit_elements_and_readonly(tmp_path):
    result = lammps(tmp_path)
    assert result.spec is not None, result.report()
    assert result.spec.method.type_elements == ["Si"]
    assert result.spec.structure.atoms.symbols == ["Si"]
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    import_native(tmp_path)
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}


def test_lammps_no_mass_guessing_or_execution(tmp_path):
    assert lammps(tmp_path, data=DATA.replace(" # Si", "")).spec is None
    result = lammps(tmp_path, extra="shell touch should-not-exist\n")
    assert result.spec is None
    assert result.unknown[0]["line"] == 8
    assert not (tmp_path / "should-not-exist").exists()


def test_ambiguous_detection_requires_explicit_code(tmp_path):
    vasp(tmp_path)
    (tmp_path / "pw.in").write_text(QE, encoding="utf-8")
    with pytest.raises(NativeImportError):
        import_native(tmp_path)


def test_sparse_oversized_input_rejected(tmp_path):
    path = tmp_path / "pw.in"
    with path.open("wb") as stream:
        stream.truncate(MAX_BYTES + 1)
    result = import_native(path)
    assert result.spec is None
    assert result.files == {}


@pytest.mark.parametrize("change", ["Velocities\n\n1 0 0 0", "Bonds\n\n1 1 1 1"])
def test_lammps_extra_sections_block_spec(tmp_path, change):
    assert lammps(tmp_path, data=DATA + "\n" + change + "\n").spec is None


def test_lammps_huge_id_blocks_ase(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        pytest.fail("ASE called with unbounded atom ID")
    monkeypatch.setattr("ase.io.read", fail)
    assert lammps(tmp_path, data=DATA.replace("1 1 0 0 0", "999999999 1 0 0 0")).spec is None


def test_generated_vasp_single_point_mapped_fields_roundtrip(tmp_path):
    from adit.codes.vasp import VaspGenerator
    from adit.spec import AtomsData, CalculationSpec, KPoints, Structure, VaspMethod
    from ase import Atoms
    source = CalculationSpec(
        structure=Structure(source="file", source_ref="example", atoms=AtomsData.from_ase(Atoms("Si", cell=[3, 3, 3], pbc=True))),
        method=VaspMethod(encut=321, ediff=2e-7, nelm=77, ismear=-1, sigma=0.05, lasph=True),
        kpoints=KPoints(mode="mesh", mesh=(3, 2, 1)))
    gen = VaspGenerator()
    for name, text in {"INCAR": gen.incar(source, None), "POSCAR": gen.poscar(source), "KPOINTS": gen.kpoints(source)}.items():
        (tmp_path / name).write_text(text, encoding="utf-8")
    result = import_native(tmp_path)
    assert result.spec is not None, result.report()
    for name in ("encut", "ediff", "nelm", "ismear", "sigma", "lasph", "kpoints_centering"):
        assert getattr(result.spec.method, name) == getattr(source.method, name)
    assert result.spec.method.extra_incar["SYSTEM"] == "adit Si"
    assert result.provenance["structure.atoms.cell"]["line"] == 3
    assert result.provenance["structure.atoms.cell"]["end_line"] == 5
    assert result.provenance["structure.atoms.positions"]["line"] == 9
    assert "method.potcar_set" in result.inferred_defaults  # review, not physics equivalence


def test_generated_qe_single_point_mapped_fields_roundtrip(tmp_path):
    from adit.codes.espresso import EspressoGenerator
    from adit.spec import AtomsData, CalculationSpec, EspressoMethod, KPoints, Structure
    from ase import Atoms
    source = CalculationSpec(
        structure=Structure(source="file", source_ref="example", atoms=AtomsData.from_ase(Atoms("Si", cell=[3, 3, 3], pbc=True))),
        method=EspressoMethod(pseudo={"Si": "Si.UPF"}, ecutwfc=31, ecutrho=124, conv_thr=2e-8, electron_maxstep=77, mixing_beta=0.4),
        kpoints=KPoints(mode="mesh", mesh=(3, 2, 1), shift=(0.5, 0, 0.5)))
    text = EspressoGenerator().pw_in(source, None)
    (tmp_path / "pw.in").write_text(text, encoding="utf-8")
    result = import_native(tmp_path)
    assert result.spec is not None, result.report()
    for name in ("ecutwfc", "ecutrho", "conv_thr", "electron_maxstep", "mixing_beta", "pseudo"):
        assert getattr(result.spec.method, name) == getattr(source.method, name)
    assert result.spec.kpoints.mesh == source.kpoints.mesh
    assert result.spec.kpoints.shift == source.kpoints.shift
    assert result.spec.method.extra["control"] == {"tstress": True, "tprnfor": True}
    assert {x["field"] for x in result.not_applied} == {"control.prefix", "control.outdir", "control.pseudo_dir"}
    provenance = result.provenance["native.ATOMIC_POSITIONS"]
    assert text.splitlines()[provenance["line"] - 1] == "ATOMIC_POSITIONS crystal"
    assert provenance["end_line"] == provenance["line"] + 1
    assert "method.pseudo_set" in result.inferred_defaults


def test_qe_output_paths_only_reported_never_opened(tmp_path):
    text = QE.replace(" calculation = 'scf'", " calculation = 'scf'\n prefix = 'other'\n outdir = '/not/read'\n pseudo_dir = '/also/not/read'")
    (tmp_path / "pw.in").write_text(text, encoding="utf-8")
    result = import_native(tmp_path)
    assert result.spec is not None, result.report()
    assert len(result.files) == 1
    assert result.not_applied[1]["value"] == "/not/read"


def generated_lammps_md(tmp_path, units="metal", ensemble="NVT"):
    from adit.codes.lammps import LammpsGenerator
    from adit.spec import AtomsData, CalculationSpec, LammpsMethod, MDSettings, Structure, Task
    from ase import Atoms
    source = CalculationSpec(
        structure=Structure(source="file", source_ref="example", atoms=AtomsData.from_ase(Atoms("Ar2", positions=[[0, 0, 0], [2, 2, 2]], cell=[5, 5, 5], pbc=True))),
        method=LammpsMethod(units=units, pair_style="lj/cut 2", pair_coeff="* * 0.01 1", seed=4567),
        task=Task(type="molecular_dynamics", md=MDSettings(ensemble=ensemble, thermostat="nose_hoover", temperature_k=250,
                                                         timestep_fs=0.5, coupling_time_fs=75, steps=17, dump_interval=3)))
    gen = LammpsGenerator()
    (tmp_path / "in.lammps").write_text(gen.in_lammps(source), encoding="utf-8")
    (tmp_path / "data.lammps").write_text(gen.data(source), encoding="utf-8")
    return source


@pytest.mark.parametrize("units", ["metal", "real"])
@pytest.mark.parametrize("ensemble", ["NVE", "NVT"])
def test_generated_lammps_md_preserves_time_and_initialization(tmp_path, units, ensemble):
    source = generated_lammps_md(tmp_path, units, ensemble)
    result = import_native(tmp_path)
    assert result.spec is not None, result.report()
    assert result.spec.task.type == "molecular_dynamics"
    for name in ("ensemble", "temperature_k", "timestep_fs", "steps", "dump_interval"):
        assert getattr(result.spec.task.md, name) == getattr(source.task.md, name)
    if ensemble == "NVT":
        assert result.spec.task.md.thermostat == "nose_hoover"
        assert result.spec.task.md.coupling_time_fs == 75
    assert result.spec.method.seed == 4567
    assert result.spec.structure.velocities is None  # initialization command retained, not invented velocities
    assert result.provenance["task.md.timestep_fs"]["text"].startswith("timestep ")


@pytest.mark.parametrize("before,after", [
    ("rot no", "rot yes"), ("dist gaussian", "dist uniform"),
    ("nvt temp 250 250", "nvt temp 250 300"),
    ("fix integ all nvt temp 250 250 0.075", "fix integ all nvt temp 250 250 0.075 tchain 1"),
    ("run 17", "run 17 every 1 NULL"),
    ("velocity all create 250 4567 mom yes rot no dist gaussian\n", ""),
])
def test_lammps_md_unsupported_semantics_block(tmp_path, before, after):
    generated_lammps_md(tmp_path)
    path = tmp_path / "in.lammps"
    text = path.read_text(encoding="utf-8")
    assert before in text
    path.write_text(text.replace(before, after), encoding="utf-8")
    result = import_native(tmp_path)
    assert result.spec is None
    assert result.unsupported


def test_vasp_logical_for_a_string_tag_is_kept_as_text(tmp_path):
    result = vasp(tmp_path, incar="ENCUT=400\nLREAL=.FALSE.\nNSW=0\nIBRION=-1\n")
    assert result.spec is not None, result.report()
    assert result.spec.method.lreal == ".FALSE."
    assert result.provenance["method.lreal"]["line"] == 2


def test_vasp_value_outside_the_shared_fields_is_reported_by_field(tmp_path):
    result = vasp(tmp_path, incar="ENCUT=400\nISPIN=3\nNSW=0\nIBRION=-1\n")
    assert result.spec is None
    reasons = [issue["reason"] for issue in result.unsupported]
    assert any("共通の欄に入れられません" in reason for reason in reasons), reasons
    assert not any("形式と数値" in reason for reason in reasons)
    assert "ispin" in result.unsupported[-1]["detail"]


@pytest.mark.parametrize("text,expected", [
    ("'mv'", "marzari-vanderbilt"), ("'m-v'", "marzari-vanderbilt"), ("'cold'", "marzari-vanderbilt"), ("'Marzari-Vanderbilt'", "marzari-vanderbilt"),
    ("'gauss'", "gaussian"), ("'mp'", "methfessel-paxton"), ("'m-p'", "methfessel-paxton"), ("'fd'", "fermi-dirac"), ("'f-d'", "fermi-dirac")])
def test_qe_smearing_abbreviations_map_to_the_long_names(tmp_path, text, expected):
    path = tmp_path / "pw.in"
    path.write_text(QE.replace("ecutwfc = 30", f"ecutwfc = 30, occupations = 'smearing', smearing = {text}, degauss = 0.01"), encoding="utf-8")
    result = import_native(path)
    assert result.spec is not None, result.report()
    assert result.spec.method.smearing == expected and result.spec.method.occupations == "smearing"
    assert f"smearing = {text}" in result.provenance["method.smearing"]["text"]


def test_qe_tetrahedra_variants_are_reported_not_mapped(tmp_path):
    path = tmp_path / "pw.in"
    path.write_text(QE.replace("ecutwfc = 30", "ecutwfc = 30, occupations = 'tetrahedra_opt'"), encoding="utf-8")
    result = import_native(path)
    assert result.spec is None
    assert any("tetrahedra_opt" in issue["reason"] for issue in result.unsupported), result.unsupported


GRO = """water
    3
    1SOL      O    1   0.126   0.639   0.322
    1SOL      H    2   0.187   0.713   0.290
    1SOL      H    3   0.108   0.577   0.248
   3.00000   3.00000   3.00000
"""


def test_gromacs_mdp_keys_with_underscores_are_read(tmp_path):
    (tmp_path / "grompp.mdp").write_text("integrator = md\nnsteps = 100\ndt = 0.002\nnstxout_compressed = 10\n"
                                         "tcoupl = v-rescale\nref_t = 300\ntau_t = 0.1\ngen_seed = 7\n", encoding="utf-8")
    (tmp_path / "conf.gro").write_text(GRO, encoding="utf-8")
    (tmp_path / "topol.top").write_text("[ system ]\nx\n[ molecules ]\nSOL 1\n", encoding="utf-8")
    result = import_native(tmp_path)
    assert result.unknown == [] and result.spec is not None, result.report()
    md = result.spec.task.md
    assert md.temperature_k == 300 and md.coupling_time_fs == pytest.approx(100.0) and md.dump_interval == 10
    assert md.thermostat == "csvr" and result.spec.method.gen_seed == 7
    assert result.provenance["task.md.temperature_k"]["text"].startswith("ref_t")
