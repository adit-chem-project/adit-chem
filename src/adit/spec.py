"""CalculationSpec: the settings shared by the user interfaces and the generators. JSON round-trips; structures are ASE Atoms."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Literal

import numpy as np
from ase import Atoms
from ase import units as _ase_units
from pydantic import BaseModel, Field, field_validator, model_validator

from adit import __version__

Vec3 = tuple[float, float, float]


class Meta(BaseModel):
    created: str = Field(default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds"))
    app_version: str = __version__
    comment: str = ""
    continued_from: dict | None = None
    template: dict | None = None
    stage: dict | None = None


class AtomsData(BaseModel):

    symbols: list[str]
    positions: list[Vec3]
    cell: list[Vec3] = Field(default_factory=lambda: [(0.0, 0.0, 0.0)] * 3)
    pbc: tuple[bool, bool, bool] = (False, False, False)

    @model_validator(mode="after")
    def _same_length(self) -> "AtomsData":
        if len(self.symbols) != len(self.positions):
            from adit.lang import L
            raise ValueError(L(f"元素記号の数 ({len(self.symbols)}) と座標の数 ({len(self.positions)}) が合いません",
                               f"there are {len(self.symbols)} chemical symbols but {len(self.positions)} positions"))
        if len(self.cell) != 3:
            from adit.lang import L
            raise ValueError(L("セルは 3 本の格子ベクトルで書いてください", "the cell needs three lattice vectors"))
        return self

    @classmethod
    def from_ase(cls, atoms: Atoms) -> "AtomsData":
        return cls(
            symbols=list(atoms.get_chemical_symbols()),
            positions=[tuple(float(x) for x in p) for p in atoms.get_positions()],
            cell=[tuple(float(x) for x in v) for v in np.asarray(atoms.cell)],
            pbc=tuple(bool(b) for b in atoms.pbc),
        )

    def to_ase(self) -> Atoms:
        return Atoms(symbols=self.symbols, positions=self.positions, cell=self.cell, pbc=self.pbc)


SPEC_VERSION = 2
HARTREE_PER_BOHR_IN_EV_PER_ANG = _ase_units.Hartree / _ase_units.Bohr  # 1 Hartree/Bohr [eV/Å] ≈ 51.42207
RY_PER_BOHR_IN_EV_PER_ANG = _ase_units.Rydberg / _ase_units.Bohr       # 1 Ry/Bohr [eV/Å] ≈ 25.71103
HARTREE_PER_BOHR3_IN_GPA = _ase_units.Hartree / _ase_units.Bohr ** 3 / _ase_units.GPa  # 1 Hartree/Bohr³ [GPa] ≈ 29421.0
EV_PER_ANG3_IN_GPA = 1.0 / _ase_units.GPa                              # 1 eV/Å³ [GPa] ≈ 160.2177


class Structure(BaseModel):
    source: Literal["preset", "smiles", "file", "bulk", "spacegroup", "surface", "mixture", "recipe"]
    source_ref: str
    atoms: AtomsData
    charge: int = 0
    multiplicity: int = 1
    fixed_atoms: list[int] = Field(default_factory=list)
    fixed_axes: dict[str, tuple[bool, bool, bool]] = Field(default_factory=dict)
    velocities: list[Vec3] | None = None
    fetched: dict | None = None  # provenance of a structure fetched from a database (adit.fetch)

    @property
    def periodic(self) -> bool:
        return any(self.atoms.pbc)

    @field_validator("multiplicity")
    @classmethod
    def _positive(cls, v: int) -> int:
        if v < 1:
            from adit.lang import L
            raise ValueError(L("スピン多重度は 1 以上にしてください", "the spin multiplicity must be at least 1"))
        return v


class KPoints(BaseModel):

    mode: Literal["gamma", "mesh", "density"] = "gamma"
    mesh: tuple[int, int, int] = (1, 1, 1)  # mode = mesh
    shift: tuple[float, float, float] = (0.0, 0.0, 0.0)
    density: float = 0.0

    def resolved_mesh(self, cell) -> tuple[int, int, int]:
        if self.mode == "gamma":
            return (1, 1, 1)
        if self.mode == "mesh":
            return self.mesh
        rec = 2 * np.pi * np.linalg.inv(np.asarray(cell, dtype=float)).T
        lengths = np.linalg.norm(rec, axis=1)
        return tuple(max(1, int(np.ceil(self.density * L))) for L in lengths)


class HubbardU(BaseModel):

    orbital: str
    u_ev: float   # U [eV]
    j_ev: float = 0.0


class DftbMethod(BaseModel):
    code: Literal["dftbplus"] = "dftbplus"
    sk_set: str
    scc: bool = True
    scc_tolerance: float = 1e-5
    max_scc_iterations: int = 100
    third_order: bool = False
    dispersion: Literal["none", "dftd3", "lennard-jones"] = "none"
    d3_params: dict[str, float] | None = None
    filling_temperature: float = 0.0  # K
    seed: int = 12345
    solvation_param_file: str = ""


class VaspMethod(BaseModel):
    code: Literal["vasp"] = "vasp"
    potcar_set: str = "potpaw_PBE"
    potcar: dict[str, str] = Field(default_factory=dict)
    binary: Literal["std", "gam", "ncl"] = "std"
    encut: float = 0.0
    ediff: float = 1e-4  # eV
    nelm: int = 60
    ismear: int = 0
    sigma: float = 0.1
    ispin: Literal[1, 2] = 1
    magmom: list[float] | None = None
    ivdw: int | None = None
    ibrion: int = 2
    algo: str = "Normal"
    prec: str = "Normal"
    lreal: str = "Auto"
    nelmin: int = 0
    lasph: bool = False
    lmaxmix: int = 0
    nbands: int = 0
    isym: int | None = None
    idipol: Literal[0, 1, 2, 3, 4] = 0
    ldipol: bool = False
    dipol: str = ""
    kpoints_centering: Literal["monkhorst-pack", "gamma"] = "monkhorst-pack"
    magmom_by_element: dict[str, float] = Field(default_factory=dict)
    hubbard: dict[str, HubbardU] = Field(default_factory=dict)
    ldau_type: Literal[1, 2, 4] = 2
    extra_incar: dict[str, str | float | int | bool] = Field(default_factory=dict)


class XtbMethod(BaseModel):
    code: Literal["xtb"] = "xtb"
    gfn: Literal["0", "1", "2", "ff"] = "2"
    accuracy: float = 1.0
    etemp: float = 300.0
    max_iterations: int = 250
    opt_level: Literal["crude", "sloppy", "loose", "normal", "tight", "verytight", "extreme"] = "normal"
    solvation: Literal["none", "alpb", "gbsa"] = "none"
    solvent: str = ""
    md_hmass: float = 4.0
    md_shake: Literal[0, 1, 2] = 2
    md_sccacc: float = 2.0


class EspressoMethod(BaseModel):
    code: Literal["espresso"] = "espresso"
    pseudo_set: str = ""
    pseudo: dict[str, str] = Field(default_factory=dict)
    ecutwfc: float = 0.0
    ecutrho: float = 0.0
    conv_thr: float = 1e-6  # Ry
    electron_maxstep: int = 100
    mixing_beta: float = 0.7
    occupations: Literal["fixed", "smearing", "tetrahedra"] = "fixed"
    smearing: Literal["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"] = "gaussian"
    degauss: float = 0.0  # Ry
    nspin: Literal[1, 2] = 1
    input_dft: str = ""
    starting_magnetization: dict[str, float] = Field(default_factory=dict)
    hubbard: dict[str, HubbardU] = Field(default_factory=dict)
    hubbard_projector: Literal["atomic", "ortho-atomic", "norm-atomic", "wf", "pseudo"] = "atomic"
    assume_isolated: Literal["", "makov-payne", "martyna-tuckerman", "esm", "2D"] = ""
    dipole_correction: bool = False
    dipole_direction: Literal[0, 1, 2, 3] = 0
    dipole_maxpos: float = 0.0
    dipole_decrease: float = 0.0
    dipole_amplitude: float = 0.0
    extra: dict[str, dict[str, str | float | int | bool]] = Field(default_factory=dict)


class OrcaMethod(BaseModel):
    code: Literal["orca"] = "orca"
    method: str = "HF"
    basis: str = "def2-SVP"
    scf_convergence: Literal["LooseSCF", "NormalSCF", "TightSCF"] = "NormalSCF"
    scf_maxiter: int = 125
    maxcore_mb: int = 0
    extra_keywords: str = ""
    extra_blocks: str = ""
    solvation: Literal["none", "cpcm", "smd"] = "none"
    solvent: str = ""
    ts_search: bool = False
    ts_calc_hess: bool = True
    ts_recalc_hess: int = 0
    ts_freq: bool = True
    irc: bool = False
    irc_max_iter: int = 0
    irc_direction: Literal["both", "forward", "backward", "down"] = "both"
    goat: bool = False
    docker_guest_file: str = ""
    docker_assume_neutral_singlet: bool = False


class DcdftbmdMethod(BaseModel):

    code: Literal["dcdftbmd"] = "dcdftbmd"
    scc: bool | None = None
    divide_and_conquer: bool | None = None
    highest_angular_momentum: dict[str, int] = Field(default_factory=dict)  # s=1, p=2, d=3, f=4
    sk_files: dict[str, str] = Field(default_factory=dict)


class NwchemMethod(BaseModel):

    code: Literal["nwchem"] = "nwchem"
    theory: Literal["scf", "dft"] | None = None
    basis: str = ""
    xc: str = ""


class GaussianMethod(BaseModel):

    code: Literal["gaussian"] = "gaussian"
    theory: str = ""
    basis: str = ""


class GamessMethod(BaseModel):

    code: Literal["gamess"] = "gamess"
    gbasis: str = ""
    ngauss: int = 0


class QchemMethod(BaseModel):

    code: Literal["qchem"] = "qchem"
    theory: str = ""  # $rem METHOD
    basis: str = ""   # $rem BASIS


class GrrmMethod(BaseModel):

    code: Literal["grrm"] = "grrm"
    theory: str = ""
    basis: str = ""


class OpenmxMethod(BaseModel):

    code: Literal["openmx"] = "openmx"
    data_path: str = ""
    pao: dict[str, str] = Field(default_factory=dict)
    vps: dict[str, str] = Field(default_factory=dict)
    valence: dict[str, float] = Field(default_factory=dict)
    xc: Literal["", "LDA", "GGA-PBE"] = ""
    energycutoff_ry: float = 0.0


class AmberMethod(BaseModel):

    code: Literal["amber"] = "amber"
    topology_file: str = ""
    coordinates_file: str = ""
    cutoff_ang: float = 0.0
    igb: int | None = None


class NamdMethod(BaseModel):

    code: Literal["namd"] = "namd"
    structure_file: str = ""  # .psf
    coordinates_file: str = ""  # .pdb
    parameter_files: list[str] = Field(default_factory=list)
    exclude: Literal["", "none", "1-2", "1-3", "1-4", "scaled1-4"] = ""
    one_four_scaling: float | None = None
    cutoff_ang: float = 0.0
    pairlistdist_ang: float = 0.0
    switching: bool | None = None
    switchdist_ang: float = 0.0
    seed: int = 12345


class Psi4Method(BaseModel):

    code: Literal["psi4"] = "psi4"
    method: str = ""
    basis: str = ""
    reference: Literal["", "rhf", "uhf", "rohf", "rks", "uks"] = ""
    memory_mb: int = 0
    extra_set: dict[str, str] = Field(default_factory=dict)


class AbinitMethod(BaseModel):

    code: Literal["abinit"] = "abinit"
    pseudos: dict[str, str] = Field(default_factory=dict)
    ecut_ha: float = 0.0
    pawecutdg_ha: float = 0.0
    tolerance: Literal["", "toldfe", "toldff", "tolrff", "tolvrs", "tolwfr"] = ""
    tolerance_value: float = 0.0
    nstep: int = 30
    ixc: int | None = None
    occopt: int | None = None
    tsmear_ha: float = 0.0
    nband: int | None = None
    extra: dict[str, str] = Field(default_factory=dict)


class OpenmmMethod(BaseModel):

    code: Literal["openmm"] = "openmm"
    input_format: Literal["", "gromacs", "amber"] = ""
    topology_file: str = ""
    coordinates_file: str = ""
    include_dir: str = ""
    defines: dict[str, str] = Field(default_factory=dict)
    nonbonded_method: Literal["", "NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME", "LJPME"] = ""
    nonbonded_cutoff_nm: float = 0.0
    constraints: Literal["none", "HBonds", "AllBonds", "HAngles"] = "none"
    rigid_water: bool = True
    barostat_interval_steps: int = 25
    platform: str = ""
    write_xyz_trajectory: bool = False
    seed: int = 12345


class Cp2kMethod(BaseModel):

    code: Literal["cp2k"] = "cp2k"
    basis_file: str = "BASIS_MOLOPT"
    potential_file: str = "GTH_POTENTIALS"
    basis: dict[str, str] = Field(default_factory=dict)
    potential: dict[str, str] = Field(default_factory=dict)
    xc: str = ""
    dispersion: Literal["none", "d3", "d3bj"] = "none"
    cutoff_ry: float = 0.0
    rel_cutoff_ry: float = 0.0
    eps_scf: float = 1e-5
    max_scf: int = 50
    uks: bool = False
    poisson_solver: Literal["", "MT", "WAVELET", "ANALYTIC", "MULTIPOLE"] = ""
    isolated_box_ang: float = 0.0
    magnetization_by_element: dict[str, float] = Field(default_factory=dict)
    hubbard: dict[str, HubbardU] = Field(default_factory=dict)
    plus_u_method: Literal["MULLIKEN", "LOWDIN", "MULLIKEN_CHARGES"] = "MULLIKEN"
    sccs_relative_permittivity: float = 0.0
    ot: bool = False
    ot_minimizer: Literal["", "SD", "CG", "DIIS", "BROYDEN", "LBFGS"] = ""
    ot_preconditioner: Literal["", "FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", "NONE"] = ""
    ot_extra: str = ""
    surface_dipole_correction: bool = False
    surf_dip_dir: Literal["", "X", "Y", "Z"] = ""
    extra_sections: dict[str, str] = Field(default_factory=dict)


class LammpsMethod(BaseModel):

    code: Literal["lammps"] = "lammps"
    units: Literal["", "metal", "real"] = ""
    atom_style: str = "atomic"
    data_file: str = ""
    type_elements: list[str] = Field(default_factory=list)
    pair_style: str = ""
    pair_coeff: str = ""
    potential_files: list[str] = Field(default_factory=list)
    style_commands: str = ""
    extra_commands: str = ""
    seed: int = 12345
    thermo_pressure_tensor: bool = False


class GromacsMethod(BaseModel):

    code: Literal["gromacs"] = "gromacs"
    topology_file: str = ""
    structure_file: str = ""
    coulombtype: Literal["Cut-off", "PME", "Reaction-Field"] = "Cut-off"
    rcoulomb_nm: float = 1.0
    rvdw_nm: float = 1.0
    constraints: Literal["none", "h-bonds", "all-bonds", "h-angles", "all-angles"] = "none"
    pcoupl: Literal["", "C-rescale", "Parrinello-Rahman"] = ""
    compressibility_per_bar: float = 0.0
    define: str = ""
    checkpoint_file: str = ""
    gen_seed: int = 12345
    extra_mdp: dict[str, str | float | int | bool] = Field(default_factory=dict)


class MlipMethod(BaseModel):

    code: Literal["mlip"] = "mlip"
    model_family: Literal["", "mace_mp", "mace_off", "chgnet"] = ""
    model: str = ""
    device: str = ""
    dtype: Literal["", "float32", "float64"] = ""
    dispersion: bool = False
    seed: int = 12345


Method = Annotated[DftbMethod | VaspMethod | XtbMethod | EspressoMethod | OrcaMethod | DcdftbmdMethod | NwchemMethod | GaussianMethod | GamessMethod | QchemMethod | GrrmMethod | OpenmxMethod | AmberMethod | NamdMethod | AbinitMethod | Psi4Method | OpenmmMethod | Cp2kMethod | LammpsMethod | GromacsMethod | MlipMethod,
                   Field(discriminator="code")]


class MDSettings(BaseModel):

    ensemble: Literal["NVE", "NVT", "NPT"] = "NVT"
    thermostat: Literal["berendsen", "andersen", "nose_hoover", "langevin", "csvr"] = "berendsen"
    temperature_k: float = 300.0
    timestep_fs: float = 1.0
    steps: int = 1000
    dump_interval: int = 10
    coupling_time_fs: float = 100.0
    pressure_bar: float = 1.0
    barostat_time_fs: float = 1000.0


class BandSettings(BaseModel):

    path: str = ""
    npoints: int = 60
    empty_bands: int = 4


class Task(BaseModel):
    type: Literal["single_point", "geometry_optimization", "molecular_dynamics", "vibrations", "band_structure"] = "single_point"
    bands: BandSettings = Field(default_factory=BandSettings)
    optimizer: Literal["Rational", "LBFGS", "FIRE", "SteepestDescent"] = "Rational"
    md: MDSettings = Field(default_factory=MDSettings)
    max_steps: int = 200
    force_tolerance_ev_per_ang: float = 1e-4 * HARTREE_PER_BOHR_IN_EV_PER_ANG
    relax_cell: Literal["no", "shape_and_volume", "volume_only"] = "no"

    @property
    def force_tolerance(self) -> float:
        return self.force_tolerance_ev_per_ang


class Runtime(BaseModel):
    profile: str = "local"
    nodes: int = 1
    ncpus: int = 1
    mpiprocs: int = 1
    omp_threads: int = 1
    walltime: str = "01:00:00"
    job_name: str = "adit"


TaskType = Literal["single_point", "geometry_optimization", "molecular_dynamics", "vibrations", "band_structure"]


class Handoff(BaseModel):

    previous_dir: str
    previous_task: TaskType
    previous_code: str
    at_run: bool = False
    velocities: bool = False
    files: dict[str, str] = Field(default_factory=dict)


class PlumedSettings(BaseModel):

    input_file: str = ""
    lines: str = ""


class CalculationSpec(BaseModel):
    version: int = SPEC_VERSION
    meta: Meta = Field(default_factory=Meta)
    structure: Structure
    method: Method
    kpoints: KPoints | None = None
    task: Task = Field(default_factory=Task)
    runtime: Runtime = Field(default_factory=Runtime)
    handoff: Handoff | None = None
    plumed: PlumedSettings | None = None

    @staticmethod
    def migrate(data: dict) -> dict:
        data = dict(data)
        version = data.get("version", 1)
        if version < 2:
            task = dict(data.get("task", {}))
            if "force_tolerance" in task:
                task["force_tolerance_ev_per_ang"] = float(task.pop("force_tolerance")) * HARTREE_PER_BOHR_IN_EV_PER_ANG
            data["task"] = task
            data["version"] = 2
        return data

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.model_dump(mode="json"), indent=indent, ensure_ascii=False)

    @classmethod
    def from_json(cls, text: str) -> "CalculationSpec":
        return cls.model_validate(cls.migrate(json.loads(text)))

    def save(self, path: Path | str) -> None:
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path | str) -> "CalculationSpec":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))

    @staticmethod
    def describe_error(ex: Exception) -> str:
        import json as _json

        from pydantic import ValidationError as _ValidationError

        from adit.lang import L

        if isinstance(ex, _json.JSONDecodeError):
            return L(f"spec.json が JSON として壊れています ({ex.lineno} 行目 {ex.colno} 文字目: {ex.msg})。"
                     "括弧・引用符・カンマを確かめてください。",
                     f"spec.json is not valid JSON (line {ex.lineno}, column {ex.colno}: {ex.msg}); "
                     "check the brackets, quotes and commas.")
        if isinstance(ex, FileNotFoundError):
            name = getattr(ex, "filename", "") or ""
            return L(f"ファイルがありません: {name}。条件を書いた JSON のファイル名を確かめてください "
                     "(サンプルは `adit-gen --list-samples` / `--sample <名前> <出力先>` で作れます)。",
                     f"no such file: {name}. Check the name of the JSON file with the settings "
                     "(`adit-gen --list-samples` / `--sample <name> <destination>` writes a sample you can start from).")
        if isinstance(ex, IsADirectoryError):
            name = getattr(ex, "filename", "") or ""
            return L(f"{name} はディレクトリです。条件を書いた JSON のファイル (spec.json) を指してください。",
                     f"{name} is a directory; point at the JSON file with the settings (spec.json).")
        if isinstance(ex, PermissionError):
            name = getattr(ex, "filename", "") or ""
            if name and Path(name).is_dir():
                return L(f"{name} はディレクトリです。条件を書いた JSON のファイル (spec.json) を指してください。",
                         f"{name} is a directory; point at the JSON file with the settings (spec.json).")
            return L(f"{name} を読む権限がありません。", f"no permission to read {name}.")
        if not isinstance(ex, _ValidationError):
            return str(ex)
        from adit.validate_types import friendly_pydantic

        head = L("spec.json の次の欄を直してください:", "fix the following fields in spec.json:")
        return head + "\n" + "\n".join("  " + line for line in friendly_pydantic(ex).splitlines())


    @property
    def atoms(self) -> Atoms:
        return self.structure.atoms.to_ase()

    @property
    def elements(self) -> list[str]:
        seen: list[str] = []
        for s in self.structure.atoms.symbols:
            if s not in seen:
                seen.append(s)
        return seen
