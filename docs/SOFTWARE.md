# Supported software

**This page lists how far ADIT supports each code.**
"Supported" means that ADIT can write the input files for that calculation.
The codes themselves, force fields, basis sets and pseudopotentials are not bundled (you provide them).
ADIT does not run or submit calculations either.

The support table in the README lists the supported calculation types. This page covers **the finer scope and limits that do not fit in that table**.

| Code | Supported task | Required `method` fields | Main input |
|---|---|---|---|
| Gaussian | Singlet molecular single point | `theory`, `basis` | `gaussian.gjf` |
| US GAMESS | Singlet molecular RHF single point | `gbasis` (`STO`, `N21`, `N31`, `N311`), `ngauss` | `gamess.inp` |
| Q-Chem | Singlet molecular single point | `theory` (`METHOD`), `basis` (`BASIS`) | `qchem.in` |
| GRRM17 | Singlet molecular MIN or FREQ | `theory`, `basis`; set `task.type` to `geometry_optimization` or `vibrations` | `grrm.com` |
| OpenMX | Neutral singlet, nonperiodic or fully periodic single point | `data_path`, per-element `pao`, `vps`, `valence`, `xc`, `energycutoff_ry` | `openmx.dat` |
| Amber (`sander`) | Nonperiodic or fully periodic energy minimization | `topology_file` (text `prmtop`), `coordinates_file` (text `rst7`/`inpcrd`), `cutoff_ang`; explicitly set `igb` for nonperiodic systems | `amber.in` plus supplied files |
| NAMD 3 | Nonperiodic or fully periodic NVE MD | `structure_file` (`.psf`), `coordinates_file` (`.pdb`), `parameter_files`, `exclude`, `one_four_scaling`, `cutoff_ang`, `pairlistdist_ang`; `switchdist_ang` when switching is enabled | `namd.conf` (with `paraTypeCharmm on`) plus supplied files |

Replace the existing `method` object in `spec.json` with one of the shapes below. The names and numbers **illustrate syntax; they are not recommended settings**. Check that the selected basis, pseudopotentials and force field are available where the calculation will run.

```json
{"code":"gaussian","theory":"HF","basis":"6-31G(d)"}
{"code":"gamess","gbasis":"STO","ngauss":3}
{"code":"qchem","theory":"HF","basis":"6-31G"}
{"code":"grrm","theory":"HF","basis":"6-31G"}
{"code":"openmx","data_path":"../DFT_DATA","pao":{"H":"H5.0-s1p1"},"vps":{"H":"H_CA19"},"valence":{"H":1},"xc":"LDA","energycutoff_ry":220}
{"code":"amber","topology_file":"system.prmtop","coordinates_file":"system.rst7","cutoff_ang":9,"igb":0}
{"code":"namd","structure_file":"system.psf","coordinates_file":"system.pdb","parameter_files":["force.prm"],"exclude":"scaled1-4","one_four_scaling":1,"cutoff_ang":10,"pairlistdist_ang":12,"switching":true,"switchdist_ang":8}
```

Each line is a separate example `method` object. The OpenMX example only shows the syntax for a system containing H. In a real calculation, supply `pao`, `vps` and `valence` for **every element in the structure**. Periodic systems also need `kpoints` in the spec, with zero shift. `DATA.PATH` must be readable **where the calculation runs**. ADIT does not bundle PAO/VPS data. The OpenMX initial spin splits the valence count given for each element equally between up and down, for a neutral, non-spin-polarized system. Magnetic and charged systems are not supported.

Amber and NAMD copy topology and coordinates prepared by the user. ADIT does not assign atom types, charges or force fields. It checks the atom count and coordinate order of the coordinate file against the spec, with a 1 mÅ tolerance. You must still check the PSF/prmtop atom types and force field, the Amber periodic cell, and the completeness of the NAMD parameters. Parameter files are read as CHARMM format (`namd.conf` sets `paraTypeCharmm on`; the NAMD default is X-PLOR format). Amber writes `task.max_steps` as `maxcyc` but does not translate the common `optimizer` or force threshold into Amber settings. NAMD writes the common MD `steps`, `temperature_k`, `timestep_fs` and `dump_interval`, and only runs NVE with initial velocities drawn from the given temperature.

This GRRM17 entry point covers only MIN/FREQ through Gaussian. Provide GRRM and a licensed Gaussian, and configure the link command as described in the official GRRM instructions. The common optimizer, maximum iteration count and force threshold are not translated into GRRM options. Using the same method and basis names across codes does not make the whole calculation setup, including implementation details and defaults, the same.

Conversion between codes records common settings that are not applied in `conversion.json` under `not_applied_by_target`. This covers the GRRM17 and Amber minimization settings, thermostat and pressure settings unused by NAMD NVE, and charge and spin multiplicity left to external Amber/NAMD topologies. Compare the converted input against the original conditions.

Without the executables and required licenses, verification is limited to the syntax of the generated files and mechanical checks. **Live runs and output analysis are not yet verified or implemented** for all seven codes, so their generated files do not include `analyze.py`. References for each input keyword and how to launch the code are in [SOURCES.md](../examples/new_engines/SOURCES.md). Materials Studio, Schrödinger Suite and Winmostar are not supported.

## NWChem: initial support

In `spec.json`, set `method` to, for example, `{"code":"nwchem","theory":"dft","basis":"cc-pvdz","xc":"b3lyp"}` and `task.type` to `single_point`. For SCF, set `theory` to `scf` and leave `xc` empty. These values illustrate syntax. They are **not ADIT recommendations**. `adit-gen spec.json out/` writes `nwchem.nw`, a `submit.sh` for submission, `spec.json`, `README.txt` and `analyze.py`, but it never submits or runs the job. For parallel runs, set the launch command explicitly in `commands.nwchem` in the preferences. After the run, `adit-analyze out/` reads only `Total SCF energy` or `Total DFT energy` from `output.log` and converts it from hartree to eV. It does not assess normal termination, geometry or frequencies. Check that the selected basis exists in the NWChem basis library where the calculation runs. Using the same method and basis names in NWChem and ORCA does not make their other default settings identical.

Input references: [geometry](https://nwchemgit.github.io/Geometry.html), [basis](https://nwchemgit.github.io/Basis.html), [charge](https://nwchemgit.github.io/Charge.html), [SCF](https://nwchemgit.github.io/Hartree-Fock-Theory-for-Molecules.html), [DFT](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html), [task](https://nwchemgit.github.io/TASK.html). Output label: [NWChem sample output](https://nwchemgit.github.io/Sample.html).

## OpenMM support

In `spec.json`, write `method` like this: `{"code":"openmm","input_format":"gromacs","topology_file":"topol.top","coordinates_file":"conf.gro","include_dir":"<GROMACS share/gromacs/top>","nonbonded_method":"PME","nonbonded_cutoff_nm":0.9,"constraints":"HBonds"}`. `input_format` is `gromacs` (`.top` plus `.gro`) or `amber` (`prmtop` plus `rst7`/`inpcrd`). Supported `task.type` values are `single_point`, `geometry_optimization` (`LocalEnergyMinimizer`) and `molecular_dynamics`. Vibrations and band structures are rejected. MD ensembles are NVE, NVT and NPT, and the thermostats are limited to those OpenMM provides: langevin (`LangevinMiddleIntegrator`), nose_hoover (`NoseHooverIntegrator`) and andersen (`AndersenThermostat`). berendsen and csvr are rejected because OpenMM has neither. **These values illustrate syntax; they are not ADIT recommendations.**

`adit-gen spec.json out/` writes `run_openmm.py`, `openmm_settings.json` and copies of the topology and coordinates. It never submits or runs the job. The environment that runs the calculation needs openmm (verified with 8.6.1 and 8.4.0.dev; ADIT itself never imports openmm). The nonbonded method, cutoff and constraints have no defaults and must be stated by the user, so the OpenMM default (`NoCutoff`) is never applied silently. The thermostat time constant is passed as a collision frequency (1 / time constant). NPT uses `MonteCarloBarostat` with a volume-move interval of `barostat_interval_steps` (the OpenMM default is 25 steps). The common barostat time constant has no mechanical conversion to that interval, so it is not used in the input, and the README says it is unused. Fixed atoms and axes, total charge and multiplicity, and continuation from checkpoints are not supported.

After the run, `adit-analyze out/` reads energies, temperature, volume and density from `md.log` (the OpenMM `StateDataReporter` CSV), and the final energy and temperature from `results.json`. The binary `trajectory.dcd` is not read (open it directly in VMD, OVITO or TRAVIS). When ADIT's RDF and MSD need a trajectory, enable `write_xyz_trajectory` to also write a text `trajectory.extxyz`. A verified run is in `examples/openmm_spce_nvt_generated` (openmm 8.4.0.dev, NVT of 216 SPC/E waters; the same input was also run with 8.6.1).

References: [running simulations](https://docs.openmm.org/latest/userguide/application/02_running_sims.html), [app layer API](https://docs.openmm.org/latest/api-python/app.html), [library API](https://docs.openmm.org/latest/api-python/library.html). Signatures were also verified against openmm 8.6.1.

## ABINIT support

Write `method` like this: `{"code":"abinit","pseudos":{"Si":"/…/14si.pspnc"},"ecut_ha":8,"tolerance":"toldfe","tolerance_value":1e-8}`, and also give `kpoints`. The supported range is only **single-point SCF and geometry optimization (fixed cell) for fully periodic systems**. Molecules (nonperiodic), MD, vibrations (anaddb), band structures, charged systems and spin polarization are rejected with a reason. The structure is written as `acell 1 1 1` with `rprim` and `xred` in bohr. Fixed atoms become `iatfix` and per-axis constraints become `iatfixx/y/z` (numbered from 1, as ABINIT expects). ABINIT interprets `iatfixx/y/z` along the lattice vectors (reduced coordinates), so Cartesian axis constraints are converted to lattice directions the same way as for VASP Selective dynamics. A cell that cannot be converted (when the fixed directions do not match the lattice vectors) is rejected. **Pseudopotentials are copied from your files. ADIT neither bundles nor validates them.** It selects no functional either (without `ixc`, ABINIT uses the functional recorded in the pseudopotential). `ecut`, the convergence criterion and its threshold, and the occupation scheme have no defaults. The criterion is one of `toldfe`, `toldff`, `tolrff`, `tolvrs` and `tolwfr` (ABINIT accepts exactly one). For a geometry optimization ABINIT rejects `toldfe`, so ADIT makes you choose `toldff`, `tolrff` or `tolvrs` (or `tolwfr`). This limit was found by running it, and the ABINIT documentation also says to avoid `toldfe` in geometry optimizations. Variables without a field on the screen can be passed as is through `extra` (a clash with a variable this generator writes is rejected).

Analysis reads the per-step `Total energy (etotal) [Ha]` from `output.log` (if absent, the last `etotal` in `input.abo`) and the final structure rebuilt from the variable table at the end of `input.abo`. **The ASE 3.29.0 `abinit-out` reader did not match the coordinates in this output, so it is not used.** A verified run is in `examples/abinit_si_generated` (ABINIT 10.0.3, diamond-structure Si).

## Psi4 support

Write `method` like this: `{"code":"psi4","method":"scf","basis":"cc-pVDZ"}`. Supported `task.type` values are `single_point` (`energy`), `geometry_optimization` (`optimize`) and `vibrations` (`frequencies`). Periodic systems, MD, band structures and fixed atoms are rejected. The input is psithon (the Psi4 input language), with the atoms in ADIT's order (`no_reorient`, `no_com`, `symmetry c1`). A multiplicity other than 1 requires an explicit open-shell reference (`uhf`, `rohf`, `uks`), because the Psi4 default RHF is closed-shell only. A mismatch between electron count and multiplicity is caught by the same shared check as for the other codes. The last few lines of the input write `results.json` and `final.xyz`. Analysis reads those plus the "Optimization Summary" table (energy per step) in `output.log`. Convergence of a geometry optimization is decided from the fact that Psi4's `optimize()` stops with an exception when it does not converge (`converged` in `results.json` means it converged). A verified run is in `examples/psi4_h2o_generated` (Psi4 1.11, RHF/cc-pVDZ optimization of water).

## PLUMED support

PLUMED is not a standalone code, so it goes in the spec's `plumed` field rather than in `method`: `{"plumed":{"lines":"d1: DISTANCE ATOMS=1,2\nPRINT ARG=d1 FILE=COLVAR STRIDE=10"}}` (the PLUMED input as text) or `{"plumed":{"input_file":"/…/plumed.dat"}}` (a file to copy), one or the other, never both. **ADIT creates no collective variables or biases.** It copies your input to the output as `plumed.dat` and adds one line to the MD input that reads it. Only LAMMPS, GROMACS and OpenMM MD runs can be connected. Other codes and non-MD calculations are rejected.

Every connection was verified by running it: `fix adit_plumed all plumed plumedfile plumed.dat outfile plumed.log` for LAMMPS, `gmx mdrun … -plumed plumed.dat` for GROMACS, and `PlumedForce` from `openmm-plumed` for OpenMM. **With GROMACS 2026.3, `-plumed` alone was not enough: the run stopped unless the environment variable `PLUMED_KERNEL` pointed at `libplumedKernel.so`.** README.txt includes this fact. When `plumed` is on PATH, ADIT has it read the input with `plumed driver --parse-only` before generating, and stops with the reason if the input cannot be read (it never comments on the physics). It also checks `lmp -h` / `gmx mdrun -h` to see whether that executable can use PLUMED. When this cannot be determined, it shows nothing and does not stop. The PLUMED default units are nm, kJ/mol and ps, separate from the unit system of the code (a verified run shows Cu's 2.556 Å written as `0.2556`). Analysis looks for COLVAR-format files that start with `#! FIELDS` and adds the column names, row count, and first, last and mean values to the summary (it does not interpret what the columns mean). A verified run is in `examples/plumed_lammps_cu_generated`.


## Additional entry points (ORCA GOAT, DOCKER and others)

### ORCA GOAT

Set `method.code="orca"`, `task.type="single_point"` and `method.goat=true` in `spec.json`, then generate with `adit-gen spec.json out/`. `method.method="XTB"` and `method.basis=""` are an input example matching the [ORCA 6.1 GOAT tutorial](https://www.faccts.de/docs/orca/6.1/tutorials/prop/goat.html). After the run, check `orca.globalminimum.xyz`, `orca.finalensemble.xyz`, and the energies and weights in `output.log`. ORCA itself is not bundled, and running this input with ORCA has not been verified.

### ORCA DOCKER / Host–guest docking

The host is the structure in `spec.json`, and the guest is the XYZ in `method.docker_guest_file`. Use `task.type="single_point"`, `method.method="XTB"` and `method.basis=""`, and write the charge and multiplicity as integers on the second line of each XYZ frame (for example `0 1`). Only if you intend the ORCA default (0, 1) for frames without these values, explicitly set `method.docker_assume_neutral_singlet=true`. ADIT copies the guest as `guest.xyz` and writes `%DOCKER GUEST "guest.xyz" END`. Following the [ORCA 6.1 manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/docker.html), DOCKER is limited to XTB/GFN-xTB/GFN-FF. Check `orca.docker.xyz` and `output.log` for the results.

### DCDFTBMD 2.0

In addition to `method.code="dcdftbmd"`, explicitly set `method.scc`, `method.divide_and_conquer`, `method.highest_angular_momentum` for each element in the structure (s=1, p=2, d=3, f=4), and `method.sk_files` (`.spl`) for every **ordered** element pair. For O and H, for example, you need `O-O`, `O-H`, `H-O` and `H-H`. `adit-gen spec.json out/` creates `dftb.inp`, copies of the `.spl` files, `submit.sh`, `spec.json` and `README.txt`. It does not convert DFTB+ `.skf` files. Supported tasks are single point, steepest-descent/FIRE geometry optimization, NVE and Berendsen NVT. NVT temperature and time constant come from the common MD settings, with only the units converted. NPT, open-shell systems, axis constraints and other thermostats stop before generation. Obtain the executable and `.spl` files yourself. See the [official 2.0 manual](https://www.chem.waseda.ac.jp/dcdftbmd/document/DCDFTBMD_2.0_en.pdf).

### DOCK6

`adit-convert dock6 prepared/dock.in ready/` packages a **completed, user-authored** `dock.in` with the ligand, spheres, grid and recognized definition files it refers to. Repeat `--asset relative/path` to add others as needed. Paths must be relative to the working directory of `dock.in`, and the output must be an empty location outside the input directory. The original `dock.in` is not changed. ADIT writes the SHA-256 of the copies and `run.sh`. This does not prepare the receptor or ligand, assign charges or atom types, or create scoring settings. Complete the input using the [DOCK 6.13 manual](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm) before using it.

### Open Babel

`adit-convert openbabel source.sdf target.mol2 --input-format sdf --output-format mol2` explicitly uses `obabel` from PATH. Input and output cannot be the same file, and an existing output is not changed unless you pass `--overwrite`. ADIT writes the result to a temporary file and then moves it into place. On failure, it keeps the previous output and shows where the partial file is. Check that bond orders, charges, cell and coordinates are what you intended. See the [Open Babel CLI documentation](https://openbabel.org/docs/Command-line_tools/babel.html).

## Comparing settings across codes

`adit-convert calculation` adds `field_audit` to the conversion record `conversion.json`.
It records the **origin and applicability** of each value. It does not certify
that different codes give the same physical model or numerical solution.

| status | Meaning |
| --- | --- |
| `preserved` | Value preserved from the source spec |
| `changed_by_request` | Changed to a target-template value at the user's request |
| `omitted_by_request` | Omitted at the user's request |
| `retained_in_spec_not_input` | Retained in the spec but not applied to the generated input |
| `not_used_by_task` | Not used by this task type |
| `from_target_template` | Taken from the target template without translation |
| `not_translated` | Code-specific method settings are not mapped |

Large arrays such as coordinates, atomic symbols and cell vectors are not copied into
`field_audit`. You can check the exact values in the generated `spec.json` and structure files.
This record does not guarantee implicit defaults not written in the input, physical
equivalence of force fields or pseudopotentials, or equivalence of thermostat algorithms.

## Mechanical checks of results

```text
adit-analyze run_lammps --audit-with run_vasp --audit-kind msd
adit-analyze run_lammps --audit-with run_vasp --audit-kind energy
```

This command only reports the check on standard output. It does not produce MSD or energy plots. It cannot be combined with ordinary analysis options such as `--msd`. If you need plots, analyze each calculation separately.

For MSD, the check covers the element-symbol sequence, periodic boundaries, the cell, trajectories and frame spacing. Symbols alone cannot tell atoms of the same element apart, so initial positions are also matched, with periodic boundaries taken into account. If they differ, the check for a comparison from the same initial structure fails. Whether a statistical comparison from different initial structures is valid is not judged.
For energy, the check covers matching composition, periodic boundaries and cell, and whether the energies can be read in eV for the whole system. Missing trajectories, frame counts and atom order do not stop this mode. The output is read-only,
and a failed check returns exit code 1. Whether energies or MSDs from different methods
are **scientifically comparable** is not judged.
For energies from different codes, a direct comparison fails the check, because the breakdown their output includes (kinetic, potential and electronic contributions, for example) has not been matched.

## Scanning output logs for warnings and failure signs

The analysis summary (`adit-analyze`, the Analysis tab, the web analysis page, `tables.diagnostics` in `analysis/summary.json`) includes the number of warnings in the output log and candidate causes of failure.
ADIT searches only for strings each code actually writes. The covered codes are DFTB+, VASP, Quantum ESPRESSO (pw.x), xtb, ORCA, CP2K, LAMMPS and GROMACS.
For other codes, it searches only for signs from the job scheduler (PBS / Slurm), MPI and the OS (time limit, out of memory, abnormal MPI termination).
It copies the matching lines and classifies them. It does not decide on remedies. The remedies given in each code's documentation are listed, with sources, at the end of the generated `README.txt`. For details, see [Analysis details](ANALYSIS.md).
