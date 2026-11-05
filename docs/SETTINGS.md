# Settings

This document describes the file that sets where the calculation codes are and where to run them (profiles).
Analysis conditions are in [Analysis details](ANALYSIS.md).

## Settings file (cluster.toml)

```toml
sk_root = "/home/<user name>/slakos"   # location of the Slater-Koster parameters (with mio-1-1/ and so on directly inside)
pseudo_root = "/home/<user name>/pseudo" # Quantum ESPRESSO pseudopotentials (with <set name>/*.UPF directly inside)
cp2k_data = ""                       # CP2K data directory (BASIS_MOLOPT and so on). If empty, ADIT looks at CP2K_DATA_DIR, then at share/cp2k/data next to cp2k
templates_dir = ""                   # location of your group's templates (spec.json without the structure). If empty, ADIT only looks at templates/ next to this file
mp_api_key = ""                      # Materials Project API key (only used to fetch structures by mp-ID; issued at https://next-gen.materialsproject.org/dashboard)
default_profile = "local"
enable_run = true                    # false means this PC does not run calculations (generation and analysis still work)
language = "en"                      # "en" (default) or "ja" (the environment variable ADIT_LANG takes precedence)
theme = "auto"                       # "auto" / "light" / "dark" (the environment variable ADIT_THEME takes precedence)
window_frame = "auto"                # window frame: "auto" (ADIT draws it on Linux) / "custom" (ADIT draws it) / "native" (left to the OS). The environment variable ADIT_FRAME takes precedence
ribbon_collapsed = true              # whether to start with the ribbon collapsed. Changing it with the ˄ / ˅ button on screen or Ctrl+F1 saves it here

[profiles.local]
kind = "direct"
description = "run bash submit.sh on this PC"
```

### Environment variables

These take precedence over the settings file. Use them to switch something temporarily.

| Environment variable | Meaning |
|---|---|
| `ADIT_CONFIG` | Path of the settings file to read (default `~/.config/adit/cluster.toml`; on Windows `%APPDATA%\adit\cluster.toml`) |
| `ADIT_LANG` | Display language `ja` / `en` |
| `ADIT_THEME` | Screen theme `auto` / `light` / `dark` |
| `ADIT_FRAME` | Window frame `auto` / `custom` / `native` |
| `ADIT_SHELL` | Shell used in the Workspace terminal (for example `ADIT_SHELL="bash --norc"`) |

### Running on a cluster

If you add a profile with `kind = "pbs"` or `kind = "slurm"` to the settings file, `submit.sh` gets that scheduler's header.
The template writes only the generic items (job name, node count, core count, MPI process count, OpenMP thread count, time limit).
All site-specific values come from the settings (empty by default).
The `<...>` values in the example below will not work as they are, so check them with your cluster administrator or a senior member of your group and replace them. If you miss one, the validation before generation reports "a placeholder is still in the settings".
Keep lines you don't use disabled with `#` at the start.

If you set `host` and `remote_dir`, the generated `transfer_and_submit.sh` contains `rsync`, `ssh` and submit commands ready to paste (**ADIT does not run them**).
If you set `cores_max` / `nodes_max` / `walltime_max`, generation stops beforehand when the screen or the runtime in `spec.json` exceeds them, in the same way as "a placeholder is still in the settings" (the default is none, meaning no check).

```toml
[profiles.remote]
kind = "pbs"                        # or "slurm"
description = "your institution's cluster"
select_extra = ""                   # string added to the end of the PBS select line (if your site requires it; for example ":jobtype=core")
header_extra = ["#PBS -q <queue>"]  # lines added to the header (queue / partition / account and so on; for Slurm, "#SBATCH --partition=<partition>")
submit_command = "qsub"             # name of the submit command (only written in README.txt; ADIT does not run it)
status_command = "qstat -u $USER"   # command that checks status (same as above)
host = "<cluster host name>"        # used in the transfer and submit commands (generation also works if empty; the commands then keep <...>)
user = "<login name>"               # if empty, your local login name is assumed; write only host
remote_dir = "<working directory on the cluster>"  # transfer destination
cores_max = 0                       # upper limit of cores available in that queue (compared with nodes × cores per node). 0 means no check
nodes_max = 0                       # upper limit of nodes. 0 means no check
walltime_max = ""                   # upper limit of the time limit (HH:MM:SS). Empty means no check. If exceeded, the validation before generation stops
modules = []                        # what to module load for every calculation code (for example ["intel-mpi"]). Loaded before code_modules below

[profiles.remote.code_modules]   # what to module load for each calculation code (type module avail on the cluster to see the list)
dftbplus = ["<DFTB+ module name>"]

# vasp = ["<MPI module name>", "<math library module name>"]

# [profiles.remote.env]          # environment variables exported at the top of submit.sh (only when using VASP; remove the # to use them)

# VASP_PP_PATH = "<parent directory of the POTCAR library>"

[profiles.remote.commands]       # run command for each calculation code. {mpiprocs} {omp_threads} {binary} {ntasks} are filled in
                                 # {ntasks} = nodes × MPI processes per node (used for mpirun -np across several nodes). Any other {...} stops the validation before generation
dftbplus = "dftb+"

# vasp = "mpirun -np {ntasks} <VASP bin>/vasp_{binary}"

# orca = "<directory where ORCA is unpacked>/orca"   # call ORCA by its absolute path, without mpirun
```

After generation there are three steps (`transfer_and_submit.sh` contains all three, with `host` and `remote_dir` filled in).

```bash
rsync -av <generated directory> <cluster>:<working directory>/       # 1. transfer (scp -r also works)
ssh <cluster> 'cd <working directory>/<name> && qsub submit.sh'     # 2. submit (sbatch on Slurm)
ssh <cluster> 'qstat -u $USER'                                      # 3. check status (squeue on Slurm)
```

On some systems the `module` command is not defined under `#!/bin/sh`, so `submit.sh` reads `/etc/profile` before using `module`. If it still cannot find the command, it shows the reason and stops.

Cluster profiles also generate `check_remote.sh`. Type `bash check_remote.sh` on your local PC, and it tries, in order, logging in with `ssh -o BatchMode=yes`, `module load`, whether the executable exists (`command -v`),
whether `remote_dir` is writable, and `sbatch --test-only submit.sh` on Slurm or `qstat -Q <queue>` on PBS. It prints one OK / NG line for each. **It does not submit anything** (it needs `host` and `remote_dir`).

## What to prepare for each code

| Calculation code | Program | Parameters |
|---|---|---|
| DFTB+ | `dftbplus` from conda-forge or similar. Put `dftb+` on PATH | Get the Slater-Koster parameters from https://dftb.org/parameters/download.html and place them under `sk_root` (CC BY-SA 4.0). The skf files for the element pairs needed, plus LICENSE / README, are copied to `skf/` in the output directory |
| VASP | Provided by a user with a license. Put the path of the executable in `commands.vasp` | Place the POTCAR library in the hierarchy `potpaw_PBE/<name>/POTCAR` and put its parent directory in `env.VASP_PP_PATH`. **POTCAR itself is not included in the output** (ADIT writes `potcar.spec` and `make_potcar.sh`, which concatenate it at run time) |
| xtb | `xtb` from conda-forge, on PATH (LGPL-3.0) | Not needed (built into the method) |
| Quantum ESPRESSO | `pw.x` on PATH, from `qe` on conda-forge or similar | Place UPF files in `pseudo_root/<set name>/` (pslibrary, SSSP https://www.materialscloud.org/discover/sssp and so on; check the license with the distributor). The UPF files for the elements needed, and the accompanying documents, are copied to `pseudo/` |
| ORCA | Register and download it, then put the path of the executable in `commands.orca` (no redistribution) | Not needed |
| CP2K | `cp2k.psmp` on PATH, from `cp2k` on conda-forge or similar (GPL-2.0-or-later) | The data directory that comes with CP2K (BASIS_MOLOPT, GTH_POTENTIALS and so on). Write its location in `cp2k_data` in the settings file (if empty, ADIT looks at the environment variable `CP2K_DATA_DIR`, then at `share/cp2k/data` next to the cp2k on PATH). Only the entries used are copied to `BASIS_adit` / `POTENTIAL_adit` |
| LAMMPS | `lmp` on PATH, from `lammps` on conda-forge or similar (GPL-2.0) | Force-field files (EAM, ReaxFF ffield, machine-learning potential models) or a data file made elsewhere. ADIT does not decide force-field coefficients |
| GROMACS | `gmx` on PATH, from `gromacs` on conda-forge or similar (LGPL-2.1) | A topology (.top / .itp) and structure (.gro / .pdb) made elsewhere, with CHARMM-GUI, acpype, pdb2gmx and so on |
| Machine-learning potential (MACE, CHGNet) | `pip install ase mace-torch` (or `chgnet`) in the environment that runs the generated `run_mlip.py`. This brings in PyTorch and takes several GB. **ADIT itself does not use it** | A trained model. If left blank, the package's default model is used (downloaded automatically on first use). If you specify a model file, it is copied to the output directory. Licenses differ by model, so check with the distributor |

### The CP2K, LAMMPS and GROMACS screens

Choose a code in the calculation code dropdown, and that code's fields appear. Blue labels with a "required" tag are required fields (hover over a label to see an explanation).

- **CP2K**: "Functional", "Cutoff [Ry]" and "Relative cutoff [Ry]" have no defaults, so decide them yourself (follow the procedure for converging CUTOFF and REL_CUTOFF in the CP2K manual). The basis set and pseudopotential for each element are offered in dropdowns, with candidates read from the files in the data directory. For an element with only one candidate, leaving the field blank uses it. If there are non-periodic directions, as with molecules and surfaces, choose a "Poisson solver". For a molecule with no box, enter "Box edge for a molecule [Å]". For settings not on the screen, write a header such as `[FORCE_EVAL/DFT/SCF]` in "Extra lines (per section)" and put the lines below it.
- **LAMMPS**: Write the units (metal / real), pair_style and pair_coeff, and put the force-field files in "Files to copy" (use only the file name in the input). If the data file is left empty, ADIT writes `data.lammps` from the structure (atomic or charge style). When you use an external data file, enter the elements of the atom types in order, and load the same system into the structure group too. The "k-points (periodic)" group does not appear.
- **GROMACS**: Specify the topology (.top) and the structure file (.gro / .pdb). **The structure file is the master copy of the structure.** Pressing "Also load into the structure group" sets the structure group's source to "File" and loads the same file (it is used to check the atom count and for the 3D view; in the web app, tick "use this file for the structure too"). For NPT, "Barostat (pcoupl)" and "Compressibility [1/bar]" are required (the water example in the GROMACS manual uses 4.5e-5 /bar). When you run stages such as NVT → NPT → production, set the previous stage's `adit.cpt` as "Previous stage .cpt". The "k-points (periodic)" group does not appear.

Calculation types that a code's generator does not have (vibrational analysis and band structures for LAMMPS and GROMACS, band structures for CP2K) show the reason in the "Cannot generate" pane when you select them. Other codes are handled the same way.

ORCA itself cannot be redistributed, so only the contents of its input files are tested. VASP has been checked by running it on a cluster, and DFTB+, xtb, pw.x, CP2K, LAMMPS, GROMACS, OpenMM, Psi4, ABINIT and PLUMED by running them on a local PC.

## Output directory

```
<output directory>/
  submit.sh      run script (local: bash submit.sh / cluster: qsub or sbatch submit.sh)
  spec.json      the full calculation settings. Opening it in ADIT restores them
  README.txt     how to run, and how to read the output files
  analyze.py     analysis script (writes figures and a summary to analysis/). Only for codes that support analysis
  transfer_and_submit.sh   transfer and submit commands. Only when the profile is pbs / slurm
  check_remote.sh          connectivity check before submitting (ssh, module, executable, working directory, queue). Same condition. It does not submit anything
  dftb_in.hsd, geometry.gen, skf/       for DFTB+
  INCAR, POSCAR, KPOINTS, potcar.spec, make_potcar.sh   for VASP
  struct.xyz, xtb.inp                   for xtb
  pw.in, pseudo/                        for Quantum ESPRESSO
  orca.inp                              for ORCA
```

