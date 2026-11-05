# examples/plumed_lammps_cu_generated —— LAMMPS MD connected to PLUMED (run on 2026-09-13)

The same Cu NVT MD as `examples/lammps_cu_nvt_generated`, with a **PLUMED input written by the user** added.
**ADIT creates neither collective variables nor biases.** The contents of `plumed.dat` were written by a person for this example.

| Item | Value |
|---|---|
| Structure and force field | Face-centered cubic Cu (EAM from `examples/lammps_cu`; source in `SOURCES.md` in that directory) |
| Calculation | NVT MD, 300 steps, time step 1 fs, 300 K, Nosé–Hoover (time constant 100 fs), recorded every 50 steps |
| PLUMED | `d1: DISTANCE ATOMS=1,2` and `cn: COORDINATION GROUPA=1-32 GROUPB=1-32 R_0=3.0`, `PRINT ... FILE=COLVAR STRIDE=50` |
| Connection | `fix adit_plumed all plumed plumedfile plumed.dat outfile plumed.log` just before `run` in `in.lammps` |
| Run | LAMMPS 22 Jul 2025 (conda-forge, with the PLUMED package), PLUMED 2.10, 1 process, 1 thread |

## Watch the units (checked by actually running it)

In LAMMPS `metal` units distances are in Å, but **the default PLUMED unit is nm**.
`d1` in `COLVAR` is around 0.2556, which is 2.556 Å (the Cu nearest-neighbor distance) written in nm.
Writing `UNITS` in the PLUMED input changes this. ADIT does not convert this difference, and README.txt also carries a warning about it.

## Files

| File | Contents |
|---|---|
| `plumed.dat` | PLUMED input written by the user (copied unchanged by ADIT) |
| `COLVAR` | Record of collective variables written by PLUMED (starts with `#! FIELDS`). ADIT's analysis puts these column names and values in the summary |
| `plumed.log` | PLUMED's own log (`outfile` of `fix plumed`) |
| `in.lammps` / `data.lammps` / `Cu_u3.eam` | LAMMPS input, structure and force field |
| `log.lammps` / `output.log` / `traj.lammpstrj` / `final.data` | LAMMPS output |

## Sources

- PLUMED: https://www.plumed.org/doc-v2.10/user-doc/html/index.html ,
  How LAMMPS connects to it: https://docs.lammps.org/fix_plumed.html
- PLUMED references are listed at the end of `plumed.log` for each feature used (cite the ones listed there)
- LAMMPS: https://docs.lammps.org/Intro_citing.html
