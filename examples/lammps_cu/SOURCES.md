# examples/lammps_cu —— the Cu EAM potential from the LAMMPS distribution

Used as input for `examples/lammps_cu_nvt_generated`. ADIT neither creates nor chooses potentials (it only copies the file the user provides).

| File | Source | License |
|---|---|---|
| `Cu_u3.eam` | `potentials/Cu_u3.eam` from the LAMMPS distribution (https://raw.githubusercontent.com/lammps/lammps/stable_22Jul2025_update5/potentials/Cu_u3.eam, retrieved 2026-09-12). sha256 `3436c491a4c75ea8b7141adbc6ee382a118f5fdb47f609c2a660fc1eb772599f` | GPL-2.0, same as LAMMPS (https://docs.lammps.org/Intro_opensource.html) |

First line of the file: `UNITS: metal`, `CITATION: Foiles et al, Phys Rev B, 33, 7983 (1986)`. Cite the paper by Foiles et al. when you use it.
The conda-forge `lammps` package (22 Jul 2025 update 5) did not include potentials/, so the file came from the tag for the same version on GitHub.
