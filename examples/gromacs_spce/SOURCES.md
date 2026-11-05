# examples/gromacs_spce —— the SPC/E water box bundled with GROMACS

Used as input for `examples/gromacs_spce_em_generated` (energy minimization) and `examples/gromacs_spce_nvt_generated` (the NVT that follows it).
ADIT does not create topologies (it takes ones made elsewhere). Here the water files that come with GROMACS are used.

| File | Source | License |
|---|---|---|
| `conf.gro` | Unmodified copy of `share/gromacs/top/spc216.gro` from GROMACS 2026.3 (conda-forge `gromacs`) (216 water molecules, 648 atoms, cubic box 1.86206 nm) | LGPL-2.1, same as GROMACS (https://www.gromacs.org/about.html) |
| `topol.top` | 8 lines written by hand on 2026-09-12. It only uses `#include` for the bundled `oplsaa.ff/forcefield.itp` and `oplsaa.ff/spce.itp`, then lists `SOL 216` | It only reads the bundled force field, so the force-field files themselves are not here (grompp reads them from the GROMACS install location) |

`spce.itp` uses SETTLE (rigid water) through `#ifndef FLEXIBLE`. SPC/E reference: Berendsen, Grigera, Straatsma, J. Phys. Chem. 91, 6269 (1987).
