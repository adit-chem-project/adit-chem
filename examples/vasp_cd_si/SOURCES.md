# Sources for examples/vasp_cd_si

The VASP wiki example "Cd Si relaxation" (relaxation of internal coordinates of diamond Si), reproduced **verbatim** as the 3 files
POSCAR / INCAR / KPOINTS. POTCAR is not included.

| Item | Details |
|---|---|
| URL | https://www.vasp.at/wiki/index.php/Cd_Si_relaxation |
| License | GNU Free Documentation License 1.2 (the wiki default). Reproduction requires stating the source and the same license; this document does that |
| Retrieved | 2026-09-10 |
| System | Diamond Si, 2 atoms. fcc lattice with lattice constant 5.5 Å (fractional coordinates). The z of one atom is shifted from 0.125 to 0.130 to break the symmetry |
| Method | ENCUT = 240, ISMEAR = 0 / SIGMA = 0.1, ISTART = 0 / ICHARG = 2 |
| Geometry optimization | IBRION = 2 (conjugate gradient), NSW = 10, ISIF = 2 (internal coordinates only; the lattice is fixed), EDIFFG = -0.0001 |
| k-points | 11×11×11 Monkhorst-Pack (odd, so Γ-centered) |
| Result given in the original | Forces at or below 1e-3 eV/Å after 10 steps, total drift 0 |

## About POTCAR

Old format without an element-name line (line 6 is `2`). In the Materials Project table (pymatgen v2023.10.11 `MPRelaxSet.yaml`, MIT), Si → `Si`.
The identifying line in the PAW PBE distribution is `PAW_PBE Si 05Jan2001`.

## Notes

- The original is written for VASP 5. A person checks on a cluster whether it works on 6.x
- Wiki descriptions of each keyword: https://www.vasp.at/wiki/index.php/ENCUT , /ISMEAR , /SIGMA , /ISTART , /ICHARG , /NSW , /IBRION , /ISIF , /EDIFFG , /POSCAR , /KPOINTS , /INCAR (existence confirmed on 2026-09-10)
