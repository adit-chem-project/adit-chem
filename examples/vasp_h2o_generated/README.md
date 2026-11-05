# examples/vasp_h2o_generated —— VASP input made by the generator (equivalent to the wiki example H2O)

The same settings as `examples/vasp_h2o` (the wiki original) copied into a Spec and generated with the ADIT generator (`codes/vasp.py`).
**POTCAR is not included.** `potcar.spec` and `make_potcar.sh` are included instead.
This machine has no POTCAR library, so TITEL / ZVAL are "unverified". A person runs it on a cluster (not done).

## Differences from the original, and why

| Item | Wiki original | Generated files | Reason |
|---|---|---|---|
| POSCAR format | Old format without an element-name line, coordinates scaled by 0.52918 | Has an element-name line (VASP 5 format), scale 1 with cell and coordinates multiplied by 0.52918 | ASE cannot read the old format. The numbers are the same |
| Coordinate style | cart (Cartesian) | Direct (fractional) | ASE writer default. Equivalent |
| Selective dynamics | O is F F F, H is T T F (z fixed) | O is F F F, H is T T T | v1 fixes whole atoms (fixed_atoms). Per-axis constraints are not supported |
| ENMAX = 400 | Old key | ENCUT = 400 | The current key (ENCUT page on the wiki). Same value |
| IBRION = 1, NFREE = 2 | DIIS | IBRION = 2 (conjugate gradient) | Generator default. IBRION can be overridden with extra_incar |
| ALGO / EDIFF / NELM / LREAL / ISPIN | Not written (VASP defaults) | Explicit (Normal / 1e-4 / 60 / Auto / 1) | Policy: always write the Spec values |
| ISIF | Not written (default 2) | ISIF = 2 | Same as above |
| KPOINTS | Monkhorst Pack 1 1 1 | Gamma 1 1 1 | Γ point only is given as Gamma (the same point) |
| POTCAR | (the wiki assumes you have it) | None. potcar.spec (O → O, H → H) and make_potcar.sh | Cannot be distributed |

## Run on a cluster (2026-09-10, VASP 6.4.2, submitted by a person)

The generated files were transferred as is, and `make_potcar.sh` assembled POTCAR (`PAW_PBE O 08Apr2002`, `PAW_PBE H 15Jun2001`) and ran.
"reached required accuracy" (converged) in 5 steps. Final energy F = -14.2226 eV. No warnings or errors.
The wiki's VASP 5 style (ENMAX) changed to ENCUT worked on 6.4.2. CPU time 5.4 s (32 cores).
