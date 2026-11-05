# examples/vasp_cd_si_generated —— VASP input made by the generator (equivalent to the wiki example Cd Si relaxation)

The same settings as `examples/vasp_cd_si` (the wiki original) copied into a Spec and generated with the ADIT generator. POTCAR is not included.
A person runs it on a cluster (not done).

## Differences from the original, and why

| Item | Wiki original | Generated files | Reason |
|---|---|---|---|
| POSCAR | Old format without an element-name line, scale 5.5 | Has an element-name line, scale 1 with lattice vectors multiplied by 5.5 | ASE cannot read the old format. The numbers are the same |
| ISTART = 0, ICHARG = 2 | Explicit in INCAR | Same values written through extra_incar | The way out for keys the GUI does not have |
| ALGO / EDIFF / NELM / LREAL / ISPIN / PREC | Not written | Explicit | Policy: always write the Spec values |
| KPOINTS | Monkhorst Pack 11 11 11 | Same | |
| POTCAR | (assumed) | None. potcar.spec (Si → Si) and make_potcar.sh | Cannot be distributed |

## Run on a cluster (2026-09-10, VASP 6.4.2, submitted by a person)

The generated files were transferred as is, and `make_potcar.sh` assembled POTCAR (`PAW_PBE Si 05Jan2001`) and ran.
It finished after 10 steps, as NSW = 10 says (it stops the same way as the wiki example; EDIFFG = -0.0001 is not met in 10 steps).
Final energy F = -10.8244 eV. The z of the shifted atom went from 0.130 to 0.1275 (the original's 4.81250 Å ≈ 0.1275). No warnings or errors.
CPU time 138 s (32 cores, 11×11×11 k-points).
