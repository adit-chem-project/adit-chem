# examples/prep_stage1_cp2k_dftu_h2o —— generated files for CP2K DFT+U (KIND/DFT_PLUS_U) and the results of actually running them

Distorted water from tests/conftest.py, PBE, DZVP-MOLOPT-SR-GTH / GTH-PBE, CUTOFF 280 Ry, REL_CUTOFF 40 Ry, box 8 Å (MT).
U = 2.0 eV on the O 2p (L 1, U_MINUS_J [eV] 2.0 in &DFT_PLUS_U of KIND O; PLUS_U_METHOD MULLIKEN in the DFT section). The values only check the mechanism. CP2K 2026.2, 2026-09-12.

| Item | Value (output.log) |
|---|---|
| U as read | "DFT+U| Method MULLIKEN", "U(eff) = (U - J) value in [eV]: 2.000" |
| DFT+U energy | 0.04085351432522 Hartree |
| Total energy | -17.084127263102282 Hartree |

## About the name (2026-09-14)

This calculation **was run before the rename (VISTA)**. For the rename (VISTA → ADIT),
file names and strings inside the files were changed to ADIT (to match what ADIT generates now).
**The calculation itself (numbers, convergence, software versions) is as it was run.**
