# examples/prep_stage1_qe_dftu_si —— generated files for pw.x DFT+U (HUBBARD card) and starting_magnetization and the results of actually running them

Si (diamond structure, pslibrary Si.pbe-n-rrkjus_psl.1.0.0.UPF), ecutwfc 20 Ry, k-points 2×2×2, nspin = 2, smearing (gaussian 0.01 Ry),
starting_magnetization(Si) = 0.2, HUBBARD {atomic} with U Si-3p 1.0 eV. The values only check the mechanism and carry no physical meaning. QE 7.5 (conda-forge), 2026-09-12.

| Item | Value (output.log) |
|---|---|
| U as read | "Hubbard projectors: atomic", "U(Si-3p) = 1.0000" (Dudarev) |
| Magnetization | Starts at 0.72 Bohr mag/cell and goes to 0.00 (Si is nonmagnetic) |
| Total energy | -22.57798452 Ry (Hubbard energy 0.07069816 Ry), converged in 6 iterations |

tmp/ holds the intermediate files of pw.x (wavefunctions and charge density).
