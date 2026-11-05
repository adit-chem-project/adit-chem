# Sources for examples/vasp_h2o

The VASP wiki example "H2O" (structural relaxation of a water molecule), reproduced **verbatim** as the 3 files POSCAR / INCAR / KPOINTS.
POTCAR is not included (it cannot be distributed to anyone without a license; docs/design_vasp_v1.md section 2).

| Item | Details |
|---|---|
| URL | https://www.vasp.at/wiki/index.php/H2O |
| License | GNU Free Documentation License 1.2 (the wiki default; each page footer says "Content is available under GNU Free Documentation License 1.2 unless otherwise noted"). Reproduction requires stating the source and the same license; this document does that |
| Retrieved | 2026-09-10 |
| System | 1 water molecule. Cubic box of 15 Å. Coordinates scaled by 0.52918 (the factor that converts coordinates written in bohr to Å) |
| Method | PAW-PBE (from POTCAR), PREC = Normal, ENMAX = 400, ISMEAR = 0 / SIGMA = 0.1 |
| Geometry optimization | IBRION = 1 (DIIS), NFREE = 2, NSW = 10, EDIFFG = -0.02 (forces on all atoms below 0.02 eV/Å). Atom 1 (at the origin) is fixed; the other 2 move only in x and y (Selective dynamics) |
| k-points | Γ point only (1×1×1) |

## About POTCAR (reference for the `potcar.spec` the generator writes)

This POSCAR is in the **old format without an element-name line**, where `1 2` on line 6 means "1 of the first species, 2 of the second".
From the coordinates, the first species reads as O fixed at the origin and the second as H placed at ±1.43.
So POTCAR is concatenated in the order **O, H**. The generator uses the new format with an element-name line (VASP 5 and later).
In the Materials Project table (POTCAR section of pymatgen v2023.10.11 `MPRelaxSet.yaml`, MIT), H → `H` and O → `O`.

## Notes

- Source for axis constraints (2026-09-12): https://vasp.at/wiki/POSCAR . The T/F of Selective dynamics refers to lattice-vector directions even when coordinates are written as Cartesian. ADIT converts Cartesian axis constraints to lattice directions only when the allowed displacements stay the same. An ISIF that changes the cell shape (https://vasp.at/wiki/ISIF) cannot preserve the directions, so the check stops on partial Cartesian axis constraints.
  Axis-constraint source (2026-09-12): https://vasp.at/wiki/POSCAR . Selective dynamics T/F flags refer to lattice-vector directions even with Cartesian positions. ADIT converts Cartesian axis constraints only when the allowed displacements remain equivalent. Partial Cartesian constraints are rejected when ISIF allows the cell shape to change (https://vasp.at/wiki/ISIF), because the fixed directions cannot be retained.

- The original is written for VASP 5 (`PREC=Normal (Default for VASP.5.X)`, the old key `ENMAX`). A person checks on a cluster whether it works on 6.x
- Wiki descriptions of each keyword: https://www.vasp.at/wiki/index.php/ENMAX , /PREC , /ISMEAR , /SIGMA , /IBRION , /NFREE , /NSW , /EDIFFG , /POSCAR , /KPOINTS (existence confirmed on 2026-09-10)
