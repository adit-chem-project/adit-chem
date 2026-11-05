# How dftb_in.hsd maps to VASP input

VASP splits its roles across 4 files, INCAR / POSCAR / POTCAR / KPOINTS,
while DFTB+ splits the same roles across blocks in one file, `dftb_in.hsd`.
The table below maps "roles". Keywords do not translate one to one.

## Mapping by file and block

| VASP | DFTB+ | Notes |
|---|---|---|
| POSCAR | `Geometry = GenFormat { <<< "geometry.gen" }` | gen format. Line 1 is "number of atoms C/S/F", line 2 the element list, then "index element-number x y z" (Å). `C` is a cluster (non-periodic), so there are no lattice-vector lines. For a periodic system, use `S` and add 1 origin line + 3 lattice-vector lines at the end |
| POTCAR | `Hamiltonian = DFTB { SlaterKosterFiles {...} }` | The biggest difference: **per element pair, not per element** (O-O, O-H, H-O, H-H). They can be combined only within the same set. The files themselves go in `~/slakos/` or similar, and the input refers to them by path (nothing corresponds to concatenating POTCARs) |
| POTCAR valence configuration (which orbitals it has) | `MaxAngularMomentum { O = "p"; H = "s" }` | The skf format does not record "up to which orbital", so this is set by hand. The value comes from the set's documentation (end of the skf). The GUI fills it in automatically (design document, section 2.5) |
| Method settings in INCAR | `Hamiltonian = DFTB { Scc = Yes ... }` | See the table below |
| Geometry optimization settings in INCAR (IBRION, NSW, EDIFFG, ISIF) | `Driver = GeometryOptimization {...}` | See the table below |
| KPOINTS | `Hamiltonian = DFTB { KPointsAndWeights {...} }` | Not needed for molecules (clusters), so it is not in this input. Required for periodic systems (v1 and later) |
| Output control in INCAR (LWAVE, LCHARG, NWRITE, etc.) | `Options {}` and `Analysis {}` | `Options` controls whether files are written (detailed.out, results.tag, charges.bin, etc.); `Analysis` sets additional quantities to compute and print (all force components, DOS, electrostatic potential, etc.) |
| (no equivalent) | `ParserOptions { ParserVersion = 12 }` | Declares the input syntax version. VASP has no equivalent. Newer DFTB+ converts older input automatically |

## Mapping of INCAR keywords (electronic structure)

| INCAR | DFTB+ (inside `Hamiltonian = DFTB {}`) | Notes |
|---|---|---|
| (SCF is always self-consistent) | `Scc = Yes` | In DFTB, self-consistent charges (SCC) can be turned off (`No` = non-SCC DFTB). Normally `Yes` |
| EDIFF (energy criterion) | `SccTolerance` (default 1e-5) | The units differ. DFTB+ judges by the **charge** difference (number of electrons), not the energy difference |
| NELM | `MaxSccIterations` (default 100) | Same meaning. When exceeded, it warns, computes the energy with the charges so far, and stops |
| NELECT / total charge | `Charge` (default 0.0) | The sign is "net charge of the system". The opposite way of expressing it from VASP NELECT (number of electrons) |
| ISPIN=2 / MAGMOM | `SpinPolarisation = Colinear { UnpairedElectrons = ... }` | Not used here. The default is no spin polarization |
| ISMEAR / SIGMA | `Filling = Fermi { Temperature [K] = ... }` | The default is a Fermi distribution at nearly zero temperature (0.1e-7 Hartree). Normally left alone for molecules |
| AMIX / BMIX / IMIX | `Mixer = Broyden { MixingParameter = 0.2 }` | The default is Broyden mixing with a mixing parameter of 0.2 |
| IVDW | `Dispersion = ...` | Not used here. `DftD3`, `LennardJones`, etc. |
| (no equivalent; DFTB3) | `ThirdOrderFull = Yes` + `HubbardDerivs {...}` | Not used here. DFTB3 is used with a matching set such as 3ob |

## Mapping of INCAR keywords (geometry optimization)

| INCAR | DFTB+ (inside `Driver = GeometryOptimization {}`) | Notes |
|---|---|---|
| IBRION=1/2/3 (algorithm) | `Optimizer = Rational {}` | The options are Rational / LBFGS / FIRE / SteepestDescent (25.1 manual 2.3.1). Both spellings `Optimiser` and `Optimizer` work |
| NSW | `MaxSteps = 100` | Same meaning. The default is 200 |
| EDIFFG (negative: force criterion) | `Convergence { GradElem = 1E-4 }` | Converged when the **largest component** of the force is below 1e-4 Hartree/Bohr (about 5.14e-3 eV/Å). Can be written with units as `GradElem [eV/AA] = 5.14e-3`. Use `Energy` for an energy criterion, or `GradNorm` for the force norm |
| EDIFFG (positive: energy criterion) | `Convergence { Energy = ... }` | Disabled by default (inf) |
| ISIF (whether the lattice moves) | `LatticeOpt = Yes/No` | Irrelevant for molecules. Used for periodic systems |
| Selective dynamics (T/F columns in POSCAR) | `MovedAtoms = 1:-1` | The atoms to move are given in the input. `1:-1` means all atoms. Selecting by element, as in `O H`, also works |
| CONTCAR | `OutputPrefix = "geom.out"` → `geom.out.gen`, `geom.out.xyz` | Where the final structure goes. You choose the name |
| ISTART / ICHARG (restart) | `charges.bin` (written automatically) + `ReadInitialCharges = Yes` | SCC charges are saved in `charges.bin` (close to WAVECAR/CHGCAR) |
| (single point) | `Driver {}` or omit the `Driver` block | Equivalent to IBRION=-1, NSW=0 |

## Mapping of output files

| VASP | DFTB+ | Notes |
|---|---|---|
| OUTCAR / standard output | Standard output (`dftb+ \| tee output`) | SCC iterations for each geometry step, total energy, largest force component |
| OSZICAR | `iSCC ...` lines in standard output | Electronic energy and charge error for each SCC iteration |
| Final results in OUTCAR | `detailed.out` | Energy breakdown of the final step, Mulliken charges, forces (when `CalculateForces`/`PrintForces` is on) |
| EIGENVAL | `band.out` | Eigenvalues and occupations |
| CONTCAR | `geom.out.gen` / `geom.out.xyz` | See above |
| CHGCAR (for restart) | `charges.bin` | Charge restart data |
| (no equivalent) | `dftb_pin.hsd` | The input with all defaults expanded. **The most reliable way to check "what you left out"** |
| Close to vasprun.xml | `results.tag` (with `Options { WriteResultsTag = Yes }`) | Machine-readable results. ASE and others read it |

## About units

The default units for numbers are atomic units (energy Hartree, length Bohr, force Hartree/Bohr),
unlike the eV / Å of VASP. Only the coordinates in gen format are in Å.
Square brackets add a unit so you can write values without converting (e.g. `GradElem [eV/AA] = 5.14e-3`).

Source: the mapping above is based on the DFTB+ 25.1 manual (2.2-2.6, 2.10, Appendix D) and
the official recipe "First calculation with DFTB+". VASP keywords are described
as known to the user.
