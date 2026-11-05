# examples/cp2k_h2o_md_generated —— CP2K generator output (a few MD steps) and the results of actually running it

NVT MD of only 5 steps with the same electronic-structure settings as `examples/cp2k_h2o_generated`. CP2K 2026.2 (conda-forge).

| Item | Details |
|---|---|
| Settings | ENSEMBLE NVT, THERMOSTAT CSVR (TIMECON 100 fs), TEMPERATURE 300 K, TIMESTEP 0.5 fs, STEPS 5, trajectory every step (XYZ) |
| Run | `mpirun -np 1 cp2k.psmp -i cp2k.inp` (OpenMP 1 thread). 2026-09-12 |
| Result | Normal termination. Run time 10.2 s. Conserved quantity (Cons Qty) in `adit-1.ener`: -17.217974062 → -17.217957779 Hartree (a change of 1.6e-5 Hartree over 5 steps). Temperature goes from initial velocities at 300 K to 104.9 K |
| Output | `adit-pos-1.xyz` (trajectory, 6 frames), `adit-1.ener` (time, kinetic energy, temperature, potential energy, conserved quantity), `adit-1.restart` |

With a 3-atom molecule run for only 2.5 fs, the temperature is a value taken while the initial velocities are still moving into vibrations. The values only confirm that it runs.

## About the name (2026-09-14)

This calculation **was run before the rename (VISTA)**. For the rename (VISTA → ADIT),
file names and strings inside the files were changed to ADIT (to match what ADIT generates now).
**The calculation itself (numbers, convergence, software versions) is as it was run.**
