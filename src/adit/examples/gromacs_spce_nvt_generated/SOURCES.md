# examples/gromacs_spce_nvt_generated —— GROMACS generator output (NVT) and the results of actually running it

NVT generated with the last structure of `examples/gromacs_spce_em_generated` (`adit.gro`) as the structure file. GROMACS 2026.3 (conda-forge).

| Item | Details |
|---|---|
| Input | Topology from `examples/gromacs_spce/topol.top`, structure from `examples/gromacs_spce_em_generated/adit.gro` |
| mdp | integrator md, dt 0.002 ps, nsteps 250 (0.5 ps), tcoupl V-rescale (thermostat csvr on screen), tau-t 0.1 ps, ref-t 300 K, gen-vel yes (gen-seed -1), PME, 0.85 nm |
| Run | grompp → `gmx mdrun -deffnm adit -ntmpi 1 -ntomp 1`. 2026-09-12 |
| Result | Normal termination. Run time 0.41 s. Averages over 250 steps: temperature 307.94 K, potential energy -9832.91 kJ/mol, pressure 140.9 bar (AVERAGES in adit.log) |
| Output | `adit.xtc` (trajectory, every 25 steps), `adit.edr` (energies), `adit.gro` (last structure), `adit.cpt` (continuation data for the next stage) |

gen-seed is -1 (the GROMACS default; different random numbers each time), so a rerun will not reproduce the values. 0.5 ps is too short for equilibration, and the values only confirm that it runs.

## About the name (2026-09-14)

This calculation **was run before the rename (VISTA)**. For the rename (VISTA → ADIT),
file names and strings inside the files were changed to ADIT (to match what ADIT generates now).
**The calculation itself (numbers, convergence, software versions) is as it was run.**
