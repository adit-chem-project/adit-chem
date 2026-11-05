# examples/prep_stage1_dftb_stages —— generated files for MD in stages (minimization → NVT → NVE) and the results of actually running them

Made with `adit.stages.write_stages` (on the command line, `adit-gen spec.json out/ --stages stages.json`). The structure is the distorted test water (water_spec in tests/conftest.py),
the SK set is mio-1-1, DFTB+ 25.1 (conda-forge, shown as a development version). `bash submit.sh` was run once (OMP 1 thread; on 2026-09-19 it was regenerated with the current generator and rerun).

| Stage | Settings | Result |
|---|---|---|
| stage_01_min | Geometry optimization (at most 100 steps) | Normal termination. geom.out.gen |
| stage_02_nvt | NVT, Berendsen 300 K (time constant 20 fs), 0.5 fs × 40 steps | The starting structure matches geom.out.gen of stage_01 to 10 digits (handoff.py rewrites geometry.gen) |
| stage_03_nve | NVE, 0.5 fs × 20 steps, velocities from stage_02 | First temperature 314.6884 K = last temperature of stage_02, 314.6884 K. Total energy change 1.5e-6 Hartree (20 steps) |

In the first generation stage_03 stopped: when velocities (Velocities) are given and `Thermostat = None { InitialTemperature … }` is written, DFTB+ 25.1
stops with it as an "item not read". The generator was fixed (`Thermostat = None {}` when velocities are present), and generating again gave the results above.
