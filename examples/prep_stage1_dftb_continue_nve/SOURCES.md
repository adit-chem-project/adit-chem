# examples/prep_stage1_dftb_continue_nve —— generated files for continuing a previous run (DFTB+ MD → NVE) and the results of actually running them

Made with `adit.continuation.continue_from("examples/dftb_md_water_generated", settings)` (on the command line, `adit-gen --continue-from …`).
The settings are the MD of the original spec.json changed to NVE, 10 steps, output every step. Positions and velocities come from the last frame of the original geo_end.xyz
(velocities converted from Å/ps to Å/fs for structure.velocities in spec.json, and written in dftb_in.hsd as Velocities [AA/ps]). DFTB+ 25.1, 2026-09-12.

| Item | Value |
|---|---|
| Last temperature of the original MD (examples/dftb_md_water_generated/md.out) | 5885.4585 K |
| First temperature of the continuation (md.out) | 5885.4585 K |

The two match, so the velocities carry over, units included. The high temperature comes from the original example (MD of distorted water) and has nothing to do with this check of the continuation.
