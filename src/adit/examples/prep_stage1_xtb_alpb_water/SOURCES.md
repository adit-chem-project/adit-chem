# examples/prep_stage1_xtb_alpb_water —— generated files for xtb solvent (ALPB, water) and the results of actually running them

Single point with XtbMethod(solvation="alpb", solvent="water"). `--alpb water` appears on the xtb line of submit.sh. xtb 6.7.1 (conda-forge), 2026-09-12.

| Item | Value |
|---|---|
| Solvent model (output.log) | ALPB, Solvent water |
| Solvation free energy Gsolv | -0.029383426097 Eh |
| Total energy | -4.997357957784 Eh (the same structure in the gas phase is -4.977864911515 Eh) |

The list of solvent names (SOLVENTS in codes/xtb.py) was made by running each of the 35 names in xtb --help and the official documentation tables once per method (GFN1 / GFN2 / GFN-FF) × model (ALPB / GBSA),
and keeping only those xtb did not stop on.
