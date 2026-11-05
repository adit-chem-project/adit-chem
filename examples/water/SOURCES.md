# Sources for examples/water

The input in this directory reproduces "First calculation with DFTB+" from the official DFTB+ recipes.
All numbers and keywords come from the official material; nothing was filled in from memory.

## URLs consulted and what was taken from them

| URL | Content taken |
|---|---|
| https://dftbplus-recipes.readthedocs.io/en/latest/basics/firstcalc.html | Full `dftb_in.hsd` for the water molecule, description of the gen format, description of each block, file inclusion with `<<<`, how to run (`dftb+ \| tee output`) |
| https://dftbplus-recipes.readthedocs.io/en/latest/_downloads/5919f4094cd60c5c70c13b47928442f5/recipes.tar.bz2 | Distributed archive. `recipes/basics/firstcalc/dftb_in.hsd` (same as the web page), `recipes/scripts/get_slakos`, `recipes/slakos/SLAKO_DOWNLOADS`, `recipes/slakos/mio-ext/` (symbolic links to mio-1-1 and tiorg-0-1) |
| https://dftbplus-recipes.readthedocs.io/en/latest/introduction.html | Recommended DFTB+ installation (conda-forge, `dftbplus=*=nompi_*`), how to get SK files (`get_slakos`) |
| https://dftbplus-recipes.readthedocs.io/en/latest/licence.html | License of the recipes (CC BY-SA 4.0) |
| https://dftb.org/parameters/download.html | List of parameter sets, elements in mio (H-C-N-O-S-P), the warning not to mix sets |
| https://dftb.org/parameters/introduction.html | SK files exist per element pair; documentation such as angular momenta is attached at the end of each skf |
| https://github.com/dftbparams/mio (releases v1.1.0) | The mio-1-1 distribution (`mio-1-1.tar.xz`), `LICENSE` (CC-BY-SA 4.0), `README` (copyright holders, required citations, elements included) |
| https://github.com/dftbplus/dftbplus/releases/tag/25.1 | DFTB+ 25.1 source and `manual.pdf` |
| DFTB+ 25.1 `manual.pdf` | 2.2.2 GenFormat, 2.3.1 GeometryOptimisation (Optimiser/Optimizer, MovedAtoms, MaxSteps, OutputPrefix, Convergence/GradElem default 1e-4 Hartree/Bohr), 2.4.11 SlaterKosterFiles (explicit list / Prefix / Type2FileNames / environment variable DFTBPLUS_PARAM_DIR), MaxAngularMomentum, 2.5 Options, 2.6 Analysis (PrintForces), 2.10 ParserOptions, Appendix D gen format |

## Changes from the original

1. `Geometry` was changed to read a file with `<<< "geometry.gen"` (a form shown on the same page)
2. The SK file path was changed from `../../slakos/mio-ext/` to `../../slakos/mio-1-1/`
   (`mio-ext` is a set of links in the distributed archive, and the 4 O/H files all point to the real files in `mio-1-1`)
3. `#` comments were added to each block

## Run environment (this machine)

- DFTB+ 25.1 (conda-forge, OpenMP build `nompi_h0f3606b_100`), parser version 14
- Executable: `/home/<user>/miniforge3/envs/dftbplus/bin/dftb+`
- The bundled tools (`gen2xyz`, `xyz2gen`, `gen2cif`) are in the same environment
- The startup header shows "DFTB+ development version (commit: a23bfb2)"
  (this is how the conda package displays itself; the version in `conda list` is 25.1)
- Location of the Slater-Koster files (planned): `/home/<user>/adit/slakos/mio-1-1/` (in `.gitignore`)

## Run results (2026-09-10)

- Run: `cd examples/water && /home/<user>/miniforge3/envs/dftbplus/bin/dftb+ | tee output.log`. Exit code 0, ran to completion with "Geometry converged"
- 10 geometry steps (step 0-9), wall time 1.35 s (OpenMP 8 threads)
- Final total energy: -4.0779379326 Hartree (-110.9663 eV)
- Optimized structure (computed from `geom.out.gen`): O-H 0.9672 Å (both bonds), H-O-H 107.19°
  (initial structure O-H 1.2701 Å, 76.13°; the tutorial starts from a deliberately distorted structure)
- Largest force component (`detailed.out`): 2.0e-5 Hartree/Bohr. This meets the convergence criterion 1e-4
- 2 warnings. Neither affects the results
  1. `CalculateForces` has been renamed to `PrintForces` (automatic conversion from parser version 12 to 14)
  2. The stack size is 8 MB, and `ulimit -s unlimited` is recommended for large systems
- Files generated: `output.log` (standard output), `dftb_pin.hsd` (input with all defaults expanded), `detailed.out`, `band.out`,
  `charges.bin`, `geom.out.gen`, `geom.out.xyz`. `output.log`, `dftb_pin.hsd`, `detailed.out` and `geom.out.*` are
  kept in the repository for reference. `charges.bin` is excluded by `.gitignore`
- Orbitals checked in the documentation at the end of the skf files: H-H.skf has `<Shells>1s </Shells>`, O-O.skf has `<Shells>2s 2p </Shells>`.
  This agrees with `MaxAngularMomentum { O = "p"; H = "s" }`
