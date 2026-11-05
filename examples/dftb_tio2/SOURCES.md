# Sources for examples/dftb_tio2 (reference example of periodic DFTB+)

Reproduces stage 1 of the official recipe "Band structure, DOS and PDOS" (self-consistent charges of anatase TiO2 at converged k-points).

| Item | Details |
|---|---|
| URL | https://dftbplus-recipes.readthedocs.io/en/latest/basics/bandstruct.html |
| Original | `recipes/basics/bandstruct/1_density/dftb_in.hsd` in the distributed archive recipes.tar.bz2 (the recipes are CC BY-SA 4.0) |
| SK sets | mio-1-1 (O-O) and tiorg-0-1 (Ti-Ti, Ti-O, O-Ti). Both CC BY-SA 4.0 (dftbparams, 2024-12-18 release). tiorg is an extension of mio, marked "Requires: mio" on dftb.org, so mixing them is intended |
| Location | `slakos/mio-1-1/`, `slakos/tiorg-0-1/`, and `slakos/mio-ext/`, a set of links with the same layout as the recipes (all in .gitignore) |
| Basis for the angular momenta | `<Shells>3d 4s 4p</Shells>` at the end of tiorg-0-1/Ti-Ti.skf → Ti = "d". `<Shells>2s 2p</Shells>` in mio-1-1/O-O.skf → O = "p" |

## Changes from the original

1. The structure is read with `<<< "geometry.gen"` (contents identical to the original)
2. The Prefix of the SK files is `../../slakos/mio-ext/`
3. `#` comments

## Run results (2026-09-10, DFTB+ 25.1, this machine)

- Exit code 0. Periodic boundaries Yes, 32 k-points (4×4×4, shift 0.5)
- SCC converged in 5 iterations. Total energy -15.6210745499 Hartree (-425.0711 eV)
- 3 warnings: parser version conversion 12→14, the dipole is not defined for a periodic system, and stack size. None affects the results
- Files generated: output.log, dftb_pin.hsd, detailed.out, band.out, dos_ti.{1,2,3}.out, dos_o.{1,2}.out (`charges.bin` is in .gitignore)
- The recipe goes on to produce the density of states with `dp_dos`; this reference stops at the SCC calculation
