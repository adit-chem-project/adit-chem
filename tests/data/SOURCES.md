# Sources of the real files used in tests

All of these are here to check formats **against real files**.

| File | Source | License | When and what was checked |
|---|---|---|---|
| `si.proj.projwfc_up` | projwfc.x from Quantum ESPRESSO 7.x run locally on 2 Si atoms | — (made here) | 2026-09-13. Order of filproj |
| `epsr_si.dat` / `epsi_si.dat` / `eels_si.dat` | QE epsilon.x run locally on Si (norm-conserving pseudopotential) | — (made here) | 2026-09-13. Column order, and that -Im(1/eps) matches QE's own value |
| `bader_ACF.dat` | bader 1.0.5 run locally on a cube with Gaussians holding 6 and 2 electrons | — (made here) | 2026-09-13. ACF.dat columns and the partitioned electron counts (5.9999 / 1.9987) |
| `PROCAR.simple` / `PROCAR.new_format_5.4.4.gz` | `test-files/io/vasp/outputs/` in [pymatgen](https://github.com/materialsproject/pymatgen) | MIT | 2026-09-13. The plain format and lm-decomposed + phase (VASP 5.4.4). The values read match pymatgen's `Procar` |
| `WAVEDER.gz` | Gzip of `test-files/io/vasp/outputs/WAVEDER` in [pymatgen](https://github.com/materialsproject/pymatgen) | MIT | 2026-09-13. Order of Fortran unformatted records. The arrays read match pymatgen's `Waveder` exactly |
| `gaussian_dvb_raman.out.gz` / `gamess_dvb_ir.out.gz` / `qchem_dvb_raman.out.gz` / `orca_dvb_raman.out.gz` | `data/` in [cclib](https://github.com/cclib/cclib) (real output from Gaussian 16, GAMESS-US 2018, Q-Chem 5.4, ORCA 6.0) | BSD-3 | 2026-09-13. Energies, structures, frequencies, IR intensities, **Raman activities**, Mulliken charges. The values read match cclib's own parser (only the GAMESS IR needs a unit conversion) |
| `water_hessian.out` | A DFTB+ run (2026-09-10) | — (made here) | Order of hessian.out |
| `fetch/pubchem_water_cids.json` / `fetch/pubchem_962_3d.sdf` / `fetch/pubchem_962_props.json` | Real responses from PubChem PUG REST (`compound/name/water/cids/JSON`, `compound/cid/962/SDF?record_type=3d`, `compound/cid/962/property/MolecularFormula,IUPACName,Title/JSON`) | NCBI policy (US government works are in the public domain) | 2026-09-19. Shape of the responses for name → CID, 3D SDF, and compound names |
| `fetch/cod_1000041.cif` / `fetch/cod_1000041_meta.json` | Real responses from COD (`https://www.crystallography.net/cod/1000041.cif`, `result?id=1000041&format=json`) | CC0 1.0 | 2026-09-19. CIF and the JSON field names (formula, title, authors, doi …) |
| `fetch/providers.json` / `fetch/idx_oqmd_links.json` / `fetch/oqmd_si.json` / `fetch/oqmd_4061352.json` | Real responses from OPTIMADE (`https://providers.optimade.org/providers.json` thinned to 3 providers, `index-metadbs/oqmd/v1/links`, `oqmd.org/optimade/v1/structures?filter=chemical_formula_reduced="Si"` thinned to 2 entries, and the first of those reshaped as a single-entry response) | providers.json is from the OPTIMADE consortium (MIT); OQMD entries follow the OQMD terms | 2026-09-19. The 3-step shape providers → links → structures, and lattice_vectors / species / cartesian_site_positions |
| `fetch/mp_summary_mp-149.json` | **Not a real response.** Written by hand to fit the Materials Project OpenAPI (`TypedStructureDict` / `TypedLatticeDict` / `TypedSiteDict` in `https://api.materialsproject.org/openapi.json`), because there was no API key to get a real one | — (made here) | 2026-09-19. Only the shape of `data[0].structure.lattice.matrix` and `sites[].species[].element` / `xyz` |

**The development machine has no VASP, so VASP output is checked with real files that are publicly available.**

## Checked against real files, but not in the repository

Their licenses do not fit ADIT (MIT), so they are **not included**. The formats were checked against real files on 2026-09-13.

| What | Source | License | How to run the test |
|---|---|---|---|
| OpenMX `.out` | `work/input_example/Methane.out` and `H2O.out` bundled with OpenMX (published at [FermiQ/openmx-square](https://github.com/FermiQ/openmx-square)) | GPL-3 | `ADIT_OPENMX_OUT=<file> pytest tests/test_readers_openmx.py` |
| VASP `vasprun.xml` + `WAVEDER` (check against the dielectric function) | `fixtures/reproduce_eps` in pymatgen | MIT, but large at 2.7 MB | Done locally only (results recorded in the notes of `analysis/waveder.py`) |
