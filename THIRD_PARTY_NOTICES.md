# Bundled third-party files and their licenses

**ADIT's own code (`src/adit/`, `tests/*.py`, documentation) is MIT** (`LICENSE`).
However, `examples/` and `tests/data/` contain **files taken from other projects**.
Those files **keep their original licenses** and do not become MIT. This document lists them.
The `SOURCES.md` in each directory gives the retrieval date, URL and sha256.

The full license texts are in `licenses/`.

## Not included in the PyPI distribution (adit-chem)

What `pip install adit-chem` installs (wheel / sdist) is **only ADIT's own code**.
The third-party files below are only in the GitHub repository.

## Files taken from calculation code distributions

| File | Source | License | Full text |
|---|---|---|---|
| `examples/lammps_cu/Cu_u3.eam` and its copies `examples/lammps_cu_nvt_generated/Cu_u3.eam`, `examples/plumed_lammps_cu_generated/Cu_u3.eam` | `potentials/` in LAMMPS | GPL-2.0 | `licenses/GPL-2.0.txt` |
| `examples/cp2k_h2o_generated/BASIS_adit`, `POTENTIAL_adit` and other CP2K examples | `data/` in CP2K (parts of BASIS_MOLOPT and GTH_POTENTIALS) | GPL-2.0-or-later | `licenses/GPL-2.0.txt` |
| `examples/gromacs_spce/conf.gro` and its copies, the `conf.gro` in `examples/gromacs_spce_em_generated/`, `examples/gromacs_spce_nvt_generated/`, `examples/openmm_spce_nvt_generated/` | `share/gromacs/top/spc216.gro` in GROMACS | LGPL-2.1 | `licenses/LGPL-2.1.txt` |
| `.skf`, `LICENSE` and `README` in `examples/water_generated/skf/`, `examples/dftb_md_water_generated/skf/`, `examples/dftb_tio2_generated/skf/`, `examples/dftb_tio2_bands_generated/skf/`, `examples/prep_stage1_dftb_continue_nve/skf/`, `examples/prep_stage1_dftb_stages/*/skf/`, `examples/prep_stage2_dftb_adsorption/*/skf/`, `examples/prep_stage3_dftb_elastic/*/skf/`, `examples/prep_stage3_dftb_neb_images/*/skf/`, `examples/prep_stage3_dftb_phonons_ase/*/skf/`, `examples/prep_stage3_dftb_phonons_phonopy/*/skf/` (92 `.skf` files in total) | Slater-Koster set mio-1-1 from dftb.org. The pairs containing Ti in `dftb_tio2*` (`Ti-Ti`, `Ti-O`, `O-Ti`) are from tiorg-0-1, also from dftb.org (https://dftb.org/parameters/download.html) | CC BY-SA 4.0 | `licenses/CC-BY-SA-4.0.txt` (each `skf/LICENSE` also has the same full text) |

`topol.top` in the same directories (`examples/gromacs_spce/` and its copies) is a short file that only `#include`s force fields shipped with GROMACS. It was written on the ADIT side.

**ADIT does not create these files.** It only copies what the user has provided. The examples are here to show a setup that works.

## Taken from documentation and recipes

| File | Source | License | Full text |
|---|---|---|---|
| `examples/vasp_h2o`, `examples/vasp_cd_si` | Examples from the official VASP wiki (reproduced verbatim) | GNU Free Documentation License 1.2 | `licenses/GNU-FDL-1.2.txt` |
| `examples/dftb_tio2`, `examples/water`, `examples/prep_stage2_dftb_adsorption` | The official DFTB+ recipes | CC BY-SA 4.0 | `licenses/CC-BY-SA-4.0.txt` |

## Test files from other projects, kept to check file formats

| File | Source | License |
|---|---|---|
| `tests/data/gaussian_dvb_raman.out.gz`, `gamess_dvb_ir.out.gz`, `qchem_dvb_raman.out.gz`, `orca_dvb_raman.out.gz` | `data/` in [cclib](https://github.com/cclib/cclib) | BSD-3-Clause |
| `tests/data/PROCAR.simple`, `PROCAR.new_format_5.4.4.gz`, `WAVEDER.gz` | `test-files/` in [pymatgen](https://github.com/materialsproject/pymatgen) | MIT |
| `tests/data/fetch/pubchem_*` | Real responses from PubChem PUG REST | NCBI policy (works of the US government are in the public domain) |
| `tests/data/fetch/cod_1000041*` | Real responses from the Crystallography Open Database | CC0 1.0 |
| `tests/data/fetch/providers.json`, `idx_oqmd_links.json`, `oqmd_*.json` | Real OPTIMADE responses (providers.json from the OPTIMADE consortium, the rest from OQMD) | providers.json is MIT; the OQMD entries are under the OQMD terms |

**These are here to check, against real file formats, that the parsers really read them correctly.**
The sources and what was checked are in `tests/data/SOURCES.md`.

## Structures fetched from databases (fetched by the user)

Files that arrive through "Fetch from a database" (`adit-gen --fetch`) are not part of the ADIT distribution. They are **what the user fetched from that site**.
ADIT records the source, ID, fetch date and time, the SHA-256 of the response and the license statement in spec.json (`provenance.fetched_structure`) and README.txt,
and copies the fetched file to the output directory. **If you publish or redistribute a file, follow the terms of that statement.**
The statements of each site were checked on 2026-09-19.

| Source | Statement | Page checked |
|---|---|---|
| PubChem | NCBI policy: information produced by the US government is in the public domain. Submitters may claim patent, copyright or other rights, and NCBI states that it cannot grant individual permissions. PubChem's guidance also tells users to check the submitter's license statement | https://www.ncbi.nlm.nih.gov/home/about/policies/ , https://pubchem.ncbi.nlm.nih.gov/docs/downloads |
| COD (Crystallography Open Database) | "All data in the COD and the database itself are dedicated to the public domain and licensed under the CC0 License" | https://www.crystallography.net/cod/ |
| Materials Project | "By downloading Content from Materials Project, User agrees to accept the Creative Commons Attribution 4.0 License … provided proper attribution is given to the Materials Project" | https://next-gen.materialsproject.org/about/terms (returns 403 to clients other than browsers, so it was checked at https://legacy.materialsproject.org/terms, which has the same text) |
| OPTIMADE (federated search) | Differs by provider. ADIT records the provider's name and home page, so check the terms there | https://providers.optimade.org/ |

PubChem's usage policy asks for no more than 5 requests per second. ADIT makes at most 3 requests per fetch (name → CID, SDF, compound name).

## Not included (cannot be distributed)

- **VASP POTCAR files**: These cannot be distributed to anyone other than license holders. ADIT **neither writes nor copies them**
- **Slater-Koster parameters (slakos/)**: Users install these from dftb.org. `slakos/`, the location for each set, is not in the repository (`.gitignore`). The only exception is the copies of mio-1-1 / tiorg-0-1 under `examples/**/skf/` listed in the table above
- **Pseudopotentials (UPF, psp8 and others)**: Same as above
- **The calculation codes themselves** (DFTB+, VASP, Gaussian, ORCA and others): Not included

## What goes into the executables (.exe / .app)

The distributions built with PyInstaller contain **a Japanese font (Noto Sans CJK JP, OFL-1.1)**.
**Its full text, `licenses/OFL-1.1-NotoSansCJK.txt`, is included with it** (`packaging/adit.spec` builds it this way.
If the full text is missing, the build stops without bundling the font).

The same distributions also contain `LICENSE`, this `THIRD_PARTY_NOTICES.md` and `licenses/*.txt`
(so that recipients can read the rights notices on their own machines).

| Included | License | Full text |
|---|---|---|
| Noto Sans CJK JP (font) | SIL Open Font License 1.1 | `licenses/OFL-1.1-NotoSansCJK.txt` (inside the executable, the same full text is also in `fonts/`) |
| Python, ASE, pydantic, Jinja2, matplotlib, SciPy, NumPy | The license of each distribution (PSF / LGPL-2.1+ / MIT / BSD-3 and others) | Included in each package's distribution |
| tomli-w | MIT | Same as above |
| PySide6 (for the GUI) | **LGPL-3.0** | See the Qt / PySide6 distribution |
| PyQtDarkTheme-fork (GUI theme) | MIT | Included in each package's distribution |
| pyte (screen of the Workspace terminal) | **LGPL-3.0** | Same as above |
| ptyprocess (terminal on Linux / macOS) | ISC | Same as above |
| pywinpty (terminal on Windows) | MIT | Same as above |
| RDKit (SMILES and Draw; only if installed at build time) | BSD-3-Clause | Same as above |
| spglib (space groups; only if installed at build time) | BSD-3-Clause | Same as above |

**PySide6 (Qt) and pyte are LGPL-3.0.** When you distribute an executable with the GUI, you must meet the LGPL requirement
that "users can replace Qt" (with PyInstaller's one-directory format,
Qt's shared libraries are included as separate files, so they can be replaced. **Avoid the single-file format (onefile)**).
`ONEFILE` in `packaging/adit.spec` is left as False.

## Trademarks

VASP, Gaussian, ORCA, GAMESS, Q-Chem, GRRM, Materials Studio and others are trademarks or registered trademarks of their respective owners.
ADIT **only writes input files for these programs and reads their output**. It is not affiliated with their providers
and has not been endorsed by them.
