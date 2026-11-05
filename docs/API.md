# Using ADIT from Python

This is a map for **calling ADIT as a library from your own scripts**, rather than through the commands
(`adit-gen` / `adit-analyze` / `adit-report` / `adit-convert`).

**Promise**: every name here can be imported directly from `adit.<module>`. Internal functions (those starting with `_`)
may change at any time, so do not use them.

## 1 The shortest path

```python
from pathlib import Path

from adit.config import load_config
from adit.spec import CalculationSpec
from adit.project import write_project
from adit.analysis.report import AnalysisOptions, run_analysis
from adit.report import load_run_report, methods_markdown

spec = CalculationSpec.load("spec.json")        # read the settings
cfg = load_config()                             # ~/.config/adit/cluster.toml (override with ADIT_CONFIG)
write_project(spec, cfg, Path("out/run1"))      # write the inputs, submit.sh and README.txt (nothing is submitted)

# --- a person runs the calculation here (ADIT does not submit jobs) ---

res = run_analysis("out/run1", AnalysisOptions(msd=True, out_dir=Path("scratch/analysis1")))
print(res.summary_text())                       # plots and tables are written to out_dir
print(methods_markdown([load_run_report("out/run1")], "en"))   # the methods section of a paper
```

## 2 Main functions

| What you want to do | What to call | Returns |
|---|---|---|
| Read or write settings | `adit.spec.CalculationSpec.load / save / from_json` | `CalculationSpec` (pydantic) |
| Check settings for mistakes | `adit.validate.validate(spec, cfg)` | `list[ValidationError]` (empty means no problems) |
| Build the inputs (without writing) | `adit.project.build_project(spec, cfg)` | `ProjectFiles` (`.texts` / `.copies`) |
| Write the inputs | `adit.project.write_project(spec, cfg, out_dir)` | `list[Path]` of the files written |
| Scan one value | `adit.scan.write_scan(spec, cfg, out, scans)` | The directories created |
| Batch over different structures | `adit.structures_batch.write_structures / write_enumerated` | The directories created |
| Read the output | `adit.analysis.readers.load_run(run_dir)` | `RunData` (energy, temperature, trajectory, final structure) |
| Analyze (plots and tables) | `adit.analysis.report.run_analysis(run_dir, opts)` | `AnalysisResult` |
| Check whether a run finished | `adit.results.summarize_run(run_dir)` / `not_run_yet` / `failure_lines` | `RunSummary` / `bool` / `list[str]` |
| Reports and reproducibility packages | `adit.report.load_run_report / methods_markdown / conditions_csv / results_csv / write_bundle / check_lines` | Strings and paths |
| Read back existing inputs | `adit.native_import.import_native(source, code=None)` | `NativeImportResult` (`.spec` is `None` when the input cannot be read; `.report()` holds the fields that could not be read) |
| Build or change a structure | `adit.structure.from_smiles / from_file / from_spacegroup`, `adit.convert` | `ase.Atoms` |
| Select atoms | `adit.analysis.select.select(atoms, "element O and z < 10")` | `numpy` indices (0-based) |
| Time series of distance, angle, dihedral, RMSD, Rg | `adit.analysis.geometry_series.*` | `numpy` arrays |
| Coordination number, centrosymmetry, Steinhardt, clusters, S(q) | `adit.analysis.local_order.*` | `numpy` arrays / `dict` |
| Count hydrogen bonds | `adit.analysis.hbond.count_series(frames, distance, angle)` | `dict` (array of counts per frame and the thresholds used; thresholds are required) |
| Hydrogen-bond lifetime and distance × angle distribution | `adit.analysis.hbond.lifetime(frames, distance, angle, dt_fs=None, tau_max=None, budget_mb=1024)` / `distance_angle_map(frames, rmax_A)` | `dict` of autocorrelation and lifetimes / `dict` of the histogram (`"distribution"` is a `Distribution2D`) |
| Compare a set (ΔE, plus ΔH, ΔS and ΔG if thermo.csv exists) | `adit.analysis.compare.analyze_compare(base, reactions=None)` / `parse_compare(text)` | `CompareResult` (`.reactions[i]["thermo"]` holds ΔH, ΔS and ΔG per temperature) / `list[Reaction]` |
| Read CREST conformers and weight them | `adit.analysis.crest.analyze_crest(crest_dir, out_dir, temperature_k=None)` / `read_ensemble` / `read_energies` / `boltzmann_weights(relative_ev, degeneracy, temperature_k)` | Conformer table (`dict`) / `list[(E [Eh], Atoms)]` / `dict` / array of weights (only when a temperature is given) |
| Velocity autocorrelation and vibrational spectrum | `adit.analysis.vacf.vacf(velocities, dt_fs, *, max_lag_fraction=0.5, source="velocities", window="hann")` (`velocities` is an array of (frames, atoms, 3), in Å/fs) | `VacfResult` |
| Effective mass | `adit.analysis.effective_mass.at_band_edges(kdist, energies, fermi)` | `list[EffectiveMass]` |
| Projected bands (QE projwfc.x) | `adit.analysis.projected_bands.read_filproj(path)` | `ProjectedBands` |
| Optics (QE epsilon.x) | `adit.analysis.optical.read_epsilon(run_dir)` | `Optical` (n, k, absorption coefficient, reflectivity, EELS) |
| Bader charges (ACF.dat from the external bader program) | `adit.analysis.bader.read_acf(path)` | `BaderCharges` (does not compute the partitioning) |
| Send a heavy analysis to a compute machine | `adit.analysis.heavy_setup.write_job(run_dir, *, dt_fs, ...)` (`dt_fs` is required) | `list[Path]` of the files written (msd_worker.py and msd_run.sh) |
| VASP PROCAR | `adit.analysis.procar.read_procar(path)` | `Procar` (projections, eigenvalues, k-points) |
| Gaussian, GAMESS and Q-Chem output | `adit.analysis.readers_qc.read_gaussian / read_gamess / read_qchem` | `QcOutput` (energy, structure, frequencies, IR, Raman, charges) |
| OpenMX output | `adit.analysis.readers_openmx.read_openmx(path)` | `OpenmxOutput` |
| VASP WAVEDER | `adit.analysis.waveder.read_waveder(path)` / `dielectric_imag(...)` | `Waveder` / array of eps2 |
| Displacement and local strain | `adit.analysis.displacement.displacements / local_strain` | `numpy` arrays / `LocalStrain` |
| Free-energy surface and joint distribution | `adit.analysis.free_energy.free_energy_surface / distribution_2d / read_colvar` | `FreeEnergySurface` / `Distribution2D` |
| Voronoi volumes | `adit.analysis.voronoi.voronoi(positions, cell)` | `VoronoiResult` |
| Solvent-accessible surface area | `adit.analysis.sasa.sasa(positions, radii, probe)` | `SasaResult` |
| Principal component analysis and clustering of frames | `adit.analysis.frames_pca.pca / kmeans` | `PcaResult` / `ClusterResult` |

## 3 Exceptions

All exceptions inherit from `adit.errors.AditError`. **This lets one line separate external causes (a missing file, out of memory) from
failures that come from ADIT.**

```python
from adit.errors import AditError, AditValueError

try:
    write_project(spec, cfg, out)
except AditError as ex:          # ADIT decided this cannot work and stopped
    print(f"ADIT: {ex}")
```

`AditValueError` inherits from both `AditError` and `ValueError`, so existing code written with `except ValueError`
keeps working.

| Exception | When | Module |
|---|---|---|
| `GenerationError` | The inputs cannot be built | `adit.codes.base` |
| `OutputNotEmpty` | The destination already has files | `adit.project` |
| `StructureError` | A structure cannot be built or read | `adit.structure` |
| `ScanError` | The values to scan cannot be read | `adit.scan` |
| `StructuresError` / `EnumerateError` | Batch generation, substituent enumeration | `adit.structures_batch` / `adit.enumerate_r` |
| `TemplateError` | Saving or loading a template | `adit.templates` |
| `ConfigMissing` / `ConfigError` | The settings file | `adit.config` |
| `ReportError` | A report cannot be made (no spec.json, for example) | `adit.report` |
| `TrajectoryTooLarge` | A trajectory exceeds the memory limit | `adit.analysis.trajectory` |
| `TransportError` / `VolumetricError` / `XrdError` / `RateError` | Viscosity, volumetric data, powder diffraction, rate constants | `adit.analysis.*` |
| `HandoffError` | Handoff between stages (**only this module does not inherit from `AditError`**) | `adit.handoff` |

`adit/handoff.py` is copied onto the compute machine and run on its own, so it does not import ADIT. It defines its exception itself.

## 4 Things to keep in mind

- **Jobs are not submitted.** `write_project` only writes `submit.sh`, and a person runs it
- **Analysis writes to the calculation directory by default.** To keep that directory clean, set `AnalysisOptions(out_dir=...)`
- **There is one settings file.** To replace it (in tests, for example), use the environment variable `ADIT_CONFIG`
 (you can also pass it directly as `load_config(path)`)
- **`adit.lang.set_language("ja" / "en")`** switches the language. Error messages and table headings read the same setting
- **Do not import the GUI.** `adit.gui` requires PySide6. The web version, the CLI and the analysis code do not import it
