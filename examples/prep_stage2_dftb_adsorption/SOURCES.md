# Run of a set of runs to compare (adsorption) — DFTB+ 25.1, mio-1-1 (2026-09-12)

## What this example is

An adsorption set with 1 water molecule placed 3 Å above 2×2 graphene (8 C atoms, 15 Å vacuum in the box, periodic).
Generated with `adit-gen spec.json out/ --compare-set set.json` (`src/adit/compare_sets.py`), run with DFTB+ in each directory,
and collected with `adit-analyze out/ --compare` (`src/adit/analysis/compare.py`).

- Settings (the same for all 3): DFTB+ 25.1, SK set mio-1-1, SCC, single point, k-points 4×4×1 (for the periodic ones)
- The molecule keeps the default `as_is` for `molecule_box`, so it is non-periodic (no k-points). ADIT does not stop on this difference; it is recorded in `compare.json` under
  `differences` (structure.periodic) and `partial` (the k-point items). This example exists to check that
- Input structures were built from ASE `graphene(size=(2, 2, 1), vacuum=7.5)` and `molecule("H2O")` (ASE 3.29.0)

## Files

| File | Contents |
|---|---|
| `inputs/` | spec.json, set.json and structures (slab.extxyz, h2o.xyz, adsorbed.extxyz) used for generation |
| `compare.json` | Definition of the set (reactions), settings of each run, differences in settings, balance of composition |
| `slab/ molecule/ adsorbed/` | Each run (input and DFTB+ output; skf files copied from mio-1-1, CC BY-SA 4.0) |
| `compare_*.csv`, `compare_summary.json`, `compare_energy.png` | Results of the analysis (adit-analyze --compare) |

## Results (numbers only; no judgment on whether they can be compared)

| Run | Atoms | Total energy [eV] |
|---|---|---|
| adsorbed | 11 | -488.6238 |
| slab | 8 | -377.7017 |
| molecule | 3 | -110.9604 |

ΔE = E(adsorbed) − E(slab) − E(molecule) = +0.0383 eV (+3.695 kJ/mol). 1 setting differs (structure.periodic),
and 6 settings are present in only some runs (the k-point items and similar).
