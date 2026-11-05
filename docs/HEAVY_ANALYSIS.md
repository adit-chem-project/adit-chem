# Heavy analysis: processing large trajectories on a compute machine

ADIT does not submit calculations. **Heavy analyses get the same treatment.**

```
ADIT prepares          a person runs             ADIT plots
msd_worker.py   ──▶   bash msd_run.sh    ──▶   adit-analyze <directory>
msd_run.sh           (→ msd_vanhove.json)     (plots, summary, report)
```

This is the same approach as phonons (`phonon_collect.py`) and elastic constants. Calculations that take tens of minutes do not run inside the GUI, so
**the GUI does not freeze and a laptop's memory is not exhausted.**

## What counts as "heavy"

Before running, `adit-analyze --vanhove` **makes an arithmetic estimate** (it does not actually compute anything).

```
displacements ≈ number of τ × frames × atoms × 0.75
time ≈ displacements × 100 ns + number of τ × 5 ms     (coefficients measured locally on 2026-09-13)
coordinate memory = frames × atoms × 24 bytes
```

If the estimate exceeds **60 seconds** (change this with `--heavy-limit`), ADIT does not compute on the spot. It writes files for running the job instead.

| System | Displacements | Estimate | What happens |
|---|---|---|---|
| 18 atoms × 1,000 frames × 100 τ | 1.35 million | 0.6 s | Computed on the spot |
| 18 atoms × 11,001 frames × 100 τ | 14.85 million | 2 s | Computed on the spot |
| 1,000 atoms × 10,000 frames × 100 τ | 750 million | 76 s | Files are written |
| 18 atoms × 11,001 frames × 3,000 τ | 450 million | 45 s + 15 s for the fit | Files are written |

## Usage

```bash
adit-analyze run/ --msd --vanhove              # estimate; if heavy, write msd_worker.py and msd_run.sh
bash run/msd_run.sh                             # a person runs this (on a cluster, submit it as a job)
adit-analyze run/ --msd --vanhove              # second pass: read msd_vanhove.json and plot it
adit-analyze run/ --msd --vanhove --vanhove-here   # small systems may also be computed on the spot
adit-analyze run/ --write-msd-job              # write the files without estimating (default location is run/)
```

The contents of `msd_run.sh` are filled in from the output: the trajectory file name (geo_end.xyz, xtb.trj, XDATCAR and so on) and the time per frame
(`--dt`; if it cannot be read, the blank stays and a comment says so). For periodic systems whose trajectory file has no cell
(DFTB+ geo_end.xyz and others), the cell of the last frame is written to `msd_cell.txt` (3 TV lines) and passed with `--cell`.
Without a cell, `msd_worker.py` cannot unwrap boundary crossings, and it says so in `warnings` in the JSON.

`msd_worker.py` **runs with only numpy and scipy** (it does not import ADIT; this is the same rule as for `handoff.py`, and the test
`tests/test_heavy_setup.py` enforces it). It runs even on a cluster where ADIT is not installed.

```bash
# reading a custom format as is, with one line per atom (label x y z)
python msd_worker.py Na.dat --natoms 18 --cell dftb.inp --dt 10 --taus 100 --out msd_vanhove.json

# formats ASE can read (extxyz, vasprun.xml and so on) work as they are
python msd_worker.py trajectory.extxyz --dt 10 --species Na
```

## What you get

`msd_vanhove.json` holds **only numbers** (no judgments).

| Key | Contents |
|---|---|
| `d_per_atom_cm2_s` | D for each atom, from the slope of the MSD using all time origins |
| `d_from_mean_msd_cm2_s` | D from the MSD averaged over all atoms |
| `vanhove.d_fit_A2_fs` | D(τ) from a Gaussian fit to the displacement distribution |
| `vanhove.d_direct_A2_fs` | D(τ) from ⟨r²⟩/(6τ) with the same displacements. For a Gaussian, this agrees with the row above |
| `vanhove.alpha2` | Non-Gaussian parameter α₂ = 3⟨r⁴⟩/(5⟨r²⟩²) − 1 (Rahman 1964). 0 means Gaussian |
| `vanhove.truncated_from_fs` | Lag time at which displacements begin to saturate at the minimum-image limit (null if none) |
| `cell_A` | Cell used for unwrapping [Å] (3 × 3; null if none) |
| `warnings` | Cautions, such as when there is no cell and nothing was unwrapped (for a periodic system, pass the cell with `--cell`) |

**Watch for minimum-image saturation.** `--vanhove-displacement mic` (the default) uses the minimum image of wrapped coordinates, so it cannot measure
displacements larger than **half the smallest cell width**. When 1 % of the displacements exceed 90 % of that limit, that lag time is recorded in
`truncated_from_fs`, and the summary also warns about it. `--vanhove-displacement unwrapped` does not saturate.

## What ADIT does not judge

- Which D to adopt (MSD slope / distribution fit / ⟨r²⟩/(6τ))
- Where "non-Gaussian" begins
- Which atoms are "fast"
- Whether the fit range is appropriate (the default is 10-50 % of the lag time; change it with `--msd-fit`)
