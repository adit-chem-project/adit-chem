# Contributing

## Running it

```bash
python -m pip install -e ".[dev,gui,smiles,analysis]"
pytest -q                     # all tests (the GUI tests run without a display)
ruff check src tests packaging   # undefined names and unused imports (also runs in CI)
```

If calculation programs (DFTB+, xtb, Quantum ESPRESSO, CP2K, LAMMPS, GROMACS...) are on PATH,
the tests that actually run them also run. If not, only those tests are skipped.

## Three principles ADIT follows

1. **No logic in the GUI.** The GUI only assembles the settings (`CalculationSpec`).
   Input files are created by `codes/` and `project.py`
2. **ADIT does not submit jobs.** A person types `qsub` / `sbatch`. The only built-in way to run is directly on the local machine
3. **ADIT does not second-guess the user's chemical judgment.** It only reports mechanical facts that mean "this will not run as it is".
   It does not say whether a run converged or whether a value is reasonable

## Writing rules

- **Check keywords and formats against the official documentation or real output before writing them.** Do not write them from memory
  (record the source of what you checked in `tests/data/SOURCES.md`)
- **Offer the same features in all three entry points (command line, desktop app, web app).**
  When you add an analysis option, change `gui/analysis_fields.py` + `gui/panels/analysis_panel.py` + `gui/help.py` +
  `gui/i18n.py` + `web/templates/analysis.html` together (`tests/test_analysis_ui_parity.py` watches this)
- **Write on-screen text in both Japanese and English.** Always add the description in `gui/help.py` and the translation in `gui/i18n.py`
- When you show a number, give its **unit and source**

## Adding a calculation code

Implement `InputGenerator` in `codes/<name>.py` and call `register()`, add the Method to `spec.py`,
and add the desktop panel (`gui/panels/`) and the web form.

## Tests

Give every new feature a test on **a small example whose answer is known**
(for example: coordination number 12 for FCC, m*/m_e = 1 for free electrons, surface area 4πr² for an isolated atom).

- **Tests do not use `examples/` in place.** When it would be the output location of an analysis or generation, copy it to `tmp_path` first
  (running in place rewrites tracked files). If `git status` shows changes in `examples/` after `pytest`,
  some test is using an example in place
- Write figures so they are reproducible (no date or time in SVG/PDF). Some tests assume that the same input gives the same bytes
