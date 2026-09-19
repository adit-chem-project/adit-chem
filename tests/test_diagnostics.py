from pathlib import Path

import pytest

from adit.analysis import diagnostics as dg
from adit.analysis.report import AnalysisResult
from adit.results import summarize_run

REPO = Path(__file__).resolve().parent.parent

# One typical message per code, in the file the code writes it to.
CASES = [
    ("dftbplus", "output.log", "SCC is NOT converged, maximal SCC iterations exceeded", dg.SCF),
    ("dftbplus", "output.log", "!!! Geometry did NOT converge!", dg.GEOMETRY),
    ("xtb", "output.log", "-1- scf: Self consistent charge iterator did not converge", dg.SCF),
    ("xtb", "output.log", "   *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 500 ITERATIONS ***", dg.GEOMETRY),
    ("espresso", "output.log", "     convergence NOT achieved after 100 iterations: stopping", dg.SCF),
    ("espresso", "output.log", "     Maximum CPU time exceeded", dg.WALLTIME),
    ("espresso", "output.log", "     S matrix not positive definite", dg.DIAGONALIZATION),
    ("espresso", "output.log", "     bfgs failed after  12 scf cycles and  10 bfgs steps, convergence not achieved", dg.GEOMETRY),
    ("vasp", "output.log", "ZBRENT: fatal error in bracketing", dg.GEOMETRY),
    ("vasp", "OUTCAR", " Error EDDDAV: Call to ZHEGV failed. Returncode = 12 2 16", dg.DIAGONALIZATION),
    ("orca", "output.log", "[file orca_tools/qcmem.cpp, line 949, Process 7]:  OUT OF MEMORY ERROR!", dg.MEMORY),
    ("orca", "output.log", "                                 INPUT ERROR\n            UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE", dg.INPUT),
    ("cp2k", "output.log", " *** SCF run NOT converged ***", dg.SCF),
    ("cp2k", "output.log", " *** MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED ***", dg.GEOMETRY),
    ("lammps", "log.lammps", "ERROR: Lost atoms: original 500 current 498 (src/thermo.cpp:488)", dg.LOST_ATOMS),
    ("lammps", "log.lammps", "ERROR on proc 0: Out of range atoms - cannot compute PPPM (src/KSPACE/pppm.cpp:1900)", dg.UNSTABLE),
    ("lammps", "log.lammps", "ERROR: Illegal velocity command (velocity.cpp:78)", dg.ERROR),
    ("gromacs", "md.log", "Too many LINCS warnings (1001)", dg.UNSTABLE),
    ("gromacs", "output.log", "Fatal error:\nThere is no domain decomposition for 8 ranks", dg.ERROR),
    ("gromacs", "grompp.log", "Error in user input:\nInvalid keyword", dg.INPUT),
]


def _run_dir(tmp_path: Path, code: str, name: str, body: str) -> Path:
    d = tmp_path / code
    d.mkdir(exist_ok=True)
    (d / name).parent.mkdir(parents=True, exist_ok=True)
    (d / name).write_text("line one\nall fine\n" + body + "\nlast line\n", encoding="utf-8")
    return d


@pytest.mark.parametrize("code,name,body,cause", CASES)
def test_known_failure_messages_are_found_with_their_line(tmp_path, code, name, body, cause):
    d = _run_dir(tmp_path, code, name, body)
    diag = dg.scan_diagnostics(d, code)
    assert diag.supported and name in diag.files
    assert diag.causes() == [cause]
    hit = diag.errors[0]
    assert hit.file == name and hit.line == 3 and hit.text == body.splitlines()[0].strip()
    lines = diag.summary_lines()
    assert len(lines) == 1 and dg.cause_label(cause) in lines[0] and f"{name}:3" in lines[0]


def test_nothing_is_reported_for_a_clean_log(tmp_path):
    d = _run_dir(tmp_path, "espresso", "output.log", "!    total energy = -22.8 Ry\n     JOB DONE.")
    diag = dg.scan_diagnostics(d, "espresso")
    assert diag.errors == [] and diag.n_warnings == 0 and diag.summary_lines() == []
    assert diag.as_dict()["causes"] == [] and diag.as_dict()["supported"] is True


def test_warnings_are_counted_and_the_first_three_quoted(tmp_path):
    body = "\n".join(f"WARNING: something {i}" for i in range(5))
    d = _run_dir(tmp_path, "lammps", "log.lammps", body)
    diag = dg.scan_diagnostics(d, "lammps")
    assert diag.n_warnings == 5 and diag.errors == []
    line = diag.summary_lines()[0]
    assert "5" in line and "log.lammps:3" in line and "log.lammps:5" in line and "log.lammps:6" not in line
    assert "2" in line.rsplit("(", 1)[-1] or "ほか 2" in line or "2 more" in line


def test_lincs_warning_is_a_warning_not_a_failure(tmp_path):
    d = _run_dir(tmp_path, "gromacs", "md.log", "Step 12, time 0.024 (ps)  LINCS WARNING\nrelative constraint deviation after LINCS:")
    diag = dg.scan_diagnostics(d, "gromacs")
    assert diag.errors == [] and diag.n_warnings == 1 and diag.warnings[0].cause == dg.UNSTABLE


def test_dftbplus_bare_marker_takes_the_next_line(tmp_path):
    d = _run_dir(tmp_path, "dftbplus", "output.log", "WARNING!\n-> Current stacksize not set to unlimited or hard limit\n\nERROR!\n-> Unable to open file")
    diag = dg.scan_diagnostics(d, "dftbplus")
    assert diag.n_warnings == 1 and "stacksize" in diag.warnings[0].text
    assert diag.causes() == [dg.ERROR] and "Unable to open file" in diag.errors[0].text


@pytest.mark.parametrize("name,body,cause", [
    ("water.o12345", "=>> PBS: job killed: walltime 3650 exceeded limit 3600", dg.WALLTIME),
    ("water.o12345", "=>> PBS: job killed: mem 9000000kb exceeded limit 8000000kb", dg.MEMORY),
    ("slurm-777.out", "slurmstepd: error: *** JOB 777 ON node1 CANCELLED AT 2026-09-19T10:00:00 DUE TO TIME LIMIT ***", dg.WALLTIME),
    ("slurm-777.out", "slurmstepd: error: Detected 1 oom_kill event in StepId=777.batch. Some of the step tasks have been OOM Killed.", dg.MEMORY),
    ("output.log", "bash: line 1: 4242 Killed                  dftb+", dg.SIGNAL),
    ("output.log", "MPI_ABORT was invoked on rank 3 in communicator MPI_COMM_WORLD", dg.MPI),
    ("output.log", "application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0", dg.MPI),
    ("output.log", "= BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES", dg.MPI),
])
def test_scheduler_and_mpi_messages_are_found_for_every_code(tmp_path, name, body, cause):
    d = _run_dir(tmp_path, "nwchem", name, body)
    diag = dg.scan_diagnostics(d, "nwchem")
    assert diag.supported is False and diag.causes() == [cause] and name in diag.files


def test_the_code_is_read_from_spec_json_or_the_inputs(tmp_path):
    d = _run_dir(tmp_path, "x", "output.log", "     convergence NOT achieved after 100 iterations: stopping")
    (d / "pw.in").write_text("&CONTROL\n/\n", encoding="utf-8")
    assert dg.scan_diagnostics(d).code == "espresso"
    (d / "spec.json").write_text('{"method": {"code": "cp2k"}}', encoding="utf-8")
    assert dg.scan_diagnostics(d).code == "cp2k"


def test_summarize_run_carries_the_diagnostics():
    s = summarize_run(REPO / "examples" / "water_generated")
    assert len(s.diagnostics) == 1 and "output.log:40" in s.diagnostics[0] and s.diagnostics[0] in s.text()


def test_analysis_summary_prints_the_diagnostics_table():
    d = dg.Diagnostics(code="lammps", supported=True, files=["log.lammps"], n_warnings=0,
                       errors=[dg.Hit("log.lammps", 9, "ERROR: Lost atoms: original 2 current 1", dg.LOST_ATOMS)])
    res = AnalysisResult(code="lammps", run_dir="/tmp/x", tables={"diagnostics": d.as_dict()})
    text = res.summary_text()
    assert dg.cause_label(dg.LOST_ATOMS) in text and "log.lammps:9" in text
    assert AnalysisResult(code="lammps", run_dir="/tmp/x").summary_text().count("\n") == 0


@pytest.mark.parametrize("code", dg.supported_codes())
def test_remedies_carry_a_source_url_on_every_item(code):
    lines = dg.remedy_lines(code)
    items = [x for x in lines if x.startswith("  - ")]
    assert len(items) >= 3
    assert all("http" in x or "出力そのもの" in x for x in items), code
    assert not any("推奨" in x and "ADIT の推奨ではありません" not in x for x in lines)
    assert dg.remedy_lines("nwchem") == []
