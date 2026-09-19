
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from adit.lang import L
from adit.spec import CalculationSpec
from adit.validate_types import ValidationError

PLUMED_FILE = "plumed.dat"
PLUMED_LOG = "plumed.log"
SUPPORTED_CODES = ("lammps", "gromacs", "openmm")
DOC = "https://www.plumed.org/doc-v2.10/user-doc/html/index.html"


def plumed_text(spec: CalculationSpec) -> str | None:
    p = spec.plumed
    if p is None:
        return None
    if p.input_file.strip():
        path = Path(p.input_file).expanduser()
        return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else None
    text = p.lines
    return text if text.endswith("\n") else text + "\n"


def engine_has_plumed(code: str) -> bool | None:
    if code == "lammps":
        exe = shutil.which("lmp") or shutil.which("lmp_serial") or shutil.which("lmp_mpi")
        args, needle = ([exe, "-h"], "PLUMED") if exe else (None, "")
    elif code == "gromacs":
        exe = shutil.which("gmx") or shutil.which("gmx_mpi")
        args, needle = ([exe, "mdrun", "-h"], "-plumed") if exe else (None, "")
    else:
        return None
    if args is None:
        return None
    try:
        out = subprocess.run(args, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return needle in (out.stdout + out.stderr)


def check_syntax(text: str, n_atoms: int) -> str | None:
    exe = shutil.which("plumed")
    if not exe:
        return None
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / PLUMED_FILE
        path.write_text(text, encoding="utf-8")
        try:
            out = subprocess.run([exe, "--no-mpi", "driver", "--plumed", str(path), "--natoms", str(n_atoms),
                                  "--parse-only", "--ixyz", "/dev/null"],
                                 capture_output=True, text=True, timeout=120, cwd=tmp)
        except (OSError, subprocess.SubprocessError):
            return None
        if out.returncode == 0:
            return None
        message = (out.stdout + out.stderr).strip().splitlines()
        detail = " / ".join(line.strip() for line in message[-6:] if line.strip())
        return detail or L("理由は分かりません", "no reason was reported")


def validate_plumed(spec: CalculationSpec) -> list[ValidationError]:
    p = spec.plumed
    if p is None:
        return []
    errs: list[ValidationError] = []
    code = spec.method.code
    if code not in SUPPORTED_CODES:
        errs.append(ValidationError("plumed", L(
            f"PLUMED をつなげるのは {' / '.join(SUPPORTED_CODES)} の MD だけです ({code} は対応していません)",
            f"PLUMED can be attached to {' / '.join(SUPPORTED_CODES)} MD only (not {code})")))
    if spec.task.type != "molecular_dynamics":
        errs.append(ValidationError("task.type", L(
            "PLUMED は MD (分子動力学) につなぎます。ほかの計算の種類では使えません",
            "PLUMED attaches to molecular dynamics; it cannot be used with other task types")))
    has_file, has_lines = bool(p.input_file.strip()), bool(p.lines.strip())
    if has_file == has_lines:
        errs.append(ValidationError("plumed.input_file", L(
            "PLUMED の入力は、ファイル (input_file) か、そのままの文 (lines) のどちらか一方を指定してください。ADIT は集合変数を作りません",
            "give the PLUMED input either as a file (input_file) or as text (lines), not both and not neither; ADIT does not invent collective variables")))
    if has_file and not Path(p.input_file).expanduser().is_file():
        errs.append(ValidationError("plumed.input_file", L(
            f"ファイルがありません: {p.input_file}", f"file not found: {p.input_file}")))
    text = plumed_text(spec)
    if text is not None and errs == []:
        reason = check_syntax(text, len(spec.structure.atoms.symbols))
        if reason:
            errs.append(ValidationError("plumed.lines", L(
                f"PLUMED がこの入力を読めません: {reason}", f"PLUMED cannot parse this input: {reason}")))
    if code in ("lammps", "gromacs") and engine_has_plumed(code) is False:
        errs.append(ValidationError("plumed", L(
            f"この PC の {code} は PLUMED を組み込んでいません (PLUMED を組み込んだビルドが要ります)。生成はできますが、この PC では走りません",
            f"the {code} on this PC was built without PLUMED (a PLUMED-enabled build is needed); the input can still be generated, but it will not run here")))
    return errs


def readme_lines(spec: CalculationSpec) -> list[str]:
    if spec.plumed is None:
        return []
    code = spec.method.code
    how = {"lammps": L("LAMMPS は PLUMED パッケージを入れてビルドしたものが要ります (fix plumed)。",
                       "LAMMPS must be built with the PLUMED package (fix plumed)."),
           "gromacs": L("GROMACS は PLUMED を使えるビルドが要ります (mdrun -plumed)。GROMACS 2026.3 では、-plumed があっても "
                        "実行時に環境変数 PLUMED_KERNEL が libplumedKernel.so を指していないと「not available」で止まることがあります。",
                        "GROMACS must support PLUMED (mdrun -plumed). With GROMACS 2026.3, having -plumed may not be enough: "
                        "it stopped with 'not available' unless the PLUMED_KERNEL environment variable pointed at libplumedKernel.so."),
           "openmm": L("実行する環境に openmm-plumed が要ります (conda install -c conda-forge openmm-plumed)。",
                       "the environment needs openmm-plumed (conda install -c conda-forge openmm-plumed).")}.get(code, "")
    return [L(f"  PLUMED の入力 {PLUMED_FILE} は利用者が書いたものです。ADIT は集合変数もバイアスも作らず、内容も確かめていません "
              "(読めるかどうかだけ、plumed があるときに確かめています)。",
              f"  The PLUMED input {PLUMED_FILE} is yours. ADIT invents no collective variables or bias and does not check the physics "
              "(only whether PLUMED can parse it, and only when plumed is installed)."),
            "  " + how,
            L("  PLUMED の既定の単位は nm・kJ/mol・ps で、計算コードの単位系とは別です (入力の UNITS で変えられます)。原子の番号は 1 始まりです。",
              "  PLUMED's default units are nm, kJ/mol and ps, independent of the engine's unit system (change them with UNITS in the input). PLUMED numbers atoms from 1."),
            L(f"  書き方と機能の説明は {DOC}。", f"  For the syntax and the available features see {DOC}.")]


def output_lines(spec: CalculationSpec) -> list[str]:
    if spec.plumed is None:
        return []
    lines = [L(f"  {PLUMED_FILE} が書き出すファイル (COLVAR、HILLS など) の名前と中身は、その入力しだいです。",
               f"  The files written by {PLUMED_FILE} (COLVAR, HILLS, ...) and their contents depend on that input.")]
    if spec.method.code == "lammps":
        lines.append(L(f"  {PLUMED_LOG}   PLUMED 自身のログ (fix plumed の outfile)",
                       f"  {PLUMED_LOG}   the PLUMED log (outfile of fix plumed)"))
    lines.append(L("  ADIT の解析は、PLUMED の COLVAR 形式のファイル (先頭が #! FIELDS) があれば、その列を読んで要約に入れます。",
                   "  ADIT's analysis reads any PLUMED COLVAR-format file (starting with #! FIELDS) and lists its columns in the summary."))
    return lines
