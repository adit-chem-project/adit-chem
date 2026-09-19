"""Small, explicit molecular-input subsets for Gaussian, US GAMESS, Q-Chem and GRRM17."""
# Official/vendor input references:
# https://www.conflex.co.jp/gaussian_support/input.php
# https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/intro.pdf
# https://www.msg.chem.iastate.edu/GAMESS/GAMESS_Manual/input.pdf
# https://manual.q-chem.com/6.0/Ch3.S2.SS1.html
# https://afir.sci.hokudai.ac.jp/manual/grrm17/grrm17_12.html
# https://afir.sci.hokudai.ac.jp/manual/grrm17/grrm17_9.html

from __future__ import annotations

import re
from pathlib import Path

from ase.data import atomic_numbers

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, GaussianMethod, GamessMethod, QchemMethod, GrrmMethod
from adit.validate import electron_parity_error
from adit.validate_types import ValidationError


_TOKEN = re.compile(r"[A-Za-z0-9_+*().,-]+\Z")


def _molecule_errors(spec: CalculationSpec, *, tasks: tuple[str, ...], code: str) -> list[ValidationError]:
    st = spec.structure
    errors: list[ValidationError] = []
    if st.periodic:
        errors.append(ValidationError("structure.atoms", L(
            f"この {code} 生成器は分子だけに対応します", f"this {code} generator supports molecules only")))
    if st.fixed_atoms or st.fixed_axes:
        errors.append(ValidationError("structure.fixed_atoms", L(
            f"この {code} 生成器は固定原子・固定軸を入力へ写せません",
            f"this {code} generator cannot write fixed atoms or axes")))
    if st.multiplicity != 1:
        errors.append(ValidationError("structure.multiplicity", L(
            f"この {code} 生成器は一重項だけに対応します。開殻条件はまだ生成しません",
            f"this {code} generator supports singlets only; open-shell input is not yet generated")))
    if spec.task.type not in tasks:
        errors.append(ValidationError("task.type", L(
            f"この {code} 生成器が対応する計算種類: {', '.join(tasks)}",
            f"this {code} generator supports only: {', '.join(tasks)}")))
    if spec.handoff is not None:
        errors.append(ValidationError("handoff", L(
            f"この {code} 生成器は前の計算を引き継げません",
            f"this {code} generator cannot carry over a previous run")))
    electrons = sum(atomic_numbers[s] for s in st.atoms.symbols) - st.charge
    parity = electron_parity_error(electrons, st.charge, st.multiplicity)
    if parity:
        errors.append(parity)
    return errors


def _word_error(value: str, location: str, label: str) -> list[ValidationError]:
    if value and _TOKEN.fullmatch(value):
        return []
    return [ValidationError(location, L(
        f"{label}を改行・空白のない 1 語で明示してください",
        f"explicitly enter {label} as one token without spaces or newlines"))]


def _parallel_error(spec: CalculationSpec, cfg: Config, code: str) -> list[ValidationError]:
    if ((spec.runtime.mpiprocs > 1 or spec.runtime.omp_threads > 1)
            and spec.runtime.profile in cfg.profiles
            and not cfg.profile(spec.runtime.profile).commands.get(code)):
        return [ValidationError("runtime.profile", L(
            f"並列実行では環境設定の commands.{code} に起動コマンドを明示してください",
            f"for parallel execution, set commands.{code} explicitly in your settings"))]
    return []


def _format_command(spec: CalculationSpec, profile: Profile, code: str, default: str) -> str:
    return profile.command_for(code, default).format(
        mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")


class GaussianGenerator(InputGenerator):
    code = "gaussian"
    uses_kpoints = False
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, GaussianMethod):
            return [ValidationError("method.code", L("Gaussian の条件ではありません", "not a Gaussian method"))]
        errors = _molecule_errors(spec, tasks=("single_point",), code="Gaussian")
        errors += _word_error(m.theory, "method.theory", L("手法", "method"))
        errors += _word_error(m.basis, "method.basis", L("基底関数", "basis set"))
        if m.basis.lower() in {"gen", "genecp", "chk"}:
            errors.append(ValidationError("method.basis", L(
                "Gen/GenECP/Chk には追加の基底データが要るため、この限定生成器では扱えません",
                "Gen/GenECP/Chk needs an additional basis section and is outside this limited generator")))
        errors += _parallel_error(spec, cfg, self.code)
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st = spec.method, spec.structure
        lines = [f"# {m.theory}/{m.basis} SP", "", "ADIT calculation", "", f"{st.charge} {st.multiplicity}"]
        lines += [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        return {"gaussian.gjf": "\n".join(lines + ["", ""])}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = _format_command(spec, profile, self.code, "g16")
        return f"{exe} < gaussian.gjf > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="g16", files=[L("  gaussian.gjf  Gaussian の分子一点計算の入力", "  gaussian.gjf  Gaussian molecular single-point input")],
                           prepare=[L("  Gaussian の利用権、指定した手法・基底関数、実行コマンドを利用者が確認してください。",
                                      "  Check your Gaussian license, selected method and basis, and launch command.")],
                           outputs=[L("  output.log    Gaussian の標準出力。ADIT はまだ内容を解析しません。",
                                      "  output.log    Gaussian output; ADIT does not yet parse it.")])


class GamessGenerator(InputGenerator):
    code = "gamess"
    uses_kpoints = False
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, GamessMethod):
            return [ValidationError("method.code", L("US GAMESS の条件ではありません", "not a US GAMESS method"))]
        errors = _molecule_errors(spec, tasks=("single_point",), code="US GAMESS")
        if m.gbasis.upper() not in {"STO", "N21", "N31", "N311"}:
            errors.append(ValidationError("method.gbasis", L(
                "GBASIS は STO、N21、N31、N311 のいずれかを明示してください",
                "explicitly choose GBASIS from STO, N21, N31, N311")))
        if not 1 <= m.ngauss <= 6:
            errors.append(ValidationError("method.ngauss", L(
                "NGAUSS を 1～6 の整数で明示してください", "explicitly set NGAUSS to an integer from 1 to 6")))
        errors += _parallel_error(spec, cfg, self.code)
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st = spec.method, spec.structure
        lines = [f" $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG={st.charge} MULT=1 UNITS=ANGS $END",
                 f" $BASIS GBASIS={m.gbasis.upper()} NGAUSS={m.ngauss} $END",
                 " $DATA", "ADIT calculation", "C1"]
        lines += [f" {s} {atomic_numbers[s]:.1f} {x:.10f} {y:.10f} {z:.10f}"
                  for s, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        return {"gamess.inp": "\n".join(lines + [" $END", ""])}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = _format_command(spec, profile, self.code, "rungms")
        return f"{exe} gamess > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="rungms", files=[L("  gamess.inp   US GAMESS の RHF 一点計算の入力", "  gamess.inp   US GAMESS RHF single-point input")],
                           prepare=[L("  US GAMESS の利用登録と基底関数の対応を確認し、必要なら commands.gamess に起動コマンドを指定してください。",
                                      "  Check your US GAMESS registration and basis support; set commands.gamess if your launch command differs.")],
                           outputs=[L("  output.log   GAMESS の標準出力。ADIT はまだ内容を解析しません。",
                                      "  output.log   GAMESS output; ADIT does not yet parse it.")])


class QchemGenerator(InputGenerator):
    code = "qchem"
    uses_kpoints = False
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, QchemMethod):
            return [ValidationError("method.code", L("Q-Chem の条件ではありません", "not a Q-Chem method"))]
        errors = _molecule_errors(spec, tasks=("single_point",), code="Q-Chem")
        errors += _word_error(m.theory, "method.theory", L("手法", "method"))
        errors += _word_error(m.basis, "method.basis", L("基底関数", "basis set"))
        errors += _parallel_error(spec, cfg, self.code)
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st = spec.method, spec.structure
        lines = ["$molecule", f"{st.charge} {st.multiplicity}"]
        lines += [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        lines += ["$end", "", "$rem", "JOBTYPE SP", f"METHOD {m.theory}", f"BASIS {m.basis}", "$end", ""]
        return {"qchem.in": "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = _format_command(spec, profile, self.code, "qchem")
        return f"{exe} qchem.in output.log"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="qchem", files=[L("  qchem.in     Q-Chem の分子一点計算の入力", "  qchem.in     Q-Chem molecular single-point input")],
                           prepare=[L("  Q-Chem のライセンスと、指定した METHOD/BASIS が利用できるか確認してください。",
                                      "  Check your Q-Chem license and availability of the selected METHOD/BASIS.")],
                           outputs=[L("  output.log   Q-Chem の出力。ADIT はまだ内容を解析しません。",
                                      "  output.log   Q-Chem output; ADIT does not yet parse it.")])


class GrrmGenerator(InputGenerator):
    code = "grrm"
    uses_kpoints = False
    supports_analysis = False
    cli_only = True

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, GrrmMethod):
            return [ValidationError("method.code", L("GRRM17 の条件ではありません", "not a GRRM17 method"))]
        errors = _molecule_errors(spec, tasks=("geometry_optimization", "vibrations"), code="GRRM17")
        if spec.task.type == "geometry_optimization" and spec.task.relax_cell != "no":
            errors.append(ValidationError("task.relax_cell", L("GRRM17 のこの生成器はセルを緩和しません", "this GRRM17 generator cannot relax a cell")))
        errors += _word_error(m.theory, "method.theory", L("手法", "method"))
        errors += _word_error(m.basis, "method.basis", L("基底関数", "basis set"))
        if spec.runtime.mpiprocs > 1:
            errors.append(ValidationError("runtime.mpiprocs", L(
                "この GRRM17 生成器は逐次実行だけに対応します", "this GRRM17 generator supports serial runs only")))
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st = spec.method, spec.structure
        job = "MIN" if spec.task.type == "geometry_optimization" else "FREQ"
        lines = [f"# {job}/{m.theory}/{m.basis}", "", f"{st.charge} {st.multiplicity}"]
        lines += [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        return {"grrm.com": "\n".join(lines + ["", ""])}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = _format_command(spec, profile, self.code, "GRRM17p")
        return f"{exe} grrm > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(program="GRRM17p", files=[L("  grrm.com     GRRM17 の MIN または FREQ 入力", "  grrm.com     GRRM17 MIN or FREQ input")],
                           prepare=[L("  GRRM17 と連携する Gaussian を各自で用意し、公式マニュアルに従って subgrr・subgau・subchk を設定してください。",
                                      "  Provide GRRM17 and a licensed Gaussian installation; configure subgrr, subgau and subchk as documented by GRRM17."),
                                    L("  共通 Spec の最適化器・最大反復数・力の閾値は grrm.com に適用していません。",
                                      "  The common optimizer, maximum step count, and force threshold are not applied to grrm.com.")],
                           outputs=[L("  grrm.log     GRRM17 の計算結果。ADIT はまだ内容を解析しません。",
                                      "  grrm.log     GRRM17 results; ADIT does not yet parse them.")])


register(GaussianGenerator())
register(GamessGenerator())
register(QchemGenerator())
register(GrrmGenerator())
