"""NWChem molecular single-point input (no implicit chemistry choices)."""
# Official syntax and scope:
# https://nwchemgit.github.io/Getting-Started.html
# https://nwchemgit.github.io/Geometry.html
# https://nwchemgit.github.io/Basis.html
# https://nwchemgit.github.io/Charge.html
# https://nwchemgit.github.io/Hartree-Fock-Theory-for-Molecules.html
# https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html
# https://nwchemgit.github.io/TASK.html

from __future__ import annotations

import re
from pathlib import Path

from ase.data import atomic_numbers

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, NwchemMethod
from adit.validate import electron_parity_error
from adit.validate_types import ValidationError


_TOKEN = re.compile(r"[A-Za-z0-9_+*().-]+\Z")


class NwchemGenerator(InputGenerator):
    code = "nwchem"
    uses_kpoints = False

    def resolve(self, spec: CalculationSpec, cfg: Config) -> None:
        return None  # The basis is looked up in the user's NWChem installation.

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m, st, task = spec.method, spec.structure, spec.task
        if not isinstance(m, NwchemMethod):
            return [ValidationError("method.code", L("NWChem の条件ではありません", "the settings are not for NWChem"))]
        errors: list[ValidationError] = []
        if st.periodic or any(st.atoms.pbc):
            errors.append(ValidationError("structure.atoms", L("この NWChem 生成器は分子の計算だけに対応します", "this NWChem generator supports molecular calculations only")))
        if st.fixed_atoms or st.fixed_axes:
            errors.append(ValidationError("structure.fixed_atoms", L("この NWChem 生成器では固定原子・固定軸を入力に書けません", "this NWChem generator cannot write fixed atoms or axes")))
        if spec.handoff is not None and spec.handoff.files:
            errors.append(ValidationError("handoff.files", L("この NWChem 生成器は前の計算のファイルを引き継げません", "this NWChem generator cannot carry over restart files")))
        if task.type != "single_point":
            errors.append(ValidationError("task.type", L("この NWChem 生成器は分子の一点計算だけに対応します", "this NWChem generator supports molecular single-point calculations only")))
        if ((spec.runtime.mpiprocs > 1 or spec.runtime.omp_threads > 1)
                and spec.runtime.profile in cfg.profiles
                and not cfg.profile(spec.runtime.profile).commands.get(self.code)):
            errors.append(ValidationError("runtime.profile", L("並列実行では環境設定の commands.nwchem に実行コマンドを明示してください", "for parallel execution, explicitly set commands.nwchem in your settings")))
        if m.theory is None:
            errors.append(ValidationError("method.theory", L("SCF または DFT を明示してください", "explicitly select SCF or DFT")))
        if not _TOKEN.fullmatch(m.basis):
            errors.append(ValidationError("method.basis", L("NWChem の basis library にある基底関数名を、空白や改行のない 1 語で指定してください", "provide a one-token basis name from your NWChem basis library, without spaces or line breaks")))
        if m.theory == "dft" and not _TOKEN.fullmatch(m.xc):
            errors.append(ValidationError("method.xc", L("DFT の XC キーワードを空白や改行のない 1 語で明示してください", "explicitly provide a one-token DFT XC keyword without spaces or line breaks")))
        if m.theory == "scf" and m.xc:
            errors.append(ValidationError("method.xc", L("SCF では XC キーワードを使用しません。空にしてください", "SCF does not use an XC keyword; leave it empty")))
        if m.theory == "scf" and st.multiplicity != 1:
            errors.append(ValidationError("structure.multiplicity", L("この生成器の SCF は閉殻一重項だけです。開殻では DFT を選ぶか、別の入力を用意してください", "this generator's SCF subset is closed-shell singlet only; select DFT for an open-shell system or prepare a separate input")))
        if m.theory is not None:
            n_electrons = sum(atomic_numbers[symbol] for symbol in st.atoms.symbols) - st.charge
            parity_error = electron_parity_error(n_electrons, st.charge, st.multiplicity)
            if parity_error:
                errors.append(parity_error)
        return errors

    def generate(self, spec: CalculationSpec, res: None) -> dict[str, str]:
        m, st = spec.method, spec.structure
        lines = ["title \"ADIT calculation\"", "geometry units angstroms noautosym nocenter noautoz"]
        lines += [f"  {symbol} {x:.10f} {y:.10f} {z:.10f}"
                  for symbol, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        lines += ["end", f"charge {st.charge}", "basis", f"  * library {m.basis}", "end"]
        if m.theory == "dft":
            lines += ["dft", f"  xc {m.xc}", f"  mult {st.multiplicity}", "end"]
        else:
            lines += ["scf", "  singlet", "end"]
        lines += [f"task {m.theory} energy", ""]
        return {"nwchem.nw": "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res: None) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, "nwchem").format(mpiprocs=spec.runtime.mpiprocs,
                                                               omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} nwchem.nw > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res: None, copies: dict[str, Path]) -> ReadmeNotes:
        return ReadmeNotes(
            program="nwchem",
            files=[L("  nwchem.nw    NWChem の分子一点計算の入力", "  nwchem.nw    NWChem molecular single-point input")],
            prepare=[L("  指定した基底関数と XC キーワードが、ご自身の NWChem で使えるか確認してください。ADIT は基底関数の内容を用意・照合しません。",
                       "  Check that your NWChem installation provides the selected basis and XC keyword. ADIT does not supply or verify basis-set contents."),
                     L("  並列実行では環境設定の commands.nwchem に、ご自身の NWChem に合う起動コマンドを指定します。",
                       "  For parallel runs, set commands.nwchem in your settings to the launch command appropriate for your NWChem installation.")],
            outputs=[L("  output.log   NWChem の標準出力。Total SCF/DFT energy は Hartree 単位です。ADIT はこの値だけを読みます。",
                       "  output.log   NWChem standard output. Total SCF/DFT energy is in hartree; ADIT reads only this value.")],
        )


register(NwchemGenerator())
