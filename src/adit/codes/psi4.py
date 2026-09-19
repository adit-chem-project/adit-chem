
from __future__ import annotations

import re
from pathlib import Path

from ase.data import atomic_numbers

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, Psi4Method
from adit.validate import electron_parity_error
from adit.validate_types import ValidationError

INPUT_FILE = "input.dat"
RESULTS_FILE = "results.json"
FINAL_XYZ = "final.xyz"
DEFAULT_COMMAND = "psi4"
DOC = "https://psicode.org/psi4manual/master/index.html"
_TOKEN = re.compile(r"[A-Za-z0-9_+*()\[\].,-]+\Z")
_KEY = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")
OPEN_SHELL_REFERENCES = ("uhf", "rohf", "uks")
TASK_CALL = {"single_point": "energy", "geometry_optimization": "optimize", "vibrations": "frequencies"}


class Psi4Generator(InputGenerator):
    code = "psi4"
    uses_kpoints = False

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, Psi4Method):
            return [ValidationError("method.code", L(f"Psi4 の生成器に {m.code!r} の手法が渡されました",
                                                     f"the Psi4 generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        st, t = spec.structure, spec.task
        if st.periodic or any(st.atoms.pbc):
            errs.append(ValidationError("structure.atoms", L(
                "この Psi4 生成器は分子 (非周期) だけに対応します", "this Psi4 generator supports molecules (nonperiodic) only")))
        if st.fixed_atoms or st.fixed_axes:
            errs.append(ValidationError("structure.fixed_atoms", L(
                "この Psi4 生成器は固定原子・固定軸を入力に書けません (optking の凍結座標は確かめていません)",
                "this Psi4 generator cannot write fixed atoms or axes (optking's frozen coordinates were not verified)")))
        if st.velocities is not None:
            errs.append(ValidationError("structure.velocities", L(
                "この Psi4 生成器は初速度を使いません (MD に対応していません)",
                "this Psi4 generator does not use initial velocities (MD is not supported)")))
        if spec.kpoints is not None:
            errs.append(ValidationError("kpoints", L(
                "Psi4 は k 点を使いません。指定を外してください", "Psi4 does not use k-points; remove this setting")))
        if spec.handoff is not None:
            errs.append(ValidationError("handoff", L(
                "この Psi4 生成器は続きの計算に対応しません", "this Psi4 generator does not support restarts")))
        if t.type not in TASK_CALL:
            errs.append(ValidationError("task.type", L(
                "この Psi4 生成器は一点計算・構造最適化・振動解析だけです",
                "this Psi4 generator supports single point, geometry optimization and vibrations only")))
        if t.relax_cell != "no":
            errs.append(ValidationError("task.relax_cell", L(
                "分子の計算にセルの緩和はありません", "a molecular calculation has no cell to relax")))
        if not m.method.strip():
            errs.append(ValidationError("method.method", L(
                "手法の名前 (scf、b3lyp、mp2 など) を指定してください。ADIT は手法を選びません",
                "give the method name (scf, b3lyp, mp2, ...); ADIT does not choose it")))
        elif not _TOKEN.fullmatch(m.method.strip()):
            errs.append(ValidationError("method.method", L(
                "手法の名前は空白や改行のない 1 語で書いてください", "the method name must be a single token without spaces or line breaks")))
        if not m.basis.strip():
            errs.append(ValidationError("method.basis", L(
                "基底関数の名前 (cc-pVDZ、6-31G* など) を指定してください", "give the basis-set name (cc-pVDZ, 6-31G*, ...)")))
        elif not _TOKEN.fullmatch(m.basis.strip()):
            errs.append(ValidationError("method.basis", L(
                "基底関数の名前は空白や改行のない 1 語で書いてください", "the basis name must be a single token without spaces or line breaks")))
        if st.multiplicity != 1 and m.reference not in OPEN_SHELL_REFERENCES:
            errs.append(ValidationError("method.reference", L(
                f"多重度が {st.multiplicity} なので、開殻を扱う参照波動関数 ({' / '.join(OPEN_SHELL_REFERENCES)}) を選んでください。Psi4 の既定 (RHF) は閉殻だけです",
                f"the multiplicity is {st.multiplicity}, so choose an open-shell reference ({' / '.join(OPEN_SHELL_REFERENCES)}); the Psi4 default (RHF) is closed-shell only")))
        if m.memory_mb < 0:
            errs.append(ValidationError("method.memory_mb", L("メモリは 0 以上にしてください", "the memory must be non-negative")))
        bad = [k for k in m.extra_set if not _KEY.fullmatch(k.strip())]
        if bad:
            errs.append(ValidationError("method.extra_set", L(
                f"Psi4 の設定の名前として読めません: {bad}", f"not valid Psi4 option names: {bad}")))
        clash = sorted({k.strip().lower() for k in m.extra_set} & {"basis", "reference"})
        if clash:
            errs.append(ValidationError("method.extra_set", L(
                f"この生成器が書く設定と重なっています: {', '.join(clash)}。専用の欄で指定してください",
                f"these clash with options this generator already writes: {', '.join(clash)}; use their own fields")))
        if any("\n" in str(v) or "\r" in str(v) for v in m.extra_set.values()):
            errs.append(ValidationError("method.extra_set", L(
                "設定の値は 1 行で書いてください", "each option value must be a single line")))
        n_electrons = sum(atomic_numbers[symbol] for symbol in st.atoms.symbols) - st.charge
        parity = electron_parity_error(n_electrons, st.charge, st.multiplicity)
        if parity:
            errs.append(parity)
        if t.type == "geometry_optimization" and t.max_steps < 1:
            errs.append(ValidationError("task.max_steps", L(
                "構造最適化の最大回数 (geom_maxiter) は 1 以上にしてください",
                "the maximum number of optimization steps (geom_maxiter) must be at least 1")))
        if ((spec.runtime.mpiprocs > 1) and spec.runtime.profile in cfg.profiles
                and not cfg.profile(spec.runtime.profile).commands.get(self.code)):
            errs.append(ValidationError("runtime.mpiprocs", L(
                "Psi4 はスレッド並列です (MPI では動きません)。プロセス数は 1 にし、スレッド数を増やしてください",
                "Psi4 is threaded, not MPI-parallel; set the number of processes to 1 and raise the thread count")))
        return errs

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        m, st, t = spec.method, spec.structure, spec.task
        method, basis = m.method.strip(), m.basis.strip()
        lines = ["# ADIT が生成した Psi4 の入力 / Psi4 input generated by ADIT",
                 "# 手法・基底・参照波動関数は利用者が指定したものです / the method, basis and reference come from the user", ""]
        if m.memory_mb:
            lines += [f"memory {m.memory_mb} MB", ""]
        lines += ["molecule adit {", f"  {st.charge} {st.multiplicity}"]
        lines += [f"  {symbol} {x:.10f} {y:.10f} {z:.10f}"
                  for symbol, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions)]
        lines += ["  units angstrom", "  no_reorient", "  no_com", "  symmetry c1", "}", "",
                  f"set basis {basis}"]
        if m.reference:
            lines.append(f"set reference {m.reference}")
        if t.type == "geometry_optimization":
            lines.append(f"set geom_maxiter {t.max_steps}")
        lines += [f"set {k.strip()} {v}" for k, v in m.extra_set.items()]
        call = TASK_CALL.get(t.type, "energy")
        lines += ["", "import json", "",
                  '# 下の 2 行の目印は、Psi4 が出力の先頭へ写す入力の行と重ならないよう、文字列を分けてあります',
                  'psi4.core.print_out("adit-" "psi4: psi4 " + psi4.__version__ + "\\n")',
                  f"energy_hartree, wfn = {call}({method!r}, return_wfn=True)",
                  "results = {\"task\": %r, \"method\": %r, \"basis\": %r, \"psi4_version\": psi4.__version__,"
                  % (t.type, method, basis),
                  "           \"energy_hartree\": float(energy_hartree)}"]
        lines += ["mol = psi4.core.get_active_molecule()", f"mol.save_xyz_file({FINAL_XYZ!r}, True)"]
        if t.type == "geometry_optimization":
            lines += ["results[\"geom_maxiter\"] = %d" % t.max_steps, 'results["converged"] = True']
        if t.type == "vibrations":
            lines += ["results[\"frequencies_cm1\"] = [float(x) for x in wfn.frequencies().to_array()]"]
        lines += ['psi4.core.print_out("adit-" "psi4: psi4 " + psi4.__version__ + "\\n")',
                  f"with open({RESULTS_FILE!r}, \"w\") as handle:",
                  "    json.dump(results, handle, indent=2)",
                  'psi4.core.print_out("adit-" "psi4: done " + json.dumps(results) + "\\n")', ""]
        return {INPUT_FILE: "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        r = spec.runtime
        exe = profile.command_for(self.code, DEFAULT_COMMAND).format(
            mpiprocs=r.mpiprocs, omp_threads=r.omp_threads, binary="")
        return f"{exe} -i {INPUT_FILE} -o output.log -n {r.omp_threads} > stdout.log 2>&1"

    def version_probe(self, spec):
        return ("output.log", "adit-psi4: psi4")

    # ---- README ----
    def readme_notes(self, spec: CalculationSpec, res, copies: dict[str, Path]) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        files = [L("  input.dat    Psi4 の入力 (psithon。分子・基底・手法と、結果を results.json に書く数行)",
                   "  input.dat    Psi4 input (psithon: the molecule, basis, method, and a few lines that write results.json)")]
        prepare = [
            L(f"  手法 {m.method} と基底 {m.basis} は利用者の指定です。ADIT は選びません。お使いの Psi4 にその基底があるか確かめてください。",
              f"  The method {m.method} and basis {m.basis} are your choices; ADIT does not select them. Check that your Psi4 provides the basis."),
            L("  Psi4 はスレッド並列です (MPI では分散しません)。submit.sh は -n にスレッド数を渡しています。",
              "  Psi4 is threaded (it does not distribute over MPI); submit.sh passes the thread count to -n."),
        ]
        if not m.reference:
            prepare.append(L("  参照波動関数 (reference) を書いていないので、Psi4 の既定 (RHF / RKS) になります。",
                             "  No reference is written, so the Psi4 default (RHF / RKS) applies."))
        if m.memory_mb:
            prepare.append(L(f"  memory {m.memory_mb} MB を書いています。足りないと Psi4 が止まります。",
                             f"  memory {m.memory_mb} MB is written; Psi4 stops if that is not enough."))
        outputs = [L("  final.xyz     計算に使った最後の構造 (最適化ならその結果)",
                     "  final.xyz     the final structure used (the optimized one for an optimization)"),
                   L("  output.log    Psi4 の出力。最後に「adit-psi4: done」と結果の要約が入ります",
                     "  output.log    Psi4 output; it ends with a 'adit-psi4: done' line and the result summary"),
                   L("  results.json  入力の最後の数行が書く結果 (エネルギー [Hartree]、Psi4 のバージョン)",
                     "  results.json  results written by the last lines of the input (energy in hartree, Psi4 version)"),
                   L("  stdout.log    画面に出る分 (optking の途中経過など)。Psi4 が止まったときの理由もここに出ます",
                     "  stdout.log    what Psi4 prints to the screen (optking progress); the reason for a crash appears here too")]
        if t == "vibrations":
            outputs.append(L("  results.json の frequencies_cm1   振動数 [cm⁻¹] (虚数は負の数で入ります)",
                             "  frequencies_cm1 in results.json   vibrational frequencies in cm^-1 (imaginary ones appear as negative)"))
        return ReadmeNotes(program=DEFAULT_COMMAND, files=files, prepare=prepare, outputs=outputs)


register(Psi4Generator())
