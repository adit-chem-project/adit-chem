"""DCDFTBMD 2.0 input generator (documented subset)."""
# Input sections, .spl order, geometry, OPT, MD, and execution:
# https://www.chem.waseda.ac.jp/dcdftbmd/document/DCDFTBMD_2.0_en.pdf
# The program and its .spl files are supplied by the user; .skf is not converted.

from __future__ import annotations

from pathlib import Path

from adit.codes.base import InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.lang import L
from adit.spec import CalculationSpec, DcdftbmdMethod, HARTREE_PER_BOHR_IN_EV_PER_ANG
from adit.validate_types import ValidationError


class DcdftbmdGenerator(InputGenerator):
    code = "dcdftbmd"
    uses_kpoints = False  # DCDFTBMD 2.0 manual uses TV vectors, not a k-point mesh
    supports_analysis = False  # no reader for dftb.out / traject yet

    def resolve(self, spec: CalculationSpec, cfg: Config) -> dict[str, Path]:
        m = spec.method
        return {pair: Path(path).expanduser() for pair, path in m.sk_files.items()}

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m, st, t = spec.method, spec.structure, spec.task
        if not isinstance(m, DcdftbmdMethod):
            return [ValidationError("method.code", L("DCDFTBMD の条件ではありません", "not a DCDFTBMD method"))]
        errs = []
        if m.scc is None or m.divide_and_conquer is None:
            errs.append(ValidationError("method.scc", L("SCC と分割統治法 (divide_and_conquer) を明示してください", "set both SCC and divide_and_conquer explicitly")))
        symbols = list(dict.fromkeys(st.atoms.symbols))
        expected_pairs = {f"{sym}-{other}" for sym in symbols for other in symbols}
        extra_pairs = sorted(set(m.sk_files) - expected_pairs)
        if extra_pairs:
            errs.append(ValidationError("method.sk_files", L(f"構造にない元素対の指定を外してください: {extra_pairs}", f"remove pairs not present in the structure: {extra_pairs}")))
        if st.multiplicity != 1:
            errs.append(ValidationError("structure.multiplicity", L("開殻の追加パラメータはこの生成器で書けません。多重度 1 の系だけに対応します", "this generator cannot write the additional open-shell parameters; only multiplicity 1 is supported")))
        if ((spec.runtime.mpiprocs > 1 or spec.runtime.omp_threads > 1)
                and spec.runtime.profile in cfg.profiles
                and not cfg.profile(spec.runtime.profile).commands.get(self.code)):
            errs.append(ValidationError("runtime.profile", L("並列実行では環境設定の commands.dcdftbmd に対応する実行コマンドを明示してください", "for parallel execution, set the appropriate commands.dcdftbmd in your settings")))
        for sym in symbols:
            if m.highest_angular_momentum.get(sym) not in (1, 2, 3, 4):
                errs.append(ValidationError("method.highest_angular_momentum", L(f"{sym} の最高角運動量を 1 (s)〜4 (f) で指定してください", f"set the highest angular momentum for {sym} from 1 (s) to 4 (f)")))
            for other in symbols:
                pair = f"{sym}-{other}"
                path = m.sk_files.get(pair, "")
                if not path:
                    errs.append(ValidationError("method.sk_files", L(f"{pair} の .spl ファイルを指定してください (.skf は使えません)", f"provide the {pair} .spl file (.skf is not supported)")))
                elif not path.lower().endswith(".spl") or not Path(path).expanduser().is_file():
                    errs.append(ValidationError("method.sk_files", L(f"{pair} の .spl ファイルがありません: {path}", f"the {pair} .spl file does not exist: {path}")))
        if st.fixed_axes:
            errs.append(ValidationError("structure.fixed_axes", L("DCDFTBMD の軸固定はこの生成器ではまだ書けません", "this generator cannot write DCDFTBMD axis constraints")))
        if st.fixed_atoms and t.type == "single_point":
            errs.append(ValidationError("structure.fixed_atoms", L("一点計算では固定原子を使いません。指定を外してください", "fixed atoms are not used in a single-point calculation; remove them")))
        if t.type not in ("single_point", "geometry_optimization", "molecular_dynamics"):
            errs.append(ValidationError("task.type", L("この生成器は一点計算・構造最適化・MD に対応します", "this generator supports single point, geometry optimization, and MD")))
        if t.type == "geometry_optimization":
            if t.optimizer not in ("SteepestDescent", "FIRE"):
                errs.append(ValidationError("task.optimizer", L("この生成器では最急降下法と FIRE のみ、同じ名前の手法として書けます", "this generator can map only SteepestDescent and FIRE to the same-named methods")))
            if t.relax_cell != "no":
                errs.append(ValidationError("task.relax_cell", L("この生成器は格子最適化を書けません", "this generator cannot write lattice optimization")))
        if t.type == "molecular_dynamics":
            if t.md.ensemble == "NPT":
                errs.append(ValidationError("task.md.ensemble", L("この生成器は NPT を書けません", "this generator cannot write NPT")))
            if t.md.ensemble == "NVT" and t.md.thermostat != "berendsen":
                errs.append(ValidationError("task.md.thermostat", L("この生成器で時定数まで書ける NVT 熱浴は Berendsen だけです", "Berendsen is the only NVT thermostat whose coupling time this generator can map")))
        if st.periodic and not all(st.atoms.pbc):
            errs.append(ValidationError("structure.atoms", L("DCDFTBMD の TV ベクトルは 3 方向の周期セルとして書きます。部分周期には対応していません", "DCDFTBMD TV vectors are written for a fully periodic cell; partial periodicity is unsupported")))
        return errs

    def generate(self, spec: CalculationSpec, res: dict[str, Path]) -> dict[str, str]:
        m, st, t = spec.method, spec.structure, spec.task
        keys = [f"SCC={'TRUE' if m.scc else 'FALSE'}", f"DC={'TRUE' if m.divide_and_conquer else 'FALSE'}"]
        if t.type == "geometry_optimization":
            opttype = {"SteepestDescent": 2, "FIRE": 5}[t.optimizer]
            gconv = t.force_tolerance_ev_per_ang / HARTREE_PER_BOHR_IN_EV_PER_ANG
            keys.append(f"OPT=(MAXITER={t.max_steps} GRADCONV={gconv:.10g} OPTTYPE={opttype})")
        if t.type == "molecular_dynamics":
            md = t.md
            fields = [f"NSTEP={md.steps}", f"DELTAT={md.timestep_fs * 1e-15:.10g}", f"PRINT={md.dump_interval}",
                      f"INITTEMP={md.temperature_k:g}"]
            if md.ensemble == "NVT":
                fields += ["NVT=TRUE", "NVTTYPE=4", f"BATHTEMP={md.temperature_k:g}", f"TAUTEMP={md.coupling_time_fs * 1e-15:.10g}"]
            else:
                fields.append("NVE=TRUE")
            keys.append("MD=(" + " ".join(fields) + ")")
        symbols = list(dict.fromkeys(st.atoms.symbols))
        lines = [" ".join(keys), "", "ADIT calculation", "", str(len(symbols))]
        for sym in symbols:
            lines.append(f"{sym} {m.highest_angular_momentum[sym]}")
            lines.append(" ".join(f"params/{sym}-{other}.spl" for other in symbols))
        lines += ["", f"{len(st.atoms.symbols)} {st.charge} {st.multiplicity}"]
        for i, (sym, pos) in enumerate(zip(st.atoms.symbols, st.atoms.positions)):
            line = f"{sym}0 " + " ".join(f"{x:.8f}" for x in pos)
            if i in st.fixed_atoms:
                line += " *"
            lines.append(line)
        if st.periodic:
            lines += ["TV " + " ".join(f"{x:.8f}" for x in vector) for vector in st.atoms.cell]
        lines.append("")
        return {"dftb.inp": "\n".join(lines)}

    def files_to_copy(self, spec: CalculationSpec, res: dict[str, Path]) -> dict[str, Path]:
        return {f"params/{pair}.spl": path for pair, path in res.items()}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, "dftb_serial.00.x").format(mpiprocs=spec.runtime.mpiprocs,
                                                                        omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, res, copies) -> ReadmeNotes:
        return ReadmeNotes(
            program="dftb_serial.00.x",
            files=[L("  dftb.inp     DCDFTBMD 2.0 の入力 (キーワード・パラメータ・座標)",
                     "  dftb.inp     DCDFTBMD 2.0 input (keywords, parameters, coordinates)"),
                   L("  params/      利用者が指定した順序つき元素対の .spl ファイル", "  params/      user-supplied .spl files for ordered element pairs")],
            prepare=[L("  DCDFTBMD 2.0 本体と、それに対応する .spl パラメータは利用者が用意します。DFTB+ の .skf は変換しません。",
                       "  Supply DCDFTBMD 2.0 and compatible .spl parameters yourself. DFTB+ .skf files are not converted."),
                     L("  実行ファイルが PATH にない場合は、環境設定の commands.dcdftbmd に指定してください。",
                       "  If the executable is not on PATH, set commands.dcdftbmd in your settings.")],
            outputs=[L("  dftb.out     標準の計算結果。MD では traject に座標が出ます。ADIT の自動解析は未対応です。",
                       "  dftb.out     standard results; MD coordinates are in traject. ADIT automatic analysis is not yet supported.")],
        )


register(DcdftbmdGenerator())
