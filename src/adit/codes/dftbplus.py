
from __future__ import annotations

import io
from pathlib import Path

from ase.io import write

from adit.bandpath import KPATH_FILE, band_path, dftb_klines, kpath_json
from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.sk_sets import SKSet, SKSetError, discover_sets
from adit.config import Config, Profile
from adit.spec import CalculationSpec, DftbMethod
from adit.validate import electron_parity_error
from adit import lang
from adit.lang import L
from adit.validate_types import ValidationError

PARSER_VERSION = 14
SKF_SUBDIR = "skf"
GEOMETRY_FILE = "geometry.gen"
INPUT_FILE = "dftb_in.hsd"
D3_KEYS = ("s6", "s8", "a1", "a2")


def _yn(b: bool) -> str:
    return "Yes" if b else "No"


def _fmt(x: float) -> str:
    return f"{x:g}"


SK_DOWNLOAD_URL = "https://dftb.org/parameters/download.html"


class DftbPlusGenerator(InputGenerator):
    code = "dftbplus"

    def resolve(self, spec: CalculationSpec, cfg: Config) -> SKSet:
        errs = self._check_sk(spec, cfg.sk_root or None)
        if errs:
            raise GenerationError("; ".join(str(e) for e in errs))
        return SKSet.from_dir(Path(cfg.sk_root).expanduser() / spec.method.sk_set)

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        errs: list[ValidationError] = []
        m = spec.method
        if not isinstance(m, DftbMethod):
            return [ValidationError("method.code", L(f"DFTB+ の生成器に {m.code!r} の手法が渡されました", f"the DFTB+ generator received a {m.code!r} method"))]
        if not m.sk_set.strip():
            errs.append(ValidationError("method.sk_set", L("Slater-Koster セットが選ばれていません", "no Slater-Koster set selected")))
        if m.scc_tolerance <= 0:
            errs.append(ValidationError("method.scc_tolerance", L("0 より大きい値が必要です", "must be greater than 0")))
        if m.max_scc_iterations < 1:
            errs.append(ValidationError("method.max_scc_iterations", L("1 以上が必要です", "must be at least 1")))
        if m.filling_temperature < 0:
            errs.append(ValidationError("method.filling_temperature", L("負の温度は指定できません", "negative temperature is not allowed")))
        from ase.data import atomic_numbers
        n_el = sum(atomic_numbers[s] for s in spec.structure.atoms.symbols) - spec.structure.charge
        pe = electron_parity_error(n_el, spec.structure.charge, spec.structure.multiplicity)
        if pe:
            errs.append(pe)
        if spec.task.type == "molecular_dynamics" and spec.task.md.ensemble != "NVE" and spec.task.md.thermostat not in ("berendsen", "andersen", "nose_hoover"):
            errs.append(ValidationError("task.md.thermostat", L(f"DFTB+ にはない熱浴です: {spec.task.md.thermostat} (berendsen / andersen / nose_hoover)", f"thermostat not available in DFTB+: {spec.task.md.thermostat} (berendsen / andersen / nose_hoover)")))
        if m.solvation_param_file.strip():
            p = Path(m.solvation_param_file).expanduser()
            if not p.is_file():
                errs.append(ValidationError("method.solvation_param_file", L(f"溶媒のパラメータファイルがありません: {p}", f"solvation parameter file not found: {p}")))
            elif p.name in (INPUT_FILE, GEOMETRY_FILE, "dftb_pin.hsd", "detailed.out") or p.name.startswith("submit"):
                errs.append(ValidationError("method.solvation_param_file", L(f"ファイル名 {p.name} は生成した名前と重なります。名前を変えてください", f"the file name {p.name} clashes with a generated file; rename it")))
            if spec.structure.periodic:
                errs.append(ValidationError("method.solvation_param_file", L("DFTB+ の GeneralisedBorn は有限の系 (分子) だけです (マニュアル 2.4.19 節)", "DFTB+ GeneralisedBorn is for finite systems only (manual 2.4.19)")))
        if m.seed < 1:
            errs.append(ValidationError("method.seed", L(
                "乱数の種は 1 以上にしてください (0 を書くと DFTB+ は毎回選び直し、MD の軌跡が再現できません)",
                "the random seed must be at least 1 (DFTB+ reads 0 as 'choose a new one each time', which makes an MD trajectory irreproducible)")))
        if m.sk_set.strip():
            errs += self._check_sk(spec, cfg.sk_root or None)
        return errs

    writes_velocities = True

    def version_probe(self, spec):
        return ("output.log", r"DFTB\+ .*(release|version)")

    def _check_sk(self, spec: CalculationSpec, sk_root: str | Path | None) -> list[ValidationError]:
        from adit.config import config_path
        where = L(f"環境設定ファイル {config_path()} の sk_root = \"...\" に書きます", f"write it in sk_root = \"...\" in the settings file {config_path()}")
        how = L(f"デスクトップ版では「Slater-Koster パラメータ」の欄の「フォルダを選ぶ…」で選び、ウェブ版とコマンド行では{where}",
                f"In the desktop app, use \"Choose folder…\" next to Slater-Koster set; in the web version and on the command line, {where}")
        if sk_root is None:
            return [ValidationError("method.sk_set", L(
                f"Slater-Koster パラメータの置き場所 (環境設定の sk_root) がまだ決まっていません。{how} "
                f"(例: sk_root = \"/home/<ユーザー名>/slakos\"。その下に mio-1-1/ などのセットのフォルダを置きます。入手先: {SK_DOWNLOAD_URL})",
                f"The Slater-Koster parameter folder (sk_root in the settings) is not set yet. {how} "
                f"(e.g. sk_root = \"/home/<user>/slakos\", containing set folders such as mio-1-1/; download: {SK_DOWNLOAD_URL})"))]
        sk_root = Path(sk_root).expanduser()
        if not sk_root.is_dir():
            return [ValidationError("method.sk_set", L(f"sk_root に書かれたフォルダがありません: {sk_root}。{how}",
                                                       f"the folder given as sk_root does not exist: {sk_root}. {how}"))]
        sets = discover_sets(sk_root)
        skset = sets.get(spec.method.sk_set)
        if skset is None:
            found = ", ".join(sorted(sets)) or L("(なし)", "(none)")
            return [ValidationError("method.sk_set", L(f"セット {spec.method.sk_set!r} が {sk_root} にありません (見つかったセット: {found})",
                                                       f"set {spec.method.sk_set!r} not found in {sk_root} (sets found: {found})"))]
        errs: list[ValidationError] = []
        missing = skset.missing_pairs(spec.elements)
        if missing:
            absent = [e for e in spec.elements if e not in skset.elements]
            pair_only = [f"{a}-{b}" for a, b in missing if a not in absent and b not in absent]
            alt = self._alternative_sets(sets, skset.name, spec.elements, periodic=spec.structure.periodic)
            if absent:
                errs.append(ValidationError("method.sk_set", L(f"セット {skset.name} に無い元素があります: {absent}。{alt}", f"elements missing from set {skset.name}: {absent}. {alt}")))
            if pair_only:
                errs.append(ValidationError("method.sk_set", L(f"セット {skset.name} に無い元素ペアがあります: {pair_only}。{alt}", f"element pairs missing from set {skset.name}: {pair_only}. {alt}")))
        else:
            for e in spec.elements:
                try:
                    skset.max_angular_momentum(e)
                except SKSetError as ex:
                    errs.append(ValidationError("method.sk_set", str(ex)))
            errs += self._check_set_data(spec, skset)
        for doc in ("LICENSE", "README"):
            if doc not in skset.doc_files():
                errs.append(ValidationError("method.sk_set", L(f"セット {skset.name} に {doc} がありません (複製に必須です)", f"set {skset.name} has no {doc} (required for copying)")))
        return errs

    @staticmethod
    def _alternative_sets(sets: dict, current: str, elements: list[str], periodic: bool = False) -> str:
        ok = sorted(n for n, s in sets.items() if n != current and not s.missing_pairs(elements))
        if ok:
            return L(f"これらのセットなら全元素 (と元素の組) がそろっています: {', '.join(ok)}",
                     f"these sets have all the elements (and pairs): {', '.join(ok)}")
        if periodic:
            return L(f"手元のセットにはそろっていないので、その元素を含むセットを入手してください ({SK_DOWNLOAD_URL})。"
                     "周期系なので xtb は使えません (この生成器の xtb は分子だけです)",
                     f"none of your sets has them; obtain a set that includes them ({SK_DOWNLOAD_URL}). "
                     "xtb is not an option here because this generator's xtb handles molecules only")
        return L(f"手元のセットにはそろっていないので、その元素を含むセットを入手する ({SK_DOWNLOAD_URL}) か、全元素に対応する xtb (GFN-xTB) を計算コードに選ぶ方法があります",
                 f"none of your sets has them; obtain a set that includes them ({SK_DOWNLOAD_URL}), or choose xtb (GFN-xTB), which covers all elements, as the code")

    @staticmethod
    def _check_set_data(spec: CalculationSpec, skset: SKSet) -> list[ValidationError]:
        errs: list[ValidationError] = []
        m = spec.method
        if spec.structure.multiplicity != 1:
            w = skset.spin_constants()
            if w is None:
                errs.append(ValidationError("structure.multiplicity", L(f"多重度 {spec.structure.multiplicity} (不対電子あり) の計算には、パラメータのセットにスピン定数の表 (spinw.txt) が要りますが、{skset.name} には入っていません",
                                                                        f"multiplicity {spec.structure.multiplicity} (unpaired electrons) needs a spin-constant table (spinw.txt) in the parameter set, but {skset.name} has none")))
            else:
                absent = [e for e in spec.elements if e not in w]
                if absent:
                    errs.append(ValidationError("structure.multiplicity", L(f"スピン定数の表 (spinw.txt) に無い元素があります: {absent}", f"elements missing from the spin-constant table (spinw.txt): {absent}")))
        if m.third_order:
            d = skset.hubbard_derivs()
            if d is None:
                errs.append(ValidationError("method.third_order", L(f"セット {skset.name} の README に Hubbard 微分値の一覧がありません", f"README of set {skset.name} has no list of Hubbard derivatives")))
            else:
                absent = [e for e in spec.elements if e not in d]
                if absent:
                    errs.append(ValidationError("method.third_order", L(f"README の Hubbard 微分値に無い元素があります: {absent}", f"elements missing from the Hubbard derivatives in README: {absent}")))
            if skset.damping_exponent() is None:
                errs.append(ValidationError("method.third_order", L(f"セット {skset.name} の README に H-X 減衰の指数 (zeta) がありません", f"README of set {skset.name} has no H-X damping exponent (zeta)")))
        if m.dispersion == "dftd3":
            absent = [k for k in D3_KEYS if k not in (m.d3_params or {})]
            if absent:
                errs.append(ValidationError("method.d3_params", L(f"DFT-D3 の係数がありません: {absent}", f"DFT-D3 parameters missing: {absent}")))
        return errs

    def generate(self, spec: CalculationSpec, skset: SKSet) -> dict[str, str]:
        out = {GEOMETRY_FILE: self.geometry_gen(spec), INPUT_FILE: self.dftb_in_hsd(spec, skset)}
        if spec.task.type == "band_structure":
            kp = band_path(spec.atoms, spec.task.bands.path, spec.task.bands.npoints)
            out["bands/" + INPUT_FILE] = self.dftb_in_hsd(spec, skset, bands=kp)
            out["bands/" + KPATH_FILE] = kpath_json(kp, spec.atoms)
        return out

    def files_to_copy(self, spec: CalculationSpec, skset: SKSet) -> dict[str, Path]:
        docs = skset.doc_files()
        missing = [d for d in ("LICENSE", "README") if d not in docs]
        if missing:
            raise GenerationError(L(f"セット {skset.name} に {missing} がありません。複製に必須なので生成しません", f"set {skset.name} lacks {missing}; they are required for copying, so nothing is generated"))
        out = {f"{SKF_SUBDIR}/{p.name}": p for p in skset.required_files(spec.elements).values()}
        out.update({f"{SKF_SUBDIR}/{n}": p for n, p in docs.items()})
        if spec.method.solvation_param_file.strip():
            p = Path(spec.method.solvation_param_file).expanduser()
            out[p.name] = p
        return out

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        cmd = profile.command_for(self.code, "dftb+").format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads)
        if spec.task.type == "band_structure":
            return f"{cmd} > output.log 2>&1 && cp charges.bin bands/ && cd bands && {cmd} > output.log 2>&1"
        return f"{cmd} > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, skset: SKSet, copies: dict[str, Path]) -> ReadmeNotes:
        t = spec.task.type
        skf_names = sorted(n.split("/", 1)[1] for n in copies if n.startswith(SKF_SUBDIR + "/"))
        files = [
            L("  dftb_in.hsd   DFTB+ の入力 (計算の設定。DFTB+ はこのディレクトリで起動すると自動でこれを読みます)",
              "  dftb_in.hsd   DFTB+ input (the settings; DFTB+ reads it automatically when started in this directory)"),
            L("  geometry.gen  構造 (原子の種類と座標。gen 形式、長さの単位は Å)",
              "  geometry.gen  structure (elements and coordinates; gen format, lengths in Å)"),
            L(f"  skf/          Slater-Koster ファイル (セット {skset.name}。DFTB+ が使う、元素の組ごとのパラメータ) と LICENSE、README: {', '.join(skf_names)}",
              f"  skf/          Slater-Koster files (set {skset.name}; DFTB+ parameters for each pair of elements) plus LICENSE and README: {', '.join(skf_names)}"),
        ]
        if t == "band_structure":
            files.append(L("  bands/        バンド計算の 2 段階目。1 段階目の電荷を読み、高対称点を結ぶ経路 (bands/kpath.json) に沿って計算します",
                           "  bands/        second stage of the band calculation: reads the charges of the first stage and follows the high-symmetry path (bands/kpath.json)"))
        out = [
            L("  output.log    実行ログ (SCC = 電荷を自己無撞着に決める反復。その様子と、各ステップのエネルギー。単位は Hartree)",
              "  output.log    run log (SCC = the self-consistent charge iterations, and the energy of each step, in Hartree)"),
            L("  detailed.out  最後のステップのエネルギーの内訳、Mulliken 電荷 (原子ごとの電荷の目安)、原子にかかる力",
              "  detailed.out  energy breakdown of the last step, Mulliken charges (approximate charge on each atom), forces"),
            L("  results.tag   全エネルギーなどを、プログラムで読みやすい形で書いたもの",
              "  results.tag   total energy etc. in a form that programs can read"),
        ]
        out += {
            "geometry_optimization": [L("  geom.out.gen / geom.out.xyz   最適化後の構造 (.xyz は分子ビューアで開けます)",
                                        "  geom.out.gen / geom.out.xyz   optimized structure (.xyz opens in molecular viewers)")],
            "molecular_dynamics": [L("  geo_end.xyz   MD の軌跡 (MDRestartFrequency ステップごとの構造と速度)",
                                     "  geo_end.xyz   MD trajectory (structure and velocities every MDRestartFrequency steps)"),
                                   L("  md.out        各ステップのエネルギーと温度", "  md.out        energy and temperature of each step")],
            "vibrations": [L("  hessian.out   ヘシアン (エネルギーの 2 階微分)。振動数は ADIT の解析タブか analyze.py で求めます",
                             "  hessian.out   Hessian (second derivatives of the energy); frequencies are computed by the ADIT analysis tab or analyze.py")],
            "band_structure": [L("  bands/band.out  経路上の各 k 点のエネルギー準位 (図は ADIT の解析タブか analyze.py で描けます)",
                                 "  bands/band.out  energy levels at each k-point on the path (plot them with the ADIT analysis tab or analyze.py)")],
        }.get(t, [])
        out.append(L("  dftb_pin.hsd  省略した設定を既定値で埋めた入力。実際に使われた設定を確かめられます",
                     "  dftb_pin.hsd  the input with omitted settings filled with defaults; shows what was actually used"))
        return ReadmeNotes(program="dftb+", files=files, outputs=out)

    def geometry_gen(self, spec: CalculationSpec) -> str:
        buf = io.StringIO()
        write(buf, spec.atoms, format="gen")
        return buf.getvalue()

    def dftb_in_hsd(self, spec: CalculationSpec, skset: SKSet, bands=None) -> str:
        blocks = [
            self._geometry_block("../" + GEOMETRY_FILE if bands else GEOMETRY_FILE),
            "Driver {}\n" if bands else self._driver_block(spec),
            self._hamiltonian_block(spec, skset, bands=bands),
            f"Options {{\n  RandomSeed = {spec.method.seed}\n  WriteResultsTag = Yes\n}}\n",
            "Analysis {\n  PrintForces = Yes\n}\n",
            f"ParserOptions {{\n  ParserVersion = {PARSER_VERSION}\n}}\n",
        ]
        return "\n".join(blocks)

    def _geometry_block(self, geometry_file: str = GEOMETRY_FILE) -> str:
        return f'Geometry = GenFormat {{\n  <<< "{geometry_file}"\n}}\n'

    def _driver_block(self, spec: CalculationSpec) -> str:
        t = spec.task
        if t.type in ("single_point", "band_structure"):
            return "Driver {}\n"
        if spec.structure.fixed_axes:
            raise GenerationError(L("DFTB+ の Driver は軸ごとの固定 (fixed_axes) に対応していません。原子ごとの固定 (fixed_atoms) にしてください", "the DFTB+ driver cannot fix individual axes (fixed_axes); fix whole atoms (fixed_atoms) instead"))
        if t.type == "molecular_dynamics":
            return self._md_block(spec)
        if t.type == "vibrations":
            return "Driver = SecondDerivatives {\n" + f"  Atoms = {self._moved_atoms(spec)}\n" + "}\n"
        lines = [
            "Driver = GeometryOptimization {",
            f"  Optimizer = {t.optimizer} {{}}",
            f"  MovedAtoms = {self._moved_atoms(spec)}",
        ]
        if t.relax_cell != "no":
            lines.append("  LatticeOpt = Yes")
            if t.relax_cell == "volume_only":
                lines.append("  Isotropic = Yes")
        lines += [
            f"  MaxSteps = {t.max_steps}",
            '  OutputPrefix = "geom.out"',
            f"  Convergence {{ GradElem [eV/AA] = {_fmt(t.force_tolerance_ev_per_ang)} }}",
            "}",
        ]
        return "\n".join(lines) + "\n"

    def _md_block(self, spec: CalculationSpec) -> str:
        md = spec.task.md
        lines = ["Driver = VelocityVerlet {", f"  MovedAtoms = {self._moved_atoms(spec)}", f"  Steps = {md.steps}",
                 f"  TimeStep [fs] = {_fmt(md.timestep_fs)}", f"  MDRestartFrequency = {md.dump_interval}", '  OutputPrefix = "geo_end"']
        vel = spec.structure.velocities
        given = vel is not None or (spec.handoff is not None and spec.handoff.at_run and spec.handoff.velocities)
        if vel is not None:
            lines += ["  Velocities [AA/ps] {", *[f"    {v[0] * 1000.0:.10f} {v[1] * 1000.0:.10f} {v[2] * 1000.0:.10f}" for v in vel], "  }"]
        elif given:
            from adit.handoff import VELOCITY_FILE
            lines += ["  Velocities [AA/ps] {", f'    <<< "{VELOCITY_FILE}"', "  }"]
        T = f"Temperature [Kelvin] = {_fmt(md.temperature_k)}"
        if md.ensemble == "NVE" and given:
            lines += ["  Thermostat = None {}"]
        elif md.ensemble == "NVE":
            lines += ["  Thermostat = None {", f"    InitialTemperature [Kelvin] = {_fmt(md.temperature_k)}", "  }"]
        elif md.thermostat == "berendsen":
            lines += ["  Thermostat = Berendsen {", f"    {T}", f"    Timescale [fs] = {_fmt(md.coupling_time_fs)}", "  }"]
        elif md.thermostat == "andersen":
            prob = min(1.0, md.timestep_fs / md.coupling_time_fs)
            lines += ["  Thermostat = Andersen {", f"    {T}", f"    ReselectProbability = {_fmt(prob)}", "    ReselectIndividually = Yes", "  }"]
        elif md.thermostat == "nose_hoover":
            lines += ["  Thermostat = NoseHoover {", f"    {T}", f"    CouplingStrength [THz] = {_fmt(1000.0 / md.coupling_time_fs)}", "  }"]
        else:
            raise GenerationError(L(f"DFTB+ にはない熱浴です: {md.thermostat} (使えるのは berendsen / andersen / nose_hoover)", f"thermostat not available in DFTB+: {md.thermostat} (use berendsen / andersen / nose_hoover)"))
        if md.ensemble == "NPT":
            lines += ["  Barostat {", f"    Pressure [Pa] = {_fmt(md.pressure_bar * 1e5)}",
                      f"    Timescale [fs] = {_fmt(md.barostat_time_fs)}", "    Isotropic = Yes", "  }"]
        lines.append("}")
        return "\n".join(lines) + "\n"

    @staticmethod
    def _moved_atoms(spec: CalculationSpec) -> str:
        fixed = set(spec.structure.fixed_atoms)
        if not fixed:
            return "1:-1"
        moved = [i + 1 for i in range(len(spec.structure.atoms.symbols)) if i not in fixed]
        return " ".join(str(i) for i in moved) if moved else "{}"

    def _kpoints_lines(self, spec: CalculationSpec) -> list[str]:
        if not spec.structure.periodic:
            return []
        kp = spec.kpoints
        n1, n2, n3 = kp.resolved_mesh(spec.structure.atoms.cell)
        s = kp.shift if kp.mode != "gamma" else (0.0, 0.0, 0.0)
        return [
            "  KPointsAndWeights = SupercellFolding {",
            f"    {n1} 0 0",
            f"    0 {n2} 0",
            f"    0 0 {n3}",
            f"    {_fmt(s[0])} {_fmt(s[1])} {_fmt(s[2])}",
            "  }",
        ]

    def _hamiltonian_block(self, spec: CalculationSpec, skset: SKSet, bands=None) -> str:
        m, st = spec.method, spec.structure
        elems = spec.elements
        lines = ["Hamiltonian = DFTB {", f"  Scc = {_yn(m.scc)}"]
        if m.scc and bands is not None:
            lines += ["  ReadInitialCharges = Yes", "  MaxSccIterations = 1", "  ConvergentSccOnly = No"]
        elif m.scc:
            lines += [f"  SccTolerance = {_fmt(m.scc_tolerance)}", f"  MaxSccIterations = {m.max_scc_iterations}"]
        lines += [
            "  SlaterKosterFiles = Type2FileNames {",
            f'    Prefix = "{"../" if bands else ""}{SKF_SUBDIR}/"',
            '    Separator = "-"',
            '    Suffix = ".skf"',
            "  }",
            "  MaxAngularMomentum {",
            *[f'    {e} = "{skset.max_angular_momentum(e)}"' for e in elems],
            "  }",
        ]
        lines += dftb_klines(bands) if bands is not None else self._kpoints_lines(spec)
        if st.charge != 0:
            lines.append(f"  Charge = {st.charge}")
        if st.multiplicity != 1:
            lines += self._spin_lines(spec, skset)
        if m.filling_temperature > 0:
            lines += ["  Filling = Fermi {", f"    Temperature [K] = {_fmt(m.filling_temperature)}", "  }"]
        if m.third_order:
            lines += self._dftb3_lines(spec, skset)
        lines += self._dispersion_lines(spec)
        if m.solvation_param_file.strip():
            lines += ["  Solvation = GeneralisedBorn {", f'    ParamFile = "{Path(m.solvation_param_file).name}"', "  }"]
        lines.append("}")
        return "\n".join(lines) + "\n"

    def _spin_lines(self, spec: CalculationSpec, skset: SKSet) -> list[str]:
        w = skset.spin_constants()
        if w is None:
            raise GenerationError(L(f"セット {skset.name} に spinw.txt が無く、スピン定数を取れません (多重度 1 以外は生成できません)", f"set {skset.name} has no spinw.txt, so spin constants are unavailable (multiplicity other than 1 cannot be generated)"))
        absent = [e for e in spec.elements if e not in w]
        if absent:
            raise GenerationError(L(f"spinw.txt に無い元素があります: {absent}", f"elements missing from spinw.txt: {absent}"))
        unpaired = spec.structure.multiplicity - 1
        out = ["  SpinPolarization = Colinear {", f"    UnpairedElectrons = {unpaired}", "  }", "  SpinConstants {"]
        out += [f"    {e} = {{ {_fmt(w[e])} }}" for e in spec.elements]
        return out + ["  }"]

    def _dftb3_lines(self, spec: CalculationSpec, skset: SKSet) -> list[str]:
        derivs = skset.hubbard_derivs()
        if derivs is None:
            raise GenerationError(L(f"セット {skset.name} の README に Hubbard 微分値の一覧がありません (DFTB3 は生成できません)", f"README of set {skset.name} has no Hubbard derivatives (DFTB3 cannot be generated)"))
        absent = [e for e in spec.elements if e not in derivs]
        if absent:
            raise GenerationError(L(f"README の Hubbard 微分値に無い元素があります: {absent}", f"elements missing from the Hubbard derivatives in README: {absent}"))
        zeta = skset.damping_exponent()
        if zeta is None:
            raise GenerationError(L(f"セット {skset.name} の README に H-X 減衰の指数 (zeta) がありません (DFTB3 は生成できません)", f"README of set {skset.name} has no H-X damping exponent (DFTB3 cannot be generated)"))
        out = ["  ThirdOrderFull = Yes", "  HubbardDerivs {"]
        out += [f"    {e} = {_fmt(derivs[e])}" for e in spec.elements]
        out += ["  }", "  HCorrection = Damping {", f"    Exponent = {_fmt(zeta)}", "  }"]
        return out

    def _dispersion_lines(self, spec: CalculationSpec) -> list[str]:
        m = spec.method
        if m.dispersion == "none":
            return []
        if m.dispersion == "lennard-jones":
            return ["  Dispersion = LennardJones {", "    Parameters = UFFParameters {}", "  }"]
        p = m.d3_params or {}
        absent = [k for k in D3_KEYS if k not in p]
        if absent:
            raise GenerationError(L(f"DFT-D3 の係数 {absent} がありません (パーサのバージョン 7 以降は既定値がありません)", f"DFT-D3 parameters {absent} are missing (no defaults since parser version 7)"))
        return [
            "  Dispersion = DftD3 {",
            "    Damping = BeckeJohnson {",
            f"      a1 = {_fmt(p['a1'])}",
            f"      a2 = {_fmt(p['a2'])}",
            "    }",
            f"    s6 = {_fmt(p['s6'])}",
            f"    s8 = {_fmt(p['s8'])}",
            "  }",
        ]


register(DftbPlusGenerator())


# References the DFTB+ developers ask for (dftbplus.org names the 2025 paper, the GitHub README the 2020 one)
CITATIONS = (
    Citation("dftbplus_hourahine2025", r"""@article{dftbplus_hourahine2025,
  author  = {Hourahine, B. and Berdakin, M. and Bich, J. A. and Bonaf{\'e}, F. P. and Camacho, C. and Cui, Q. and Deshaye, M. Y. and D{\'\i}az Mir{\'o}n, G. and Ehlert, S. and Elstner, M. and Frauenheim, T. and Goldman, N. and Gonz{\'a}lez Le{\'o}n, R. A. and van der Heide, T. and Irle, S. and Kowalczyk, T. and Kuba{\v{r}}, T. and Lee, I. S. and Lien-Medrano, C. R. and Maryewski, A. and Melson, T. and Min, S. K. and Niehaus, T. and Niklasson, A. M. N. and Pecchia, A. and Reuter, K. and S{\'a}nchez, C. G. and Scheurer, C. and Sentef, M. A. and Stishenko, P. V. and Vuong, V. Q. and Aradi, B.},
  title   = {Recent Developments in {DFTB}+, a Software Package for Efficient Atomistic Quantum Mechanical Simulations},
  journal = {The Journal of Physical Chemistry A},
  volume  = {129},
  number  = {24},
  pages   = {5373--5390},
  year    = {2025},
  doi     = {10.1021/acs.jpca.5c01146}
}""", doi="10.1021/acs.jpca.5c01146", source="https://dftbplus.org/about/index.html"),
    Citation("dftbplus_hourahine2020", r"""@article{dftbplus_hourahine2020,
  author  = {Hourahine, B. and Aradi, B. and Blum, V. and Bonaf{\'e}, F. and Buccheri, A. and Camacho, C. and Cevallos, C. and Deshaye, M. Y. and Dumitric{\u{a}}, T. and Dominguez, A. and Ehlert, S. and Elstner, M. and van der Heide, T. and Hermann, J. and Irle, S. and Kranz, J. J. and K{\"o}hler, C. and Kowalczyk, T. and Kuba{\v{r}}, T. and Lee, I. S. and Lutsker, V. and Maurer, R. J. and Min, S. K. and Mitchell, I. and Negre, C. and Niehaus, T. A. and Niklasson, A. M. N. and Page, A. J. and Pecchia, A. and Penazzi, G. and Persson, M. P. and {\v{R}}ez{\'a}{\v{c}}, J. and S{\'a}nchez, C. G. and Sternberg, M. and St{\"o}hr, M. and Stuckenberg, F. and Tkatchenko, A. and Yu, V. W.-z. and Frauenheim, T.},
  title   = {{DFTB}+, a software package for efficient approximate density functional theory based atomistic simulations},
  journal = {The Journal of Chemical Physics},
  volume  = {152},
  number  = {12},
  pages   = {124101},
  year    = {2020},
  doi     = {10.1063/1.5143190}
}""", doi="10.1063/1.5143190", source="https://github.com/dftbplus/dftbplus#citing"),
)
