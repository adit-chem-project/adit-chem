
from __future__ import annotations

import io
import math
from pathlib import Path

from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.spec import CalculationSpec, OrcaMethod, HARTREE_PER_BOHR_IN_EV_PER_ANG
from adit.validate import electron_parity_error
from adit import lang
from adit.lang import L
from adit.validate_types import ValidationError

INPUT_FILE = "orca.inp"
DEFAULT_COMMAND = "orca"


def _guest_xyz_error(path: Path, *, allow_default_charge: bool) -> bool:
    # Stream a single- or multi-frame GUEST XYZ without changing it.
    from ase.data import atomic_numbers

    try:
        with path.open(encoding="utf-8") as stream:
            frames = 0
            while header := stream.readline():
                count = int(header.strip())
                if count < 1:
                    return True
                comment = stream.readline().split()
                try:
                    if len(comment) == 2:
                        int(comment[0])
                        explicit_charge = int(comment[1]) >= 1
                    else:
                        explicit_charge = False
                except ValueError:
                    explicit_charge = False
                if not explicit_charge and not allow_default_charge:
                    return True
                for _ in range(count):
                    words = stream.readline().split()
                    if len(words) < 4 or words[0] not in atomic_numbers or not all(math.isfinite(float(x)) for x in words[1:4]):
                        return True
                frames += 1
            return frames == 0
    except (OSError, UnicodeError, ValueError, IndexError):
        return True


def solvent_names(model: str) -> list[str]:
    from adit.codes.orca_solvents import SOLVENTS

    col = 1 if model == "cpcm" else 2
    out = []
    for row in SOLVENTS:
        if row[col]:
            out.append(next((a for a in row[0] if " " not in a), row[0][0]))
    return out


def _check_solvent(m: OrcaMethod) -> list[ValidationError]:
    from adit.codes.orca_solvents import SOLVENTS

    name = m.solvent.strip()
    if m.solvation == "none":
        if name:
            return [ValidationError("method.solvent", L(f"溶媒モデルが「なし」なので、溶媒 {name!r} は使われません (モデルを CPCM か SMD にするか、溶媒を空にしてください)",
                                                        f"the solvation model is 'none', so the solvent {name!r} would not be used (choose CPCM or SMD, or clear the solvent)"))]
        return []
    model = m.solvation.upper()
    if not name:
        return [ValidationError("method.solvent", L(f"溶媒の名前を入れてください (ORCA 6.1 マニュアルの溶媒の表の名前。例 water)", "enter the solvent name (from the solvent table of the ORCA 6.1 manual, e.g. water)"))]
    row = next((r for r in SOLVENTS if name.lower() in r[0]), None)
    if row is None:
        return [ValidationError("method.solvent", L(f"溶媒 {name!r} は ORCA 6.1 マニュアルの溶媒の表にありません", f"the solvent {name!r} is not in the solvent table of the ORCA 6.1 manual"))]
    if not row[1 if m.solvation == "cpcm" else 2]:
        return [ValidationError("method.solvent", L(f"溶媒 {name!r} は ORCA 6.1 の表で {model} に使えない印になっています", f"the ORCA 6.1 table marks the solvent {name!r} as not available for {model}"))]
    if " " in name:
        alt = [a for a in row[0] if " " not in a]
        return [ValidationError("method.solvent", L(f"空白を含む名前は ! 行に書けません。同じ溶媒の別名を使ってください: {', '.join(alt) or '(なし)'}",
                                                    f"names with spaces cannot go on the ! line; use another name of the same solvent: {', '.join(alt) or '(none)'}"))]
    return []


class OrcaGenerator(InputGenerator):
    code = "orca"

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, OrcaMethod):
            return [ValidationError("method.code", L(f"ORCA の生成器に {m.code!r} の手法が渡されました", f"the ORCA generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        if spec.structure.periodic:
            errs.append(ValidationError("structure.atoms", L("ORCA は分子系だけです (周期セルは扱いません)", "ORCA is for molecules only (no periodic cells)")))
        if not m.method.strip():
            errs.append(ValidationError("method.method", L("手法のキーワード (例 HF, B3LYP, PBE) が空です", "method keyword (e.g. HF, B3LYP, PBE) is empty")))
        if not m.basis.strip() and not ((m.goat or m.docker_guest_file) and m.method.strip().upper() in {"XTB", "GFN2", "GFN1", "GFN2-XTB", "GFN1-XTB", "GFN0-XTB", "GFN-FF"}):
            errs.append(ValidationError("method.basis", L("基底関数のキーワード (例 def2-SVP) が空です", "basis set keyword (e.g. def2-SVP) is empty")))
        if m.scf_maxiter < 1:
            errs.append(ValidationError("method.scf_maxiter", L("1 以上が必要です", "must be at least 1")))
        if m.maxcore_mb < 0:
            errs.append(ValidationError("method.maxcore_mb", L("負の値は指定できません", "negative values are not allowed")))
        if spec.structure.fixed_axes:
            errs.append(ValidationError("structure.fixed_axes", L("ORCA の生成器は軸ごとの固定に対応していません (原子ごとの固定は Constraints で行います)", "the ORCA generator cannot fix individual axes (whole atoms are fixed via Constraints)")))
        t = spec.task
        if t.type == "molecular_dynamics":
            if t.md.ensemble == "NPT":
                errs.append(ValidationError("task.md.ensemble", L("ORCA の MD に NPT はありません", "ORCA MD has no NPT")))
            elif t.md.ensemble == "NVT" and t.md.thermostat not in ("berendsen", "csvr", "nose_hoover"):
                errs.append(ValidationError("task.md.thermostat", L(f"ORCA にはない熱浴です: {t.md.thermostat} (berendsen / csvr / nose_hoover)", f"thermostat not available in ORCA: {t.md.thermostat} (berendsen / csvr / nose_hoover)")))
        bad = [w for w in (m.method, m.basis, m.extra_keywords) if "\n" in w]
        if bad:
            errs.append(ValidationError("method.extra_keywords", L("キーワードは 1 行で書いてください", "keywords must be a single line")))
        if spec.task.type == "band_structure":
            errs.append(ValidationError("task.type", L("ORCA にバンド計算はありません (周期系のコードで行ってください)", "ORCA has no band structure (use a periodic code)")))
        if spec.structure.fixed_atoms and t.type == "vibrations":
            errs.append(ValidationError("structure.fixed_atoms", L(
                "ORCA の振動解析 (Freq) には固定原子を反映できません (部分ヘシアンは %freq の Partial_Hess を extra_blocks に書いてください)",
                "fixed atoms cannot be applied to an ORCA frequency calculation (for a partial Hessian write %freq Partial_Hess in extra_blocks)")))
        if spec.structure.fixed_atoms and t.type == "single_point" and m.irc:
            errs.append(ValidationError("structure.fixed_atoms", L("ORCA の IRC には固定原子を反映できません", "fixed atoms cannot be applied to an ORCA IRC calculation")))
        errs += _check_solvent(m)
        errs += self._check_ts_irc(spec)
        if m.goat:
            if t.type != "single_point":
                errs.append(ValidationError("method.goat", L("GOAT は計算の種類を一点計算にして使います", "use GOAT with the single-point task")))
            if m.ts_search or m.irc:
                errs.append(ValidationError("method.goat", L("GOAT と遷移状態探索・IRC は別の計算にしてください", "run GOAT and transition-state search or IRC as separate calculations")))
            if spec.structure.fixed_atoms:
                errs.append(ValidationError("structure.fixed_atoms", L("GOAT の入力では固定原子を反映できません", "GOAT input cannot apply fixed atoms")))
        if m.docker_guest_file:
            if t.type != "single_point" or m.goat or m.ts_search or m.irc:
                errs.append(ValidationError("method.docker_guest_file", L("ORCA DOCKER は一点計算として、GOAT・遷移状態探索・IRC とは別に生成してください", "generate ORCA DOCKER as a single-point task, separate from GOAT, transition-state search, and IRC")))
            if m.method.strip().upper() not in {"XTB", "GFN2-XTB", "GFN1-XTB", "GFN0-XTB", "GFN-FF"} or m.basis.strip():
                errs.append(ValidationError("method.method", L("ORCA 6.1 DOCKER は XTB/GFN-xTB/GFN-FF のみです。基底関数欄は空にしてください", "ORCA 6.1 DOCKER supports only XTB/GFN-xTB/GFN-FF; leave the basis field empty")))
            if m.solvation != "none":
                errs.append(ValidationError("method.solvation", L("ORCA DOCKER の溶媒モデルは ALPB のみです。この欄の CPCM/SMD は使えません", "ORCA DOCKER supports only ALPB solvation; CPCM/SMD in this field cannot be used")))
            if spec.structure.fixed_atoms:
                errs.append(ValidationError("structure.fixed_atoms", L("この生成器は DOCKER のホスト固定を反映できません", "this generator cannot apply host constraints in DOCKER")))
            guest = Path(m.docker_guest_file).expanduser()
            if guest.suffix.lower() != ".xyz" or not guest.is_file():
                errs.append(ValidationError("method.docker_guest_file", L("GUEST には既存の .xyz ファイルを指定してください", "GUEST requires an existing .xyz file")))
            elif _guest_xyz_error(guest, allow_default_charge=m.docker_assume_neutral_singlet):
                errs.append(ValidationError("method.docker_guest_file", L("GUEST の各 XYZ フレームには原子数、元素と有限の座標が必要です。2 行目に電荷と多重度の整数を書くか、ORCA の既定 (0, 1) を使うことを指定してください",
                                                                       "each GUEST XYZ frame needs an atom count, elements, and finite coordinates. Put integer charge and multiplicity on the second line, or explicitly choose the ORCA default (0, 1)")))
        from ase.data import atomic_numbers
        n_el = sum(atomic_numbers[s] for s in spec.structure.atoms.symbols) - spec.structure.charge
        pe = electron_parity_error(n_el, spec.structure.charge, spec.structure.multiplicity)
        if pe:
            errs.append(pe)
        return errs

    @staticmethod
    def _check_ts_irc(spec: CalculationSpec) -> list[ValidationError]:
        m, t, errs = spec.method, spec.task, []
        if m.ts_search and t.type != "geometry_optimization":
            errs.append(ValidationError("method.ts_search", L("遷移状態の探索 (OptTS) は、計算の種類を構造最適化にして使います",
                                                              "the transition-state search (OptTS) is used with the geometry-optimization task")))
        if m.irc and t.type != "single_point":
            errs.append(ValidationError("method.irc", L("IRC は、計算の種類を一点計算にして使います (! Freq IRC を書きます)",
                                                        "IRC is used with the single-point task (! Freq IRC is written)")))
        if m.ts_search and m.irc:
            errs.append(ValidationError("method.irc", L("OptTS と IRC を 1 つのジョブにまとめる形は ORCA 6.1 のマニュアルで確かめていないので、段階に分けて生成してください (adit-gen --ts ts+irc)",
                                                        "combining OptTS and IRC in one job was not verified in the ORCA 6.1 manual; generate them as stages (adit-gen --ts ts+irc)")))
        if m.ts_recalc_hess < 0:
            errs.append(ValidationError("method.ts_recalc_hess", L("0 以上にしてください (0 は書かない)", "must be 0 or more (0 = not written)")))
        if m.irc_max_iter < 0:
            errs.append(ValidationError("method.irc_max_iter", L("0 以上にしてください (0 は ORCA の既定)", "must be 0 or more (0 = ORCA default)")))
        return errs

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        return {INPUT_FILE: self.orca_inp(spec)}

    def files_to_copy(self, spec, res) -> dict[str, Path]:
        if spec.method.docker_guest_file:
            return {"guest.xyz": Path(spec.method.docker_guest_file).expanduser()}
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        exe = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        return f"{exe} {INPUT_FILE} > output.log 2>&1"

    def readme_notes(self, spec, res, copies) -> ReadmeNotes:
        t = spec.task.type
        files = [L("  orca.inp      ORCA の入力 (! で始まるキーワード行、%pal / %scf / %geom などの設定、座標)",
                   "  orca.inp      ORCA input (the keyword line starting with !, %pal / %scf / %geom blocks, coordinates)")]
        prep = [
            L("  ORCA 本体は https://orcaforum.kofo.mpg.de/ で利用登録をして入手します (このディレクトリには入っていません)。",
              "  ORCA itself is obtained by registering at https://orcaforum.kofo.mpg.de/ (it is not in this directory)."),
            L("  並列で実行するときも mpirun は付けません (ORCA が自分で並列の部分を起動します)。そのかわり ORCA を",
              "  Do not prepend mpirun even in parallel (ORCA starts its own parallel parts). Instead call ORCA"),
            L("  絶対パス (/ で始まる、省略しない場所の書き方) で呼ぶ必要があり、ADIT の環境設定 (cluster.toml の commands.orca) に書きます。",
              "  by its absolute path (the full location starting with /), written in the ADIT settings (commands.orca in cluster.toml)."),
        ]
        out = [L("  output.log    ORCA の出力 (FINAL SINGLE POINT ENERGY の行が全エネルギー。単位は Eh = ハートリー)",
                 "  output.log    ORCA output (the FINAL SINGLE POINT ENERGY line is the total energy, in Eh = hartree)")]
        out += {
            "geometry_optimization": [L("  orca.xyz      最適化後の構造 (分子ビューアで開けます)。各ステップの構造は orca_trj.xyz",
                                        "  orca.xyz      optimized structure (opens in molecular viewers); orca_trj.xyz has each step")],
            "molecular_dynamics": [L("  trajectory.xyz  MD の軌跡 (dump の間隔ごとの座標)", "  trajectory.xyz  MD trajectory (coordinates every dump interval)")]
            + ([L("                固定原子は %md の Constraint Add Cartesian で止めています (自由度が 1 原子あたり 3 減ります)",
                  "                fixed atoms are held by Constraint Add Cartesian in %md (3 degrees of freedom fewer per atom)")]
               if spec.structure.fixed_atoms else []),
            "vibrations": [L("                振動数は output.log の VIBRATIONAL FREQUENCIES の節 (単位 cm⁻¹)。orca.hess はヘシアン",
                             "                the frequencies are in the VIBRATIONAL FREQUENCIES section of output.log (cm⁻¹); orca.hess holds the Hessian")],
        }.get(t, [])
        m = spec.method
        if m.ts_search and t == "geometry_optimization":
            out.append(L("                遷移状態の探索 (OptTS) です。orca.xyz が探索の終わりの構造、orca.hess が最後のヘシアン"
                         + (" (Freq の振動数は output.log の VIBRATIONAL FREQUENCIES の節。虚振動は負の数で出ます)" if m.ts_freq else ""),
                         "                transition-state search (OptTS); orca.xyz is the final structure and orca.hess the last Hessian"
                         + (" (Freq frequencies are in VIBRATIONAL FREQUENCIES of output.log; imaginary ones are printed as negative)" if m.ts_freq else "")))
        if m.irc and t == "single_point":
            out.append(L("  orca_IRC_Full_trj.xyz   IRC の経路全体 (orca_IRC_F.xyz / orca_IRC_B.xyz が前向き / 後ろ向きの終わりの構造。ORCA 6.1 マニュアルの IRC の節)",
                         "  orca_IRC_Full_trj.xyz   the whole IRC path (orca_IRC_F.xyz / orca_IRC_B.xyz are the forward / backward end structures; IRC section of the ORCA 6.1 manual)"))
        if m.goat and t == "single_point":
            out += [L("  orca.globalminimum.xyz   GOAT が見つけた最低エネルギーの構造", "  orca.globalminimum.xyz   lowest-energy structure found by GOAT"),
                    L("  orca.finalensemble.xyz   GOAT の配座集合。エネルギーと重みは output.log で確認してください", "  orca.finalensemble.xyz   GOAT conformer ensemble; check output.log for energies and weights")]
        if m.docker_guest_file and t == "single_point":
            files.append(L("  guest.xyz     利用者が指定した GUEST 構造 (コメント行の電荷・多重度も保持)", "  guest.xyz     user-supplied GUEST structure (including any charge and multiplicity in the comment line)"))
            out.append(L("  orca.docker.xyz   DOCKER の配置結果。相互作用エネルギーは output.log で確認してください", "  orca.docker.xyz   DOCKER poses; check output.log for interaction energies"))
            if m.docker_assume_neutral_singlet:
                prep.append(L("  GUEST の XYZ コメント行に電荷・多重度がないフレームは、ORCA の既定 (0, 1) として扱われます (明示指定で許可)。",
                              "  GUEST XYZ frames without charge and multiplicity in the comment use ORCA's default (0, 1), as explicitly allowed."))
        return ReadmeNotes(program="orca", files=files, prepare=prep, outputs=out)

    def orca_inp(self, spec: CalculationSpec) -> str:
        m, t, st, r = spec.method, spec.task, spec.structure, spec.runtime
        keys = [m.method.strip(), m.basis.strip()]
        if t.type == "geometry_optimization":
            keys.append("OptTS" if m.ts_search else "Opt")
            if m.ts_search and m.ts_freq:
                keys.append("Freq")
        elif t.type == "single_point" and m.irc:
            keys += ["Freq", "IRC"]
        elif t.type == "single_point" and m.goat:
            keys.append("GOAT")
        elif t.type == "vibrations":
            keys.append("Freq")
        elif t.type == "molecular_dynamics":
            keys.append("MD")
        if m.scf_convergence != "NormalSCF":
            keys.append(m.scf_convergence)
        if m.solvation != "none":
            keys.append(f"{m.solvation.upper()}({m.solvent.strip()})")
        if m.extra_keywords.strip():
            keys.append(m.extra_keywords.strip())
        lines = ["# ADIT が生成した ORCA の入力", "! " + " ".join(k for k in keys if k)]
        if r.mpiprocs > 1:
            lines += [f"%pal nprocs {r.mpiprocs} end"]
        if m.maxcore_mb > 0:
            lines += [f"%maxcore {m.maxcore_mb}"]
        if m.docker_guest_file:
            lines += ['%DOCKER GUEST "guest.xyz" END']
        lines += ["%scf", f"   MaxIter {m.scf_maxiter}", "end"]
        if t.type == "geometry_optimization":
            tol = t.force_tolerance_ev_per_ang / HARTREE_PER_BOHR_IN_EV_PER_ANG  # eV/Å → Eh/Bohr
            lines += ["%geom", f"   MaxIter {t.max_steps if t.max_steps > 0 else 200}", f"   TolMaxG {tol:.3e}"]
            if m.ts_search and m.ts_calc_hess:
                lines.append("   Calc_Hess true")
            if m.ts_search and m.ts_recalc_hess > 0:
                lines.append(f"   Recalc_Hess {m.ts_recalc_hess}")
            if st.fixed_atoms:
                lines += ["   Constraints"] + [f"      {{ C {i} C }}" for i in sorted(st.fixed_atoms)] + ["   end"]
            lines.append("end")
        if t.type == "single_point" and m.irc and (m.irc_max_iter > 0 or m.irc_direction != "both"):
            lines.append("%irc")
            if m.irc_max_iter > 0:
                lines.append(f"   MaxIter {m.irc_max_iter}")
            if m.irc_direction != "both":
                lines.append(f"   Direction {m.irc_direction}")
            lines.append("end")
        if t.type == "molecular_dynamics":
            md = t.md
            lines += ["%md", f"   Timestep {md.timestep_fs:g}_fs", f"   Initvel {md.temperature_k:g}_K"]
            if md.ensemble == "NVT":
                th = {"berendsen": "Berendsen", "csvr": "CSVR", "nose_hoover": "NHC"}.get(md.thermostat)
                if th is None:
                    raise GenerationError(L(f"ORCA にはない熱浴です: {md.thermostat} (berendsen / csvr / nose_hoover)", f"thermostat not available in ORCA: {md.thermostat} (berendsen / csvr / nose_hoover)"))
                lines.append(f"   Thermostat {th} {md.temperature_k:g}_K Timecon {md.coupling_time_fs:g}_fs")
            elif md.ensemble == "NPT":
                raise GenerationError(L("ORCA の MD に NPT はありません", "ORCA MD has no NPT"))
            # ORCA 6.1 manual, MD "Constraint" command: zero-based atom index, one atom per line.
            lines += [f"   Constraint Add Cartesian {i}" for i in sorted(st.fixed_atoms)]
            lines += [f'   Dump Position Stride {md.dump_interval} Filename "trajectory.xyz"', f"   Run {md.steps}", "end"]
        if m.extra_blocks.strip():
            lines += ["# 追加のブロック (extra_blocks。そのまま書く)", m.extra_blocks.strip()]
        lines.append(f"* xyz {st.charge} {st.multiplicity}")
        for sym, (x, y, z) in zip(st.atoms.symbols, st.atoms.positions):
            lines.append(f"  {sym:2s} {x:14.8f} {y:14.8f} {z:14.8f}")
        lines.append("*")
        return "\n".join(lines) + "\n"


register(OrcaGenerator())
