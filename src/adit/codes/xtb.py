
from __future__ import annotations

import io
from pathlib import Path

from ase.io import write

from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.config import Config, Profile
from adit.spec import CalculationSpec, XtbMethod
from adit.validate import electron_parity_error
from adit import lang
from adit.lang import L
from adit.validate_types import ValidationError

GEOMETRY_FILE = "struct.xyz"
CONTROL_FILE = "xtb.inp"
DEFAULT_COMMAND = "xtb"
OPT_LEVELS = ["crude", "sloppy", "loose", "normal", "tight", "verytight", "extreme"]

_ALPB = ("acetone acetonitrile aniline benzaldehyde benzene ch2cl2 chcl3 chloroform cs2 dioxane dmf dmso ether ethanol ethylacetate furane "
         "h2o hexadecane hexane methanol n-hexane nitromethane octanol phenol thf toluene water woctanol").split()
_GBSA_GFN1 = "acetone acetonitrile benzene ch2cl2 chcl3 chloroform cs2 dmso ether h2o methanol thf toluene water".split()
_GBSA_GFN2 = "acetone acetonitrile benzene ch2cl2 chcl3 chloroform cs2 dmf dmso ether h2o hexane methanol n-hexane thf toluene water".split()
SOLVENTS: dict[tuple[str, str], list[str]] = {
    ("alpb", "1"): _ALPB, ("alpb", "2"): _ALPB, ("alpb", "ff"): _ALPB,
    ("gbsa", "1"): _GBSA_GFN1, ("gbsa", "2"): _GBSA_GFN2, ("gbsa", "ff"): _GBSA_GFN2,
}


def solvent_names(model: str, gfn: str) -> list[str]:
    return list(SOLVENTS.get((model, gfn), []))


class XtbGenerator(InputGenerator):
    code = "xtb"

    def resolve(self, spec: CalculationSpec, cfg: Config):
        return None

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, XtbMethod):
            return [ValidationError("method.code", L(f"xtb の生成器に {m.code!r} の手法が渡されました", f"the xtb generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        if spec.structure.periodic:
            errs.append(ValidationError("structure.atoms", L("xtb の生成器 (v1) は分子系だけです。周期セルは扱いません", "the xtb generator (v1) is for molecules only; no periodic cells")))
        if m.accuracy <= 0:
            errs.append(ValidationError("method.accuracy", L("0 より大きい値が必要です", "must be greater than 0")))
        if m.etemp < 0:
            errs.append(ValidationError("method.etemp", L("負の温度は指定できません", "negative temperature is not allowed")))
        if m.max_iterations < 1:
            errs.append(ValidationError("method.max_iterations", L("1 以上が必要です", "must be at least 1")))
        if spec.structure.fixed_axes:
            errs.append(ValidationError("structure.fixed_axes", L("xtb の生成器は軸ごとの固定に対応していません (原子ごとの固定は $fix で行います)", "the xtb generator cannot fix individual axes (whole atoms are fixed via $fix)")))
        t = spec.task
        if t.type == "molecular_dynamics":
            if t.md.ensemble == "NPT":
                errs.append(ValidationError("task.md.ensemble", L("xtb の MD に NPT はありません (NVE か NVT にしてください)", "xtb MD has no NPT (NVE or NVT)")))
            if t.md.ensemble == "NVT" and t.md.thermostat != "berendsen":
                errs.append(ValidationError("task.md.thermostat", L("xtb の熱浴は Berendsen 固定です (公式文書の $md)。thermostat を berendsen にしてください", "the xtb thermostat is fixed to Berendsen ($md in the docs); set thermostat to berendsen")))
        if spec.task.type == "band_structure":
            errs.append(ValidationError("task.type", L("xtb にバンド計算はありません (周期系のコードで行ってください)", "xtb has no band structure (use a periodic code)")))
        errs += self._check_solvent(m)
        from ase.data import atomic_numbers
        n_el = sum(atomic_numbers[s] for s in spec.structure.atoms.symbols) - spec.structure.charge
        pe = electron_parity_error(n_el, spec.structure.charge, spec.structure.multiplicity)
        if pe:
            errs.append(pe)
        return errs

    @staticmethod
    def _check_solvent(m: XtbMethod) -> list[ValidationError]:
        name = m.solvent.strip()
        if m.solvation == "none":
            if name:
                return [ValidationError("method.solvent", L(f"溶媒モデルが「なし」なので、溶媒 {name!r} は使われません (モデルを alpb か gbsa にするか、溶媒を空にしてください)",
                                                            f"the solvation model is 'none', so the solvent {name!r} would not be used (choose alpb or gbsa, or clear the solvent)"))]
            return []
        method = {"0": "GFN0-xTB", "1": "GFN1-xTB", "2": "GFN2-xTB", "ff": "GFN-FF"}[m.gfn]
        names = SOLVENTS.get((m.solvation, m.gfn), [])
        if not names:
            return [ValidationError("method.solvation", L(f"{method} では --{m.solvation} を使えません (xtb が止まります。ALPB と GBSA は GFN1・GFN2・GFN-FF 用)",
                                                          f"--{m.solvation} cannot be used with {method} (xtb stops; ALPB and GBSA are for GFN1, GFN2 and GFN-FF)"))]
        if not name:
            return [ValidationError("method.solvent", L(f"溶媒の名前を入れてください (xtb が {method} の --{m.solvation} で受け付ける名前: {', '.join(names)})",
                                                        f"enter the solvent name (names xtb accepts for --{m.solvation} with {method}: {', '.join(names)})"))]
        if name.lower() not in names:
            return [ValidationError("method.solvent", L(f"xtb は {method} の --{m.solvation} で溶媒 {name!r} を受け付けません (止まります)。使える名前: {', '.join(names)}",
                                                        f"xtb does not accept the solvent {name!r} for --{m.solvation} with {method} (it stops). Available: {', '.join(names)}"))]
        return []

    def version_probe(self, spec):
        return ("output.log", "xtb version")

    def generate(self, spec: CalculationSpec, res) -> dict[str, str]:
        return {GEOMETRY_FILE: self.xyz(spec), CONTROL_FILE: self.control(spec)}

    def files_to_copy(self, spec, res) -> dict[str, Path]:
        return {}

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        m, st, t = spec.method, spec.structure, spec.task
        exe = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        args = [exe, GEOMETRY_FILE, "--gfnff" if m.gfn == "ff" else f"--gfn {m.gfn}", f"--chrg {st.charge}", f"--uhf {st.multiplicity - 1}",
                f"--acc {m.accuracy:g}", f"--etemp {m.etemp:g}", f"--input {CONTROL_FILE}", f"--parallel {spec.runtime.omp_threads}"]
        if m.solvation != "none":
            args.append(f"--{m.solvation} {m.solvent.strip().lower()}")
        if t.type == "geometry_optimization":
            args.append(f"--opt {m.opt_level}")
        elif t.type == "molecular_dynamics":
            args.append("--md")
        elif t.type == "vibrations":
            args.append("--hess")
        args.append("--json")
        return " ".join(args) + " > output.log 2>&1"

    def readme_notes(self, spec, res, copies) -> ReadmeNotes:
        t = spec.task.type
        files = [
            L("  struct.xyz    構造 (原子の種類と座標。xyz 形式、長さの単位は Å)", "  struct.xyz    structure (elements and coordinates; xyz format, lengths in Å)"),
            L("  xtb.inp       xtb への追加の指示 (SCC の最大反復回数、最適化の最大サイクル数、固定原子)",
              "  xtb.inp       extra instructions for xtb (max SCC iterations, max optimization cycles, fixed atoms)"),
            L("                計算手法・電荷・多重度・精度・温度は submit.sh の xtb の行に引数として書いてあります",
              "                method, charge, multiplicity, accuracy and temperature are arguments on the xtb line of submit.sh"),
        ]
        out = [
            L("  output.log    xtb が画面に出す文字を保存したもの (TOTAL ENERGY の行が全エネルギー、単位は Eh = ハートリー)",
              "  output.log    what xtb prints to the screen (the TOTAL ENERGY line is the total energy, in Eh = hartree)"),
            L("  xtbout.json   結果をプログラムで読みやすい形で書いたもの", "  xtbout.json   results in a form that programs can read"),
        ]
        out += {
            "geometry_optimization": [L("  xtbopt.xyz    最適化後の構造 (分子ビューアで開けます)", "  xtbopt.xyz    optimized structure (opens in molecular viewers)")],
            "molecular_dynamics": [L("  xtb.trj       MD の軌跡 (xyz 形式、dump の間隔ごと)。mdrestart は続きから実行するためのもの",
                                     "  xtb.trj       MD trajectory (xyz format, every dump interval); mdrestart is for continuing the run")],
            "vibrations": [L("  vibspectrum / g98.out   振動数 (cm⁻¹) と IR 強度", "  vibspectrum / g98.out   frequencies (cm⁻¹) and IR intensities")],
        }.get(t, [])
        prepare = []
        if t == "molecular_dynamics":
            m = spec.method
            prepare = [L(f"  xtb.inp には水素の質量 hmass = {m.md_hmass:g} u、結合拘束 shake = {m.md_shake}、MD 中の SCC 精度 sccacc = {m.md_sccacc:g} を明記しています。",
                         f"  xtb.inp explicitly records the hydrogen mass hmass = {m.md_hmass:g} u, bond constraint setting shake = {m.md_shake}, and MD SCC accuracy sccacc = {m.md_sccacc:g}.")]
        return ReadmeNotes(program="xtb", files=files, prepare=prepare, outputs=out)

    def xyz(self, spec: CalculationSpec) -> str:
        buf = io.StringIO()
        write(buf, spec.atoms, format="xyz")
        return buf.getvalue()

    def control(self, spec: CalculationSpec) -> str:
        m, t, st = spec.method, spec.task, spec.structure
        lines = ["# ADIT が生成した xcontrol の指示", "$scc", f"   maxiterations={m.max_iterations}"]
        if t.type == "geometry_optimization":
            lines += ["$opt", f"   maxcycle={t.max_steps if t.max_steps > 0 else 1000}"]
        if t.type == "molecular_dynamics":
            md = t.md  # https://xtb-docs.readthedocs.io/en/latest/md.html: temp [K], time [ps], dump [fs], step [fs], nvt
            lines += ["$md", f"   temp={md.temperature_k:g}", f"   time={md.steps * md.timestep_fs / 1000.0:g}",
                      f"   dump={md.dump_interval * md.timestep_fs:g}", f"   step={md.timestep_fs:g}",
                      f"   nvt={'true' if md.ensemble == 'NVT' else 'false'}", f"   hmass={m.md_hmass:g}",
                      f"   shake={m.md_shake}", f"   sccacc={m.md_sccacc:g}"]
            if spec.handoff is not None and spec.handoff.velocities:
                lines.append("   restart=true")
        if st.fixed_atoms:
            lines += ["$fix", "   atoms: " + ",".join(str(i + 1) for i in sorted(st.fixed_atoms))]
        lines.append("$end")
        return "\n".join(lines) + "\n"


register(XtbGenerator())


# References from the "Citations" section of the xtb README; GFN_CITATIONS is keyed by method.gfn
_XTB_CITE_URL = "https://github.com/grimme-lab/xtb#citations"
CITATIONS = (
    Citation("xtb_bannwarth2020", r"""@article{xtb_bannwarth2020,
  author  = {Bannwarth, Christoph and Caldeweyher, Eike and Ehlert, Sebastian and Hansen, Andreas and Pracht, Philipp and Seibert, Jakob and Spicher, Sebastian and Grimme, Stefan},
  title   = {Extended tight-binding quantum chemistry methods},
  journal = {WIREs Computational Molecular Science},
  volume  = {11},
  number  = {2},
  pages   = {e1493},
  year    = {2020},
  doi     = {10.1002/wcms.1493}
}""", doi="10.1002/wcms.1493", source=_XTB_CITE_URL),
)
GFN_CITATIONS: dict[str, tuple[Citation, ...]] = {
    "1": (Citation("xtb_grimme2017", r"""@article{xtb_grimme2017,
  author  = {Grimme, Stefan and Bannwarth, Christoph and Shushkov, Philip},
  title   = {A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements ({Z} = 1--86)},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {13},
  number  = {5},
  pages   = {1989--2009},
  year    = {2017},
  doi     = {10.1021/acs.jctc.7b00118}
}""", doi="10.1021/acs.jctc.7b00118", source=_XTB_CITE_URL),),
    "2": (Citation("xtb_bannwarth2019", r"""@article{xtb_bannwarth2019,
  author  = {Bannwarth, Christoph and Ehlert, Sebastian and Grimme, Stefan},
  title   = {{GFN2-xTB}---An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {15},
  number  = {3},
  pages   = {1652--1671},
  year    = {2019},
  doi     = {10.1021/acs.jctc.8b01176}
}""", doi="10.1021/acs.jctc.8b01176", source=_XTB_CITE_URL),),
    "ff": (Citation("xtb_spicher2020", r"""@article{xtb_spicher2020,
  author  = {Spicher, Sebastian and Grimme, Stefan},
  title   = {Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems},
  journal = {Angewandte Chemie International Edition},
  volume  = {59},
  number  = {36},
  pages   = {15665--15673},
  year    = {2020},
  doi     = {10.1002/anie.202004239}
}""", doi="10.1002/anie.202004239", source=_XTB_CITE_URL),),
}
