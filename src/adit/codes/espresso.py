
from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np

from ase.io.espresso import write_espresso_in

from adit.bandpath import KPATH_FILE, band_path, kpath_json, qe_crystal_b
from adit.citations import Citation
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.upf import UpfError, UpfLibrary
from adit.config import Config, Profile
from adit.spec import CalculationSpec, EspressoMethod
from adit.validate import electron_parity_error
from adit import lang
from adit.lang import L
from adit.validate_types import ValidationError

INPUT_FILE = "pw.in"
PH_FILE = "ph.in"
DYNMAT_FILE = "dynmat.in"
PSEUDO_SUBDIR = "pseudo"
DEFAULT_COMMAND = "mpirun -np {mpiprocs} pw.x"
RY_PER_BOHR_IN_EV_PER_ANG = 25.71104309541616  # 1 Ry/Bohr = 25.711 eV/Å
CELL_DOFREE = {"shape_and_volume": "all", "volume_only": "volume"}


def profile_cmd(cfg: Config, spec: CalculationSpec) -> str:
    p = cfg.profiles.get(spec.runtime.profile)
    return p.command_for("espresso", DEFAULT_COMMAND) if p else DEFAULT_COMMAND


def upf_name(m: EspressoMethod, lib: UpfLibrary | None, element: str) -> str | None:
    if element in m.pseudo:
        return m.pseudo[element]
    if lib is not None:
        cands = lib.files_for(element)
        if len(cands) == 1:
            return cands[0]
    return None


SSSP_URL = "https://www.materialscloud.org/discover/sssp"


def _pseudo_root_problem(pseudo_root: str, pseudo_set: str) -> str:
    from adit.config import config_path
    how = L(f"デスクトップ版では「擬ポテンシャルのセット」の欄の「フォルダを選ぶ…」で選び、ウェブ版とコマンド行では環境設定ファイル {config_path()} の "
            f"pseudo_root = \"...\" に書きます",
            f"In the desktop app, use \"Choose folder…\" next to Pseudopotential set; in the web version and on the command line, write it in "
            f"pseudo_root = \"...\" in the settings file {config_path()}")
    if not pseudo_root:
        return L(f"擬ポテンシャル (UPF) の置き場所 (環境設定の pseudo_root) がまだ決まっていません。{how} "
                 f"(例: pseudo_root = \"/home/<ユーザー名>/pseudo\"。その下にセットのフォルダ SSSP_efficiency/ などを置き、中に *.UPF。入手先の例: {SSSP_URL})",
                 f"The pseudopotential (UPF) folder (pseudo_root in the settings) is not set yet. {how} "
                 f"(e.g. pseudo_root = \"/home/<user>/pseudo\", containing set folders such as SSSP_efficiency/ with *.UPF inside; e.g. from {SSSP_URL})")
    root = Path(pseudo_root).expanduser()
    if not root.is_dir():
        return L(f"pseudo_root に書かれたフォルダがありません: {root}。{how}", f"the folder given as pseudo_root does not exist: {root}. {how}")
    found = ", ".join(sorted(d.name for d in root.iterdir() if d.is_dir() and any(f.suffix.lower() == ".upf" for f in d.iterdir()))) or L("(なし)", "(none)")
    return L(f"{root} にセット {pseudo_set!r} (UPF を入れたフォルダ) がありません (見つかったセット: {found})",
             f"{root} has no set {pseudo_set!r} (a folder of UPF files) (sets found: {found})")


def _cutoff_hint(spec: CalculationSpec, cfg: Config, ja: bool) -> str:
    lib = UpfLibrary.open_if_present(cfg.pseudo_root or None, spec.method.pseudo_set)
    if lib is None:
        return ""
    vals = []
    for e in spec.elements:
        name = upf_name(spec.method, lib, e)
        try:
            h = lib.header(name) if name and lib.has(name) else None
        except UpfError:
            h = None
        if h is not None and h.wfc_cutoff:
            vals.append(f"{e} {h.wfc_cutoff:g} Ry ({name})")
    if not vals:
        return ""
    return (" (この UPF の推奨値: " if ja else " (suggested in these UPFs: ") + ", ".join(vals) + ")"


class EspressoGenerator(InputGenerator):
    code = "espresso"

    def resolve(self, spec: CalculationSpec, cfg: Config) -> UpfLibrary:
        lib = UpfLibrary.open_if_present(cfg.pseudo_root or None, spec.method.pseudo_set)
        if lib is None:
            raise GenerationError(L(f"擬ポテンシャルのライブラリ (pseudo_root) に {spec.method.pseudo_set!r} がありません", f"pseudopotential library (pseudo_root) has no {spec.method.pseudo_set!r}"))
        return lib

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, EspressoMethod):
            return [ValidationError("method.code", L(f"pw.x の生成器に {m.code!r} の手法が渡されました", f"the pw.x generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        if not spec.structure.periodic:
            errs.append(ValidationError("structure.atoms", L("pw.x は周期セルが必要です (分子は構造の「周期セルに入れる」に印を付けてください)", "pw.x needs a periodic cell (for a molecule, tick \"Put in a periodic cell\" under Structure)")))
        if m.ecutwfc <= 0:
            errs.append(ValidationError("method.ecutwfc", L("値を入れてください (pw.x の必須項目。平面波の打ち切りエネルギー [Ry])。目安は UPF ファイルに書かれた推奨値 wfc_cutoff です"
                                                             + _cutoff_hint(spec, cfg, True),
                                                             "enter a value (required by pw.x: plane-wave cutoff in Ry). The value suggested in the UPF file (wfc_cutoff) is a guide"
                                                             + _cutoff_hint(spec, cfg, False))))
        if m.ecutrho < 0:
            errs.append(ValidationError("method.ecutrho", L("負の値は指定できません", "negative values are not allowed")))
        if m.conv_thr <= 0:
            errs.append(ValidationError("method.conv_thr", L("0 より大きい値が必要です", "must be greater than 0")))
        if m.electron_maxstep < 1:
            errs.append(ValidationError("method.electron_maxstep", L("1 以上が必要です", "must be at least 1")))
        if m.occupations == "smearing" and m.degauss <= 0:
            errs.append(ValidationError("method.degauss", L("smearing のときは 0 より大きい値が必要です", "must be greater than 0 when occupations = smearing")))
        bad_keys = [f"{ns}.{k}" for ns, kv in m.extra.items() for k in kv if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_()]*", k.strip())]
        if bad_keys:
            errs.append(ValidationError("method.extra", L(f"変数の名前として読めません: {bad_keys}", f"not valid variable names: {bad_keys}")))
        if m.dipole_correction:
            if m.dipole_direction not in (1, 2, 3):
                errs.append(ValidationError("method.dipole_direction", L(
                    "双極子補正を使うときは、真空をはさむ方向を 1 (a) / 2 (b) / 3 (c) から選んでください "
                    "(どの方向が真空かは構造しだいなので、ADIT は決めません)",
                    "choose the direction across the vacuum as 1 (a), 2 (b) or 3 (c) when using the dipole correction "
                    "(which direction holds the vacuum depends on the structure, so ADIT does not choose it)")))
            for field, value in (("dipole_maxpos", m.dipole_maxpos), ("dipole_decrease", m.dipole_decrease)):
                if not 0.0 <= value < 1.0:
                    errs.append(ValidationError(f"method.{field}", L(
                        "セルの分数座標なので 0 以上 1 未満です", "this is a fractional coordinate, so it must be at least 0 and less than 1")))
            if m.dipole_decrease <= 0:
                errs.append(ValidationError("method.dipole_decrease", L(
                    "ポテンシャルが下がる区間の幅 (eopreg) が 0 です。真空の中の、原子の無い区間を指定してください",
                    "the width over which the potential decreases (eopreg) is 0; give an interval inside the vacuum with no atoms")))
        if spec.structure.multiplicity != 1 and m.nspin != 2:
            errs.append(ValidationError("method.nspin", L("多重度が 1 でないなら nspin = 2 が必要です", "nspin = 2 is required when the multiplicity is not 1")))
        t = spec.task
        if t.type == "vibrations" and "pw.x" not in profile_cmd(cfg, spec):
            errs.append(ValidationError("runtime.profile", L("振動解析は ph.x を pw.x と同じ形で呼ぶので、commands.espresso に 'pw.x' の文字が必要です (既定のままで使えます)", "the vibrational analysis calls ph.x in the same way as pw.x, so commands.espresso must contain 'pw.x' (the default is fine)")))
        if t.type == "molecular_dynamics":
            if t.md.ensemble in ("NVT", "NPT") and t.md.thermostat not in ("berendsen", "andersen", "csvr"):
                errs.append(ValidationError("task.md.thermostat", L(f"pw.x に対応していない熱浴: {t.md.thermostat} (berendsen / andersen / csvr)", f"thermostat not supported for pw.x: {t.md.thermostat} (berendsen / andersen / csvr)")))
        lib = UpfLibrary.open_if_present(cfg.pseudo_root or None, m.pseudo_set)
        if lib is None:
            errs.append(ValidationError("method.pseudo_set", _pseudo_root_problem(cfg.pseudo_root, m.pseudo_set)))
            return errs
        zsum = 0.0
        for e in spec.elements:
            name = upf_name(m, lib, e)
            if name is None:
                cands = lib.files_for(e)
                msg = (
                    L(f"元素 {e} の UPF がライブラリにありません", f"no UPF for element {e} in the library")
                    if not cands
                    else L(f"元素 {e} の UPF が複数あります。どれを使うか指定してください: {cands}", f"several UPFs for element {e}; choose one: {cands}")
                )
                errs.append(ValidationError("method.pseudo", msg))
                continue
            if not lib.has(name):
                errs.append(ValidationError("method.pseudo", L(f"UPF がありません: {lib.set_dir / name}", f"UPF not found: {lib.set_dir / name}")))
                continue
            try:
                h = lib.header(name)
                if h.element.lower() != e.lower():
                    errs.append(ValidationError("method.pseudo", L(f"{name} は元素 {h.element} の UPF で、{e} に使えません", f"{name} is a UPF for {h.element} and cannot be used for {e}")))
                zsum += h.z_valence * spec.structure.atoms.symbols.count(e)
            except UpfError as ex:
                errs.append(ValidationError("method.pseudo", str(ex)))
        if not any(e.location == "method.pseudo" for e in errs):
            pe = electron_parity_error(int(round(zsum)) - spec.structure.charge, spec.structure.charge, spec.structure.multiplicity)
            if pe:
                errs.append(pe)
        errs += self._check_magnetism(spec, lib)
        return errs

    @staticmethod
    def _check_magnetism(spec: CalculationSpec, lib: UpfLibrary | None) -> list[ValidationError]:
        from adit.validate import check_hubbard

        m, errs = spec.method, []
        absent = [e for e in m.starting_magnetization if e not in spec.elements]
        if absent:
            errs.append(ValidationError("method.starting_magnetization", L(f"構造に無い元素です: {absent}", f"elements not in the structure: {absent}")))
        if m.nspin != 2 and any(v != 0 for v in m.starting_magnetization.values()):
            errs.append(ValidationError("method.nspin", L("starting_magnetization を入れたときは nspin = 2 が必要です (nspin = 1 では使われません)",
                                                          "nspin = 2 is required when starting_magnetization is given (it is not used with nspin = 1)")))
        errs += check_hubbard(m.hubbard, spec.elements, "method.hubbard")
        withj = [e for e, h in m.hubbard.items() if h.j_ev != 0]
        if withj:
            errs.append(ValidationError("method.hubbard", L(f"QE の生成器は U だけを書きます (J は書きません)。J を 0 にしてください: {withj}",
                                                            f"the QE generator writes U only (not J); set J to 0: {withj}")))
        if lib is not None and not errs:
            for e, h in m.hubbard.items():
                name = upf_name(m, lib, e)
                if not name or not lib.has(name):
                    continue
                labels = lib.chi_labels(name)
                want = h.orbital.strip().upper()
                if labels and want not in labels:
                    errs.append(ValidationError("method.hubbard", L(
                        f"{name} には {want} の原子軌道 (PP_CHI) がありません (ある軌道: {', '.join(labels)})。pw.x は HUBBARD の射影にこの軌道を使います",
                        f"{name} has no {want} atomic wavefunction (PP_CHI) (available: {', '.join(labels)}); pw.x uses it for the HUBBARD projectors")))
        return errs

    def version_probe(self, spec):
        return ("output.log", r"Program (PWSCF|NEB) v\.")

    def generate(self, spec: CalculationSpec, lib: UpfLibrary) -> dict[str, str]:
        out = {INPUT_FILE: self.pw_in(spec, lib)}
        if spec.task.type == "band_structure":
            kp = band_path(spec.atoms, spec.task.bands.path, spec.task.bands.npoints)
            out["bands/" + INPUT_FILE] = self.pw_in(spec, lib, bands=kp)
            out["bands/" + KPATH_FILE] = kpath_json(kp, spec.atoms)
        if spec.task.type == "vibrations":
            out[PH_FILE] = self.ph_in(spec)
            out[DYNMAT_FILE] = self.dynmat_in(spec)
        return out

    @staticmethod
    def ph_in(spec: CalculationSpec) -> str:
        return "phonons at Gamma (adit)\n&inputph\n   prefix = 'adit'\n   outdir = './tmp'\n   tr2_ph = 1.0d-14\n   fildyn = 'adit.dyn'\n/\n0.0 0.0 0.0\n"

    @staticmethod
    def dynmat_in(spec: CalculationSpec) -> str:
        atoms = spec.atoms
        pos = atoms.get_positions(); lengths = atoms.cell.lengths()
        spread = pos.max(axis=0) - pos.min(axis=0)
        isolated = all(L - s > 5.0 for L, s in zip(lengths, spread))
        asr = "zero-dim" if isolated else "crystal"
        return f"&input\n   fildyn = 'adit.dyn'\n   asr = '{asr}'\n   filout = 'dynmat.out'\n/\n"

    def files_to_copy(self, spec: CalculationSpec, lib: UpfLibrary) -> dict[str, Path]:
        out = {}
        for e in spec.elements:
            name = upf_name(spec.method, lib, e)
            if name is None or not lib.has(name):
                raise GenerationError(L(f"元素 {e} の UPF を決められません", f"cannot decide the UPF for element {e}"))
            out[f"{PSEUDO_SUBDIR}/{name}"] = lib.path(name)
        for n, p in lib.doc_files().items():
            out[f"{PSEUDO_SUBDIR}/{n}"] = p
        return out

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        cmd = profile.command_for(self.code, DEFAULT_COMMAND).format(mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary="")
        if spec.task.type == "vibrations":
            ph = cmd.replace("pw.x", "ph.x")
            return f"{cmd} -in {INPUT_FILE} > output.log 2>&1 && {ph} -in {PH_FILE} > ph.log 2>&1 && dynmat.x < {DYNMAT_FILE} > dynmat.log 2>&1"
        if spec.task.type == "band_structure":
            return f"{cmd} -in {INPUT_FILE} > output.log 2>&1 && cd bands && {cmd} -in {INPUT_FILE} > output.log 2>&1"
        return f"{cmd} -in {INPUT_FILE} > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, lib: UpfLibrary, copies) -> ReadmeNotes:
        t = spec.task.type
        names = ", ".join(sorted(n.split("/", 1)[1] for n in copies))
        docs = [n for n in ("LICENSE", "LICENSE.txt", "README", "README.md") if f"{PSEUDO_SUBDIR}/{n}" in copies]
        files = [
            L("  pw.in         pw.x の入力 (&CONTROL &SYSTEM などの設定と、元素・座標・k 点 = 周期系で電子の状態を計算する波数空間の点)",
              "  pw.in         pw.x input (&CONTROL, &SYSTEM etc., and the species, positions and k-points = points in reciprocal space where the electronic states are computed)"),
            L(f"  pseudo/       この計算に必要な擬ポテンシャル (内殻電子の効果をまとめた、元素ごとのファイル。UPF 形式): {names}",
              f"  pseudo/       pseudopotentials needed for this calculation (one file per element that stands in for the core electrons; UPF format): {names}"),
        ]
        if t == "vibrations":
            files.append(L("  ph.in / dynmat.in   振動解析の 2 段階目と 3 段階目。submit.sh が pw.x → ph.x → dynmat.x の順に実行します",
                           "  ph.in / dynmat.in   second and third stages of the vibrational analysis; submit.sh runs pw.x, ph.x and dynmat.x in that order"))
        if t == "band_structure":
            files.append(L("  bands/        バンド計算の 2 段階目。1 段階目の電荷密度 (tmp/) を読み、高対称点を結ぶ経路 (bands/kpath.json) に沿って計算します",
                           "  bands/        second stage of the band calculation: reads the charge density of the first stage (tmp/) and follows the high-symmetry path (bands/kpath.json)"))
        out = [
            L("  output.log    pw.x が画面に出す文字を保存したもの ('!    total energy' の行が全エネルギー。単位は Ry = リュードベリ)",
              "  output.log    what pw.x prints to the screen (the '!    total energy' lines give the total energy, in Ry = rydberg)"),
        ]
        out += {
            "geometry_optimization": [L("                最適化後の構造は output.log の末尾、Begin final coordinates から End final coordinates まで",
                                        "                the optimized structure is at the end of output.log, from Begin final coordinates to End final coordinates")],
            "molecular_dynamics": [L("                各ステップの座標 (ATOMIC_POSITIONS) と温度も output.log に出ます",
                                     "                the coordinates (ATOMIC_POSITIONS) and temperature of each step are also in output.log")],
            "vibrations": [L("  ph.log / dynmat.log   ph.x と dynmat.x のログ", "  ph.log / dynmat.log   logs of ph.x and dynmat.x"),
                           L("  dynmat.out    振動数 (cm⁻¹)", "  dynmat.out    frequencies (cm⁻¹)")],
            "band_structure": [L("  bands/output.log   経路上の各 k 点のエネルギー準位 (図は ADIT の解析タブか analyze.py で描けます)",
                                 "  bands/output.log   energy levels at each k-point on the path (plot them with the ADIT analysis tab or analyze.py)")],
        }.get(t, [])
        out.append(L("  tmp/          波動関数と電荷密度 (続きの計算に使う途中のファイル。大きくなることがあります)",
                     "  tmp/          wavefunctions and charge density (intermediate files for follow-up runs; can be large)"))
        if spec.method.dipole_correction:
            out.append(L(f"                双極子補正を入れています (edir = {spec.method.dipole_direction})。"
                         "output.log の「Computed dipole along edir」の行に、各ステップの双極子 [Debye] が出ます。"
                         "真空の区間 (emaxpos / eopreg) に原子が入っていないことを、構造で確かめてください。",
                         f"                the dipole correction is on (edir = {spec.method.dipole_direction}); the lines "
                         "'Computed dipole along edir' in output.log give the dipole per step in Debye. "
                         "Check on the structure that the vacuum region (emaxpos / eopreg) contains no atoms."))
        return ReadmeNotes(program="pw.x", files=files, outputs=out)

    # ---- pw.in ----
    def pw_in(self, spec: CalculationSpec, lib: UpfLibrary, bands=None) -> str:
        m, t, st = spec.method, spec.task, spec.structure
        atoms = spec.atoms
        pseudos = {e: upf_name(m, lib, e) for e in spec.elements}
        control = {"calculation": "scf", "prefix": "adit", "pseudo_dir": f"./{PSEUDO_SUBDIR}", "outdir": "./tmp", "tprnfor": True, "tstress": True}
        if bands is not None:
            control.update({"calculation": "bands", "pseudo_dir": f"../{PSEUDO_SUBDIR}", "outdir": "../tmp", "verbosity": "high"})
            control.pop("tprnfor"); control.pop("tstress")
        system = {"ecutwfc": m.ecutwfc, "occupations": m.occupations, "nspin": m.nspin}
        if m.assume_isolated:
            system["assume_isolated"] = m.assume_isolated
        if m.dipole_correction:
            control["tefield"] = True
            control["dipfield"] = True
            system["edir"] = int(m.dipole_direction)
            system["emaxpos"] = float(m.dipole_maxpos)
            system["eopreg"] = float(m.dipole_decrease)
            system["eamp"] = float(m.dipole_amplitude)
        if m.ecutrho > 0:
            system["ecutrho"] = m.ecutrho
        if m.occupations == "smearing":
            system["smearing"] = m.smearing; system["degauss"] = m.degauss
        if st.charge != 0:
            system["tot_charge"] = float(st.charge)
        if st.multiplicity != 1:
            system["tot_magnetization"] = float(st.multiplicity - 1)
        if m.input_dft:
            system["input_dft"] = m.input_dft
        if bands is not None:
            zsum = sum(lib.header(pseudos[e]).z_valence * st.atoms.symbols.count(e) for e in spec.elements)
            system["nbnd"] = int(np.ceil((zsum - st.charge) / 2.0)) + t.bands.empty_bands
        electrons = {"conv_thr": m.conv_thr, "electron_maxstep": m.electron_maxstep, "mixing_beta": m.mixing_beta}
        ions, cell = {}, {}
        if t.type == "molecular_dynamics":
            md = t.md  # INPUT_PW: calculation='md', dt [Ry a.u.] (1 a.u. = 4.8378e-17 s), nstep, ion_temperature, tempw, nraise
            control["calculation"] = "md"
            control["nstep"] = md.steps
            # INPUT_PW: &CONTROL iprint is the trajectory output interval in MD steps.
            control["iprint"] = md.dump_interval
            system["nosym"] = True
            control["dt"] = md.timestep_fs * 1e-15 / 4.8378e-17
            ions["ion_dynamics"] = "verlet"
            ions["tempw"] = md.temperature_k
            nraise = max(1, int(round(md.coupling_time_fs / md.timestep_fs)))
            if md.ensemble == "NVE":
                ions["ion_temperature"] = "initial"
            elif md.thermostat == "berendsen":
                ions["ion_temperature"] = "berendsen"; ions["nraise"] = nraise
            elif md.thermostat == "andersen":
                ions["ion_temperature"] = "andersen"; ions["nraise"] = nraise
            elif md.thermostat == "csvr":
                ions["ion_temperature"] = "svr"; ions["nraise"] = nraise
            else:
                raise GenerationError(L(f"pw.x の MD に対応していない熱浴: {md.thermostat} (berendsen / andersen / csvr)", f"thermostat not supported for pw.x MD: {md.thermostat} (berendsen / andersen / csvr)"))
            if md.ensemble == "NPT":
                control["calculation"] = "vc-md"
                ions["ion_dynamics"] = "beeman"
                cell["cell_dynamics"] = "pr"
                cell["press"] = md.pressure_bar / 1000.0
        if t.type == "geometry_optimization":
            control["calculation"] = "vc-relax" if t.relax_cell != "no" else "relax"
            control["nstep"] = t.max_steps if t.max_steps > 0 else 50
            control["forc_conv_thr"] = t.force_tolerance_ev_per_ang / RY_PER_BOHR_IN_EV_PER_ANG  # eV/Å → Ry/Bohr
            if t.relax_cell != "no":
                cell["cell_dofree"] = CELL_DOFREE[t.relax_cell]
        for ns, extra in m.extra.items():
            target = {"control": control, "system": system, "electrons": electrons, "ions": ions, "cell": cell}.get(ns.lower())
            if target is None:
                raise GenerationError(L(f"extra の名前空間 {ns!r} は control / system / electrons / ions / cell のどれかにしてください", f"extra namelist {ns!r} must be one of control / system / electrons / ions / cell"))
            target.update(extra)
        kp = spec.kpoints
        if kp is None or kp.mode == "gamma":
            kpts, koffset = (1, 1, 1), (0, 0, 0)
        else:
            kpts = kp.resolved_mesh(st.atoms.cell); koffset = tuple(1 if s == 0.5 else 0 for s in kp.shift)
        buf = io.StringIO()
        write_espresso_in(buf, atoms, input_data={"control": control, "system": system, "electrons": electrons, "ions": ions, "cell": cell},
                          pseudopotentials=pseudos, kpts=kpts, koffset=koffset, crystal_coordinates=True)
        text = self._with_if_pos(buf.getvalue(), spec)
        if m.assume_isolated and "assume_isolated" not in text:
            lines = text.splitlines(); start = next(i for i, line in enumerate(lines) if line.strip().upper() == "&SYSTEM")
            end = next(i for i in range(start + 1, len(lines)) if lines[i].strip() == "/")
            lines.insert(end, f"   assume_isolated = '{m.assume_isolated}'"); text = "\n".join(lines) + "\n"
        species = self._species_order(text)
        if species != spec.elements:
            raise GenerationError(L(f"ATOMIC_SPECIES の順 {species} が想定 {spec.elements} と違います", f"ATOMIC_SPECIES order {species} differs from the expected {spec.elements}"))
        mags = [(i, float(m.starting_magnetization[e])) for i, e in enumerate(species, 1) if m.starting_magnetization.get(e)]
        if mags:
            lines = text.splitlines()
            k = next(i for i, l in enumerate(lines) if l.strip().upper() == "&SYSTEM")
            end = next(i for i in range(k + 1, len(lines)) if lines[i].strip() == "/")
            own = {f"starting_magnetization({i})" for i, _ in mags}
            body = [l for l in lines[k + 1:end] if l.split("=")[0].strip().lower() not in own]
            lines = lines[:k + 1] + body + [f"   starting_magnetization({i}) = {v:g}" for i, v in mags] + lines[end:]
            text = "\n".join(lines) + "\n"
        if bands is not None:
            head, _, rest = text.partition("K_POINTS")
            after = rest.split("\n\n", 1)[1] if "\n\n" in rest else ""
            text = head + qe_crystal_b(bands) + "\n" + after
        if m.hubbard:
            text = text.rstrip("\n") + "\n\n" + f"HUBBARD {{{m.hubbard_projector}}}\n" + "".join(
                f"U {e}-{m.hubbard[e].orbital.strip().lower()} {float(m.hubbard[e].u_ev):g}\n" for e in spec.elements if e in m.hubbard)
        return text

    @staticmethod
    def _species_order(text: str) -> list[str]:
        lines = text.splitlines()
        k = next(i for i, l in enumerate(lines) if l.startswith("ATOMIC_SPECIES"))
        out = []
        for l in lines[k + 1:]:
            if not l.strip():
                break
            out.append(l.split()[0])
        return out

    @staticmethod
    def _with_if_pos(text: str, spec: CalculationSpec) -> str:
        st = spec.structure
        if not st.fixed_atoms and not st.fixed_axes:
            return text
        lines = text.splitlines()
        k = next(i for i, l in enumerate(lines) if l.startswith("ATOMIC_POSITIONS"))
        n = len(st.atoms.symbols)
        fixed = set(st.fixed_atoms)
        for j in range(n):
            move = (False, False, False) if j in fixed else tuple(st.fixed_axes.get(str(j), (True, True, True)))
            lines[k + 1 + j] = lines[k + 1 + j].rstrip() + " " + " ".join("1" if m else "0" for m in move)
        return "\n".join(lines) + "\n"


register(EspressoGenerator())


# References the Quantum ESPRESSO user guide asks for (its "Terms of use" section)
_QE_CITE_URL = "https://www.quantum-espresso.org/Doc/user_guide/node6.html"
CITATIONS = (
    Citation("qe_giannozzi2009", r"""@article{qe_giannozzi2009,
  author  = {Giannozzi, Paolo and Baroni, Stefano and Bonini, Nicola and Calandra, Matteo and Car, Roberto and Cavazzoni, Carlo and Ceresoli, Davide and Chiarotti, Guido L. and Cococcioni, Matteo and Dabo, Ismaila and Dal Corso, Andrea and de Gironcoli, Stefano and Fabris, Stefano and Fratesi, Guido and Gebauer, Ralph and Gerstmann, Uwe and Gougoussis, Christos and Kokalj, Anton and Lazzeri, Michele and Martin-Samos, Layla and Marzari, Nicola and Mauri, Francesco and Mazzarello, Riccardo and Paolini, Stefano and Pasquarello, Alfredo and Paulatto, Lorenzo and Sbraccia, Carlo and Scandolo, Sandro and Sclauzero, Gabriele and Seitsonen, Ari P. and Smogunov, Alexander and Umari, Paolo and Wentzcovitch, Renata M.},
  title   = {{QUANTUM ESPRESSO}: a modular and open-source software project for quantum simulations of materials},
  journal = {Journal of Physics: Condensed Matter},
  volume  = {21},
  number  = {39},
  pages   = {395502},
  year    = {2009},
  doi     = {10.1088/0953-8984/21/39/395502}
}""", doi="10.1088/0953-8984/21/39/395502", source=_QE_CITE_URL),
    Citation("qe_giannozzi2017", r"""@article{qe_giannozzi2017,
  author  = {Giannozzi, P. and Andreussi, O. and Brumme, T. and Bunau, O. and Buongiorno Nardelli, M. and Calandra, M. and Car, R. and Cavazzoni, C. and Ceresoli, D. and Cococcioni, M. and Colonna, N. and Carnimeo, I. and Dal Corso, A. and de Gironcoli, S. and Delugas, P. and DiStasio, R. A. and Ferretti, A. and Floris, A. and Fratesi, G. and Fugallo, G. and Gebauer, R. and Gerstmann, U. and Giustino, F. and Gorni, T. and Jia, J. and Kawamura, M. and Ko, H.-Y. and Kokalj, A. and K{\"u}{\c{c}}{\"u}kbenli, E. and Lazzeri, M. and Marsili, M. and Marzari, N. and Mauri, F. and Nguyen, N. L. and Nguyen, H.-V. and Otero-de-la-Roza, A. and Paulatto, L. and Ponc{\'e}, S. and Rocca, D. and Sabatini, R. and Santra, B. and Schlipf, M. and Seitsonen, A. P. and Smogunov, A. and Timrov, I. and Thonhauser, T. and Umari, P. and Vast, N. and Wu, X. and Baroni, S.},
  title   = {Advanced capabilities for materials modelling with {Quantum ESPRESSO}},
  journal = {Journal of Physics: Condensed Matter},
  volume  = {29},
  number  = {46},
  pages   = {465901},
  year    = {2017},
  doi     = {10.1088/1361-648X/aa8f79}
}""", doi="10.1088/1361-648X/aa8f79", source=_QE_CITE_URL),
)

# Pseudopotential families, looked up by a substring of method.pseudo_set (lower case)
PSEUDO_SET_CITATIONS: dict[str, tuple[Citation, ...]] = {
    "sssp": (Citation("sssp_prandini2018", r"""@article{sssp_prandini2018,
  author  = {Prandini, Gianluca and Marrazzo, Antimo and Castelli, Ivano E. and Mounet, Nicolas and Marzari, Nicola},
  title   = {Precision and efficiency in solid-state pseudopotential calculations},
  journal = {npj Computational Materials},
  volume  = {4},
  number  = {1},
  pages   = {72},
  year    = {2018},
  doi     = {10.1038/s41524-018-0127-2}
}""", doi="10.1038/s41524-018-0127-2", source="https://sssp.materialscloud.org/"),),
    "pslibrary": (Citation("pslibrary_dalcorso2014", r"""@article{pslibrary_dalcorso2014,
  author  = {Dal Corso, Andrea},
  title   = {Pseudopotentials periodic table: From {H} to {Pu}},
  journal = {Computational Materials Science},
  volume  = {95},
  pages   = {337--350},
  year    = {2014},
  doi     = {10.1016/j.commatsci.2014.07.043}
}""", doi="10.1016/j.commatsci.2014.07.043", source="https://dalcorso.github.io/pslibrary/"),),
}
