
from __future__ import annotations

import io
from pathlib import Path

from ase.io.vasp import write_vasp

from adit.bandpath import KPATH_FILE, band_path, kpath_json, vasp_line_mode
from adit.codes.base import GenerationError, InputGenerator, ReadmeNotes, register
from adit.codes.potcar import PotcarError, PotcarLibrary
from adit.codes.potcar_names import MP_POTCAR_NAMES
from adit.config import Config, Profile
from adit.spec import CalculationSpec, VaspMethod
from adit.validate import electron_parity_error
from adit import lang
from adit.lang import L
from adit.validate_types import ValidationError
from adit.vasp_constraints import VaspConstraintError, cartesian_to_direct_mask

PP_ENV = "VASP_PP_PATH"
POTCAR_SPEC = "potcar.spec"
MAKE_POTCAR = "make_potcar.sh"
ISIF_OF = {"no": 2, "shape_and_volume": 3, "volume_only": 7}
DEFAULT_COMMAND = "mpirun -np {mpiprocs} vasp_{binary}"


def potcar_name(m: VaspMethod, element: str) -> str:
    return m.potcar.get(element) or MP_POTCAR_NAMES.get(element) or element


def _inc(v) -> str:
    if isinstance(v, bool):
        return ".TRUE." if v else ".FALSE."
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def md_unapplied_settings(spec: CalculationSpec) -> dict:
    if spec.task.type != "molecular_dynamics":
        return {}
    md = spec.task.md
    extra = {key.strip().upper(): value for key, value in spec.method.extra_incar.items()}
    unused = {}
    if md.ensemble != "NVE" and md.thermostat == "nose_hoover":
        unused["task.md.coupling_time_fs"] = {
            "value": md.coupling_time_fs,
            "incar_parameters": {"SMASS": extra.get("SMASS")},
            "reason": L(
                "Nosé–Hoover の時定数は SMASS に換算していません。SMASS は method.extra_incar の値を使います。",
                "The Nosé–Hoover coupling time is not converted to SMASS; SMASS comes from method.extra_incar."),
        }
    if md.ensemble == "NPT":
        unused["task.md.barostat_time_fs"] = {
            "value": md.barostat_time_fs,
            "incar_parameters": {key: extra.get(key) for key in ("PMASS", "LANGEVIN_GAMMA_L")},
            "reason": L(
                "圧力浴の時定数は PMASS や LANGEVIN_GAMMA_L に換算していません。これらは method.extra_incar の値を使います。",
                "The barostat coupling time is not converted to PMASS or LANGEVIN_GAMMA_L; both come from method.extra_incar."),
        }
    return unused


class VaspGenerator(InputGenerator):
    code = "vasp"

    def resolve(self, spec: CalculationSpec, cfg: Config) -> PotcarLibrary | None:
        profile = cfg.profiles.get(spec.runtime.profile)
        root = (profile.env.get(PP_ENV) if profile else None) or None
        return PotcarLibrary.open_if_present(root, spec.method.potcar_set)

    def validate(self, spec: CalculationSpec, cfg: Config) -> list[ValidationError]:
        m = spec.method
        if not isinstance(m, VaspMethod):
            return [ValidationError("method.code", L(f"VASP の生成器に {m.code!r} の手法が渡されました", f"the VASP generator received a {m.code!r} method"))]
        errs: list[ValidationError] = []
        if not spec.structure.periodic:
            errs.append(ValidationError("structure.atoms", L("VASP は周期セルが必要です (分子は構造の「周期セルに入れる」に印を付けてください)", "VASP needs a periodic cell (for a molecule, tick \"Put in a periodic cell\" under Structure)")))
        try:
            self._move_flags(spec, list(range(len(spec.structure.atoms.symbols))))
        except VaspConstraintError as ex:
            errs.append(ValidationError("structure.fixed_axes", str(ex)))
        if m.encut < 0:
            errs.append(ValidationError("method.encut", L("負の値は指定できません", "negative values are not allowed")))
        if m.ediff <= 0:
            errs.append(ValidationError("method.ediff", L("0 より大きい値が必要です", "must be greater than 0")))
        if m.nelm < 1:
            errs.append(ValidationError("method.nelm", L("1 以上が必要です", "must be at least 1")))
        if m.sigma < 0:
            errs.append(ValidationError("method.sigma", L("負の値は指定できません", "negative values are not allowed")))
        if m.magmom is not None and len(m.magmom) != len(spec.structure.atoms.symbols):
            errs.append(ValidationError("method.magmom", L(f"MAGMOM の数 {len(m.magmom)} が原子数 {len(spec.structure.atoms.symbols)} と違います", f"MAGMOM has {len(m.magmom)} values but there are {len(spec.structure.atoms.symbols)} atoms")))
        if spec.structure.multiplicity != 1 and m.ispin != 2:
            errs.append(ValidationError("method.ispin", L("多重度が 1 でないなら ISPIN = 2 が必要です", "ISPIN = 2 is required when the multiplicity is not 1")))
        errs += self._check_magnetism(spec)
        t = spec.task
        if t.type == "molecular_dynamics" and t.md.ensemble != "NVE":
            if t.md.thermostat not in ("andersen", "langevin", "nose_hoover", "csvr", "berendsen"):
                errs.append(ValidationError("task.md.thermostat", L(f"VASP に対応していない熱浴: {t.md.thermostat} (andersen / langevin / nose_hoover / csvr)", f"thermostat not supported for VASP: {t.md.thermostat} (andersen / langevin / nose_hoover / csvr)")))
            if t.md.thermostat == "berendsen":
                errs.append(ValidationError("task.md.thermostat", L("VASP に Berendsen 熱浴はありません (andersen / langevin / nose_hoover / csvr から選んでください)", "VASP has no Berendsen thermostat (andersen / langevin / nose_hoover / csvr)")))
            if t.md.thermostat == "nose_hoover" and "SMASS" not in {k.upper() for k in m.extra_incar}:
                errs.append(ValidationError("method.extra_incar", L("nose_hoover (MDALGO = 2) には SMASS が必要です。追加の INCAR に SMASS = <値> を書いてください (時定数からは機械的に決められません)", "nose_hoover (MDALGO = 2) needs SMASS; add SMASS = <value> to the extra INCAR (it cannot be derived from the time constant)")))
            if t.md.ensemble == "NPT":
                if t.md.thermostat != "langevin":
                    errs.append(ValidationError("task.md.thermostat", L("VASP の NPT は langevin (MDALGO = 3) と ISIF = 3 の組み合わせです", "VASP NPT requires langevin (MDALGO = 3) with ISIF = 3")))
                missing = [k for k in ("LANGEVIN_GAMMA_L", "PMASS") if k not in {kk.upper() for kk in m.extra_incar}]
                if missing:
                    errs.append(ValidationError("method.extra_incar", L(f"NPT (MDALGO = 3, ISIF = 3) には {missing} が必要です。追加の INCAR に書いてください", f"NPT (MDALGO = 3, ISIF = 3) needs {missing}; add them to the extra INCAR")))
        bad_keys = [k for k in m.extra_incar if not k.strip() or " " in k.strip()]
        if bad_keys:
            errs.append(ValidationError("method.extra_incar", L(f"キーとして解釈できません: {bad_keys}", f"not valid keys: {bad_keys}")))
        lib = self.resolve(spec, cfg)
        if lib is None:
            if spec.structure.charge != 0:
                errs.append(ValidationError("structure.charge", L(f"電荷 {spec.structure.charge} を NELECT に直すには POTCAR ライブラリ (プロファイルの env.{PP_ENV}) が必要ですが、この PC にはありません", f"converting charge {spec.structure.charge} to NELECT needs the POTCAR library (env.{PP_ENV} of the profile), which is not on this machine")))
            return errs
        zsum = 0.0
        for e in spec.elements:
            name = potcar_name(m, e)
            if not lib.has(name):
                errs.append(ValidationError("method.potcar", L(f"POTCAR ライブラリ {lib.set_dir} に {name!r} がありません (元素 {e})", f"POTCAR library {lib.set_dir} has no {name!r} (element {e})")))
                continue
            try:
                zsum += lib.header(name).zval * spec.structure.atoms.symbols.count(e)
            except PotcarError as ex:
                errs.append(ValidationError("method.potcar", str(ex)))
        if not any(e.location == "method.potcar" for e in errs):
            n_el = int(round(zsum)) - spec.structure.charge
            pe = electron_parity_error(n_el, spec.structure.charge, spec.structure.multiplicity)
            if pe:
                errs.append(pe)
        return errs

    @staticmethod
    def _check_magnetism(spec: CalculationSpec) -> list[ValidationError]:
        from adit.validate import check_hubbard

        m, errs = spec.method, []
        if m.magmom is not None and m.magmom_by_element:
            errs.append(ValidationError("method.magmom_by_element", L("原子ごとの MAGMOM と元素ごとの磁気モーメントは、どちらか一方にしてください",
                                                                      "give either per-atom MAGMOM or per-element magnetic moments, not both")))
        absent = [e for e in m.magmom_by_element if e not in spec.elements]
        if absent:
            errs.append(ValidationError("method.magmom_by_element", L(f"構造に無い元素です: {absent}", f"elements not in the structure: {absent}")))
        if m.ispin != 2 and (m.magmom is not None or m.magmom_by_element):
            errs.append(ValidationError("method.ispin", L("磁気モーメントを入れたときは ISPIN = 2 が必要です (ISPIN = 1 では MAGMOM が使われません)",
                                                          "ISPIN = 2 is required when magnetic moments are given (with ISPIN = 1, MAGMOM is not used)")))
        errs += check_hubbard(m.hubbard, spec.elements, "method.hubbard")
        return errs

    def generate(self, spec: CalculationSpec, lib: PotcarLibrary | None) -> dict[str, str]:
        out = {
            "POSCAR": self.poscar(spec),
            "INCAR": self.incar(spec, lib),
            "KPOINTS": self.kpoints(spec),
            POTCAR_SPEC: self.potcar_spec(spec, lib),
            MAKE_POTCAR: self.make_potcar_sh(spec),
        }
        if spec.task.type == "band_structure":
            kp = band_path(spec.atoms, spec.task.bands.path, spec.task.bands.npoints)
            nseg = max(1, len(kp.segments()))
            out["bands/POSCAR"] = out["POSCAR"]
            out["bands/INCAR"] = self.incar(spec, lib) + "\n# バンド計算 (親の CHGCAR を読む。wiki ICHARG = 11)\nICHARG = 11\n"
            out["bands/KPOINTS"] = vasp_line_mode(kp, max(2, spec.task.bands.npoints // nseg))
            out["bands/" + KPATH_FILE] = kpath_json(kp, spec.atoms)
        return out

    def files_to_copy(self, spec: CalculationSpec, lib) -> dict[str, Path]:
        return {}

    writes_velocities = True

    def version_probe(self, spec):
        return ("vasprun.xml", "name=.version.")

    def run_command(self, spec: CalculationSpec, profile: Profile) -> str:
        cmd = profile.command_for(self.code, DEFAULT_COMMAND).format(
            mpiprocs=spec.runtime.mpiprocs, omp_threads=spec.runtime.omp_threads, binary=spec.method.binary)
        if spec.task.type == "band_structure":
            return f"bash {MAKE_POTCAR} && {cmd} > output.log 2>&1 && cp CHGCAR POTCAR bands/ && cd bands && {cmd} > output.log 2>&1"
        return f"bash {MAKE_POTCAR} && {cmd} > output.log 2>&1"

    def readme_notes(self, spec: CalculationSpec, lib, copies) -> ReadmeNotes:
        m, t = spec.method, spec.task.type
        names = ", ".join(f"{e} → {potcar_name(m, e)}" for e in spec.elements)
        files = [
            L("  INCAR         VASP の入力 (計算の設定)", "  INCAR         VASP input (the settings)"),
            L("  KPOINTS       k 点 (周期系で電子の状態を計算する、波数空間の点の取り方)",
              "  KPOINTS       k-points (the points in reciprocal space where the electronic states of a periodic system are computed)"),
            L(f"  potcar.spec   使う POTCAR (元素ごとの擬ポテンシャル) の一覧 ({names})。POTCAR 自体は入っていません",
              f"  potcar.spec   list of POTCARs (pseudopotentials per element) to use ({names}); the POTCAR files are not included"),
            L("  make_potcar.sh  上の一覧の順に POTCAR をつなげて 1 つのファイルにします (submit.sh が最初に呼びます)",
              "  make_potcar.sh  joins the POTCARs in the above order into one file (called first by submit.sh)"),
        ]
        if "IMAGES" not in {k.strip().upper() for k in m.extra_incar}:
            files.insert(1, L("  POSCAR        構造 (元素ごとにまとめ直してあり、順番は potcar.spec と同じです)",
                              "  POSCAR        structure (regrouped by element, in the same order as potcar.spec)"))
        if t == "band_structure":
            files.append(L("  bands/        バンド計算の 2 段階目。1 段階目の電荷密度 (CHGCAR) を読み、高対称点を結ぶ経路に沿って計算します",
                           "  bands/        second stage of the band calculation: reads the charge density (CHGCAR) of the first stage and follows the high-symmetry path"))
        prep = [
            L("  POTCAR は VASP のライセンス保持者にしか配布できないため、このディレクトリには入っていません。",
              "  POTCAR may only be distributed to VASP license holders, so it is not in this directory."),
            L(f"  計算を実行する計算機で、環境変数 {PP_ENV} (POTCAR ライブラリの親ディレクトリ。その下に {m.potcar_set}/<名前>/POTCAR がある場所) を設定しておいてください。",
              f"  On the machine that runs the job, set the environment variable {PP_ENV} (the parent of the POTCAR library, containing {m.potcar_set}/<name>/POTCAR)."),
            L("  ADIT の環境設定 (cluster.toml のプロファイルの env) に書いておくと、submit.sh が設定します。",
              "  If it is written in the ADIT settings (env of the profile in cluster.toml), submit.sh sets it."),
            (L("  この PC では POTCAR ライブラリを確認できなかったため、potcar.spec の TITEL / ZVAL は未確認です (実行時に make_potcar.sh が potcar.used に記録します)。",
               "  The library was not available on this machine, so TITEL / ZVAL in potcar.spec are unverified (make_potcar.sh records them in potcar.used at run time).")
             if lib is None else L(f"  この PC のライブラリ {lib.set_dir} でヘッダを確認済みです (TITEL / ZVAL を potcar.spec に記録)。",
                                   f"  Headers were checked against the library {lib.set_dir} on this machine (TITEL / ZVAL recorded in potcar.spec).")),
        ]
        for field, detail in md_unapplied_settings(spec).items():
            prep.append(L(
                f"  {field} = {detail['value']:g} fs は spec.json に残りますが、VASP 入力には適用されません。{detail['reason']}",
                f"  {field} = {detail['value']:g} fs remains in spec.json but is not applied to VASP input. {detail['reason']}"))
        out = [
            L("  output.log    VASP が画面に出す文字を保存したもの (電子の反復 1 回ごとの行)", "  output.log    what VASP prints to the screen (one line per electronic iteration)"),
            L("  OSZICAR       各ステップのエネルギーの要約 (単位は eV)", "  OSZICAR       summary of the energy of each step, in eV"),
            L("  OUTCAR        詳しい出力 (エネルギー、力、各ステップ)", "  OUTCAR        detailed output (energy, forces, each step)"),
            L("  vasprun.xml   結果をプログラムで読みやすい形で書いたもの (ASE や pymatgen で読めます)",
              "  vasprun.xml   results in a form that programs can read (ASE, pymatgen)"),
        ]
        out += {
            "geometry_optimization": [L("  CONTCAR       最適化後の構造 (POSCAR と同じ形式)", "  CONTCAR       optimized structure (same format as POSCAR)")],
            "molecular_dynamics": [L("  XDATCAR       MD の軌跡 (NBLOCK ステップごとの座標)。NBLOCK は task.md.dump_interval から設定。CONTCAR は最後の構造",
                                     "  XDATCAR       MD trajectory (coordinates every NBLOCK steps). NBLOCK is set from task.md.dump_interval; CONTCAR is the last structure")],
            "vibrations": [L("  OUTCAR の後ろのほうの「Eigenvectors and eigenvalues of the dynamical matrix」の節に振動数があります",
                             "  the frequencies are in the 'Eigenvectors and eigenvalues of the dynamical matrix' section near the end of OUTCAR")],
            "band_structure": [L("  bands/EIGENVAL, bands/vasprun.xml   経路上の各 k 点のエネルギー準位 (図は ADIT の解析タブか analyze.py で描けます)",
                                 "  bands/EIGENVAL, bands/vasprun.xml   energy levels at each k-point on the path (plot them with the ADIT analysis tab or analyze.py)")],
        }.get(t, [])
        return ReadmeNotes(program=f"vasp_{m.binary}", files=files, prepare=prep, outputs=out)

    @staticmethod
    def _ordered_atoms(spec: CalculationSpec):
        atoms = spec.atoms
        symbols = list(atoms.get_chemical_symbols())
        order = [i for e in spec.elements for i, s in enumerate(symbols) if s == e]
        return atoms[order], order

    @staticmethod
    def _move_flags(spec: CalculationSpec, order: list[int]) -> list[tuple[bool, bool, bool]] | None:
        st = spec.structure
        if not st.fixed_atoms and not st.fixed_axes:
            return None
        fixed = set(st.fixed_atoms)
        out = []
        for i in order:
            if i in fixed:
                out.append((False, False, False))
            else:
                movable = st.fixed_axes.get(str(i), (True, True, True))
                mask = cartesian_to_direct_mask(st.atoms.cell, [not value for value in movable])
                out.append(tuple(bool(value) for value in ~mask))
        if any(any(move) and not all(move) for move in out):
            t = spec.task
            isif = (ISIF_OF[t.relax_cell] if t.type == "geometry_optimization" else
                    3 if t.type == "molecular_dynamics" and t.md.ensemble == "NPT" else 2)
            extra = {key.strip().upper(): value for key, value in spec.method.extra_incar.items()}
            isif = _inc(extra.get("ISIF", isif)).split("#", 1)[0].strip()
            if isif in {"3", "4", "5", "6"}:
                raise VaspConstraintError(L(
                    "セル形状を動かす ISIF では Cartesian 軸の固定を保持できません。"
                    "POSCAR の固定方向は変形する格子ベクトルに従います。",
                    "Cartesian axis constraints cannot be preserved when ISIF allows the cell shape to change. "
                    "POSCAR fixed directions follow the deforming lattice vectors."))
        return out

    def poscar(self, spec: CalculationSpec) -> str:
        atoms, order = self._ordered_atoms(spec)
        buf = io.StringIO()
        write_vasp(buf, atoms, direct=True, sort=False, vasp5=True, ignore_constraints=True,
                   symbol_count=[(e, spec.structure.atoms.symbols.count(e)) for e in spec.elements])
        text = buf.getvalue().splitlines()
        text[0] = f"adit: {' '.join(spec.elements)} ({len(atoms)} atoms)"
        flags = self._move_flags(spec, order)
        if flags is not None:
            k = next(i for i, l in enumerate(text) if l.strip().lower().startswith(("direct", "cartesian")))
            text.insert(k, "Selective dynamics")
            coords = text[k + 2: k + 2 + len(atoms)]
            text[k + 2: k + 2 + len(atoms)] = [f"{c.rstrip()}   {'T' if fx else 'F'}   {'T' if fy else 'F'}   {'T' if fz else 'F'}"
                                               for c, (fx, fy, fz) in zip(coords, flags)]
        text = [l for l in text if l.strip()]
        vel = spec.structure.velocities
        if vel is not None:
            text.append("")
            text += [f"  {vel[i][0]:.10e}  {vel[i][1]:.10e}  {vel[i][2]:.10e}" for i in order]
        return "\n".join(text) + "\n"

    def incar(self, spec: CalculationSpec, lib: PotcarLibrary | None) -> str:
        m, t, st = spec.method, spec.task, spec.structure
        lines = [f"SYSTEM = adit {' '.join(spec.elements)}", "", "# 電子状態"]
        if m.encut > 0:
            lines.append(f"ENCUT = {_inc(m.encut)}")
        lines += [f"PREC = {m.prec}", f"ALGO = {m.algo}", f"EDIFF = {_inc(m.ediff)}", f"NELM = {m.nelm}"]
        if m.nelmin > 0:
            lines.append(f"NELMIN = {m.nelmin}")
        lines += [f"ISMEAR = {m.ismear}", f"SIGMA = {_inc(m.sigma)}", f"LREAL = {m.lreal}", f"ISPIN = {m.ispin}"]
        if m.lasph:
            lines.append("LASPH = .TRUE.")
        if m.nbands > 0:
            lines.append(f"NBANDS = {m.nbands}")
        if m.isym is not None:
            lines.append(f"ISYM = {m.isym}")
        _, order = self._ordered_atoms(spec)
        if m.ispin == 2 and m.magmom is not None:
            mm = list(m.magmom)
            lines.append("MAGMOM = " + " ".join(_inc(mm[i]) for i in order))
        elif m.ispin == 2 and m.magmom_by_element:
            sym = spec.structure.atoms.symbols
            lines.append("MAGMOM = " + " ".join(_inc(float(m.magmom_by_element.get(sym[i], 0.0))) for i in order))
        if m.hubbard:
            from adit.validate import orbital_l
            lu = [m.hubbard.get(e) for e in spec.elements]
            lines += ["LDAU = .TRUE.", f"LDAUTYPE = {m.ldau_type}",
                      "LDAUL = " + " ".join(str(orbital_l(h.orbital)) if h else "-1" for h in lu),
                      "LDAUU = " + " ".join(_inc(float(h.u_ev)) if h else "0" for h in lu),
                      "LDAUJ = " + " ".join(_inc(float(h.j_ev)) if h else "0" for h in lu)]
            if t.type == "band_structure" and m.lmaxmix <= 0:
                lmax = max(orbital_l(h.orbital) for h in m.hubbard.values())
                if lmax >= 2 and "LMAXMIX" not in {k.upper() for k in m.extra_incar}:
                    lines.append(f"LMAXMIX = {2 * lmax}")
        if m.lmaxmix > 0:
            lines.append(f"LMAXMIX = {m.lmaxmix}")
        if st.multiplicity != 1:
            lines.append(f"NUPDOWN = {st.multiplicity - 1}")
        if st.charge != 0:
            if lib is None:
                raise GenerationError(L("電荷が 0 ではないので NELECT を書く必要がありますが、POTCAR ライブラリが無いため価電子数を決められません", "a non-zero charge needs NELECT, but no POTCAR library is available to count valence electrons"))
            zsum = sum(lib.header(potcar_name(m, e)).zval * st.atoms.symbols.count(e) for e in spec.elements)
            lines.append(f"NELECT = {_inc(zsum - st.charge)}")
        if m.ivdw is not None:
            lines.append(f"IVDW = {m.ivdw}")
        if m.idipol:
            lines.append(f"IDIPOL = {m.idipol}")
        if m.ldipol:
            lines.append("LDIPOL = .TRUE.")
        if m.dipol.strip():
            lines.append(f"DIPOL = {' '.join(m.dipol.split())}")
        lines += ["", "# 構造"]
        if t.type == "single_point":
            lines += ["NSW = 0", "IBRION = -1"]
        elif t.type == "geometry_optimization":
            lines += [f"NSW = {t.max_steps}", f"IBRION = {m.ibrion}", f"ISIF = {ISIF_OF[t.relax_cell]}",
                      f"EDIFFG = -{_inc(t.force_tolerance_ev_per_ang)}"]
        elif t.type == "vibrations":
            lines += ["NSW = 1", "IBRION = 5", "NFREE = 2", "POTIM = 0.015", "ISIF = 2"]
        elif t.type == "molecular_dynamics":
            md_lines = self._md_lines(t.md)
            md_lines = [(" ".join([l.split("=")[1].strip()] * len(spec.elements)).join([l.split("=")[0] + "= ", ""]) if l.startswith("LANGEVIN_GAMMA =") else l) for l in md_lines]
            lines += md_lines
        elif t.type == "band_structure":
            lines += ["NSW = 0", "IBRION = -1"]
        if m.extra_incar:
            lines += ["", "# 追加 (extra_incar。そのまま書く)"]
            lines += [f"{k.strip().upper()} = {_inc(v)}" for k, v in m.extra_incar.items()]
        return "\n".join(lines) + "\n"

    @staticmethod
    def _md_lines(md) -> list[str]:
        # NBLOCK: write an ionic configuration to XDATCAR every NBLOCK ionic steps (wiki URL above).
        out = ["IBRION = 0", f"NSW = {md.steps}", f"NBLOCK = {md.dump_interval}", f"POTIM = {_inc(md.timestep_fs)}", f"TEBEG = {_inc(md.temperature_k)}", f"TEEND = {_inc(md.temperature_k)}"]
        if md.ensemble == "NVE":
            return out + ["SMASS = -3", "ISIF = 2"]
        if md.thermostat == "andersen":
            out += ["MDALGO = 1", f"ANDERSEN_PROB = {_inc(min(1.0, md.timestep_fs / md.coupling_time_fs))}"]
        elif md.thermostat == "langevin":
            out += ["MDALGO = 3", f"LANGEVIN_GAMMA = {_inc(1000.0 / md.coupling_time_fs)}"]
        elif md.thermostat == "nose_hoover":
            out += ["MDALGO = 2"]
        elif md.thermostat == "csvr":
            # CSVR_PERIOD counts MD steps (interpreted together with POTIM), not fs: wiki CSVR_PERIOD page.
            out += ["MDALGO = 5", f"CSVR_PERIOD = {max(1, round(md.coupling_time_fs / md.timestep_fs))}"]
        else:
            raise GenerationError(L(f"VASP に対応していない熱浴: {md.thermostat} (andersen / langevin / nose_hoover / csvr)", f"thermostat not supported for VASP: {md.thermostat} (andersen / langevin / nose_hoover / csvr)"))
        if md.ensemble == "NPT":
            if md.thermostat != "langevin":
                raise GenerationError(L("VASP の NPT は langevin 熱浴 (MDALGO = 3) と ISIF = 3 の組み合わせです (wiki の MDALGO ページ)", "VASP NPT is the combination of the langevin thermostat (MDALGO = 3) and ISIF = 3 (wiki MDALGO page)"))
            out += ["ISIF = 3", f"PSTRESS = {_inc(md.pressure_bar / 1000.0)}"]  # kbar
        else:
            out.append("ISIF = 2")
        return out

    def kpoints(self, spec: CalculationSpec) -> str:
        kp = spec.kpoints
        if kp is None or kp.mode == "gamma":
            return "Gamma-point only\n 0\nGamma\n 1 1 1\n 0 0 0\n"
        n1, n2, n3 = kp.resolved_mesh(spec.structure.atoms.cell)
        s = kp.shift
        scheme = "Gamma" if getattr(spec.method, "kpoints_centering", "monkhorst-pack") == "gamma" else "Monkhorst Pack"
        return f"adit {kp.mode}\n 0\n{scheme}\n {n1} {n2} {n3}\n {_inc(s[0])} {_inc(s[1])} {_inc(s[2])}\n"

    def potcar_spec(self, spec: CalculationSpec, lib: PotcarLibrary | None) -> str:
        m = spec.method
        lines = [f"# 使う POTCAR。順番は POSCAR の元素の順。set = {m.potcar_set}",
                 "# element  name  TITEL  ZVAL   (TITEL と ZVAL はライブラリがあれば生成時に、無ければ make_potcar.sh が実行時に potcar.used へ記録)"]
        for e in spec.elements:
            name = potcar_name(m, e)
            if lib is not None and lib.has(name):
                h = lib.header(name)
                lines.append(f"{e}  {name}  {h.titel!r}  {_inc(h.zval)}")
            else:
                lines.append(f"{e}  {name}  (未確認)  (未確認)")
        return "\n".join(lines) + "\n"

    def make_potcar_sh(self, spec: CalculationSpec) -> str:
        m = spec.method
        names = " ".join(potcar_name(m, e) for e in spec.elements)
        return f'''#!/bin/bash
# potcar.spec の順に $VASP_PP_PATH/{m.potcar_set}/<名前>/POTCAR を連結して POTCAR を作る (ADIT が生成)
# Concatenates the POTCARs listed in potcar.spec, in that order, into one POTCAR (generated by ADIT)
set -e
cd "$(dirname "$0")"
: "${{{PP_ENV}:?{PP_ENV} が空です。POTCAR ライブラリの親ディレクトリを環境変数で指定してから実行してください / {PP_ENV} is empty: set it to the parent directory of the POTCAR library first}}"
: > POTCAR
: > potcar.used
for name in {names}; do
  f="${PP_ENV}/{m.potcar_set}/$name/POTCAR"
  [ -s "$f" ] || {{ echo "POTCAR がありません / no POTCAR at: $f" >&2; exit 1; }}
  cat "$f" >> POTCAR
  printf '%s  %s  %s\\n' "$name" "$(grep -m1 TITEL "$f" | sed 's/^ *//')" "$(grep -m1 ZVAL "$f" | sed 's/^ *//')" >> potcar.used
done
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum POTCAR > potcar.sha256  # 作成時の記録: 組み立てた POTCAR の SHA-256 (中身の照合用のハッシュ。POTCAR 自体は配れないので値だけを残す)
fi
echo "POTCAR を作りました ($(grep -c TITEL POTCAR) 個。内訳は potcar.used) / wrote POTCAR ($(grep -c TITEL POTCAR) entries; see potcar.used)"
'''


register(VaspGenerator())
