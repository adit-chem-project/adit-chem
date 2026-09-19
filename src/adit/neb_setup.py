
from __future__ import annotations

from adit.errors import AditValueError
import io
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from adit import batch
from adit.lang import L
from adit.spec import AtomsData, CalculationSpec

NATIVE_CODES = ("vasp", "espresso", "cp2k")
MODES = ("native", "images")
INTERPOLATIONS = ("idpp", "linear")
QE_OPT_SCHEMES = ("sd", "broyden", "broyden2", "quick-min", "langevin")
MAX_IMAGES = 100
CELL_TOL = 1e-4
NEB_FILE = "neb.json"
IMAGES_FILE = "images.extxyz"


class NebError(AditValueError):
    pass


@dataclass
class NebOptions:
    images: int
    mode: str = "native"
    interpolation: str = "idpp"
    climb: bool = False
    mic: bool | None = None
    vasp_spring: float | None = None
    cp2k_k_spring: float | None = None
    qe_opt_scheme: str | None = None


def read_endpoint(ref: str, spec: CalculationSpec, where: Path):
    from adit.structure import StructureError, from_file

    if ref == "spec":
        return spec.atoms
    p = Path(ref).expanduser()
    p = p if p.is_absolute() else Path(where) / p
    try:
        return from_file(p)
    except StructureError as ex:
        raise NebError(str(ex)) from ex


def check_endpoints(a, b) -> list[str]:
    errs = []
    if len(a) != len(b):
        return [L(f"原子の数が違います (始状態 {len(a)}、終状態 {len(b)})", f"different numbers of atoms (initial {len(a)}, final {len(b)})")]
    sa, sb = a.get_chemical_symbols(), b.get_chemical_symbols()
    bad = next((i for i, (x, y) in enumerate(zip(sa, sb)) if x != y), None)
    if bad is not None:
        errs.append(L(f"原子の並びが違います: {bad + 1} 番目が始状態 {sa[bad]}、終状態 {sb[bad]} (NEB は同じ番号の原子どうしを結びます)",
                      f"atom order differs: atom {bad + 1} is {sa[bad]} in the initial and {sb[bad]} in the final state (NEB connects atoms with the same index)"))
    if list(a.pbc) != list(b.pbc):
        errs.append(L(f"周期が違います (始状態 {list(a.pbc)}、終状態 {list(b.pbc)})", f"periodicity differs (initial {list(a.pbc)}, final {list(b.pbc)})"))
    elif any(a.pbc):
        d = float(np.abs(np.asarray(a.cell) - np.asarray(b.cell)).max())
        if d > CELL_TOL:
            errs.append(L(f"セルが違います (成分の差の最大 {d:.4g} Å)。NEB の像はすべて同じセルです (xyz にはセルが無いので、extxyz・POSCAR・cif などを使います)",
                          f"the cells differ (max component difference {d:.4g} Å); all NEB images share one cell (xyz has no cell; use extxyz, POSCAR, cif, ...)"))
    return errs


def interpolate(a, b, n: int, method: str = "idpp", mic: bool = False) -> list:
    from ase.mep import NEB

    images = [a.copy()] + [a.copy() for _ in range(n)] + [b.copy()]
    NEB(images, method="improvedtangent").interpolate(method=method, mic=mic)
    return images


def _spec_for(spec: CalculationSpec, atoms, task_type: str | None = None, comment: str = "") -> CalculationSpec:
    data = spec.model_dump(mode="json")
    st = dict(data["structure"])
    ad = AtomsData.from_ase(atoms)
    if not any(atoms.pbc) and not any(spec.structure.atoms.pbc):
        ad = ad.model_copy(update={"cell": spec.structure.atoms.cell})
    st.update(atoms=ad.model_dump(mode="json"), velocities=None)
    data.update(structure=st, handoff=None)
    if task_type:
        data["task"] = dict(data["task"], type=task_type)
    if comment:
        data["meta"] = dict(data["meta"], comment=(data["meta"].get("comment", "") + " " + comment).strip())
    return CalculationSpec.model_validate(data)


def plan(spec: CalculationSpec, start_ref: str, end_ref: str, opts: NebOptions, where: Path | str = ".") -> tuple[list, list[str]]:
    from adit.validate import MIN_DISTANCE, close_pairs_atoms

    if opts.mode not in MODES:
        raise NebError(L(f"NEB の作り方は {' / '.join(MODES)} のどちらかです", f"NEB mode must be {' / '.join(MODES)}"))
    if opts.interpolation not in INTERPOLATIONS:
        raise NebError(L(f"補間は {' / '.join(INTERPOLATIONS)} のどちらかです", f"interpolation must be {' / '.join(INTERPOLATIONS)}"))
    if not 1 <= opts.images <= MAX_IMAGES:
        raise NebError(L(f"中間の像の数は 1〜{MAX_IMAGES} にしてください", f"the number of intermediate images must be 1..{MAX_IMAGES}"))
    code = spec.method.code
    if opts.mode == "native" and code not in NATIVE_CODES:
        raise NebError(L(f"計算コードの NEB を書けるのは VASP / QE / CP2K です ({code} は --neb-mode images で、像ごとの一点計算にできます)",
                         f"code-native NEB is written for VASP / QE / CP2K ({code}: use --neb-mode images for one single point per image)"))
    if opts.mode == "native" and code == "vasp" and opts.climb:
        raise NebError(L("VASP の climbing image (LCLIMB) は VTST のタグで、VASP 本体の文書にありません。ADIT は書きません (--climb を外してください)",
                         "the VASP climbing image (LCLIMB) is a VTST tag, not in the VASP documentation; ADIT does not write it (drop --climb)"))
    if opts.qe_opt_scheme is not None and opts.qe_opt_scheme not in QE_OPT_SCHEMES:
        raise NebError(L(f"opt_scheme は {' / '.join(QE_OPT_SCHEMES)} のどれかです (INPUT_NEB)", f"opt_scheme must be one of {' / '.join(QE_OPT_SCHEMES)} (INPUT_NEB)"))
    if opts.mode == "native" and spec.task.relax_cell != "no":
        raise NebError(L("NEB の像はすべて同じセルなので、格子を動かす指定 (task.relax_cell) は使えません", "all NEB images share the cell, so cell relaxation (task.relax_cell) cannot be used"))
    a, b = read_endpoint(start_ref, spec, Path(where)), read_endpoint(end_ref, spec, Path(where))
    if not any(a.pbc) and not any(b.pbc) and spec.structure.periodic:
        a.cell, a.pbc = spec.atoms.cell, spec.atoms.pbc
        b.cell, b.pbc = spec.atoms.cell, spec.atoms.pbc
    errs = check_endpoints(a, b)
    if not errs and len(a) != len(spec.structure.atoms.symbols):
        errs.append(L(f"spec.json の構造 ({len(spec.structure.atoms.symbols)} 原子) と端の構造 ({len(a)} 原子) の原子数が違います (電荷・固定原子は spec.json のものを使うため)",
                      f"spec.json has {len(spec.structure.atoms.symbols)} atoms but the end points have {len(a)} (charge and fixed atoms come from spec.json)"))
    if errs:
        raise NebError("\n".join(errs))
    mic = opts.mic if opts.mic is not None else bool(any(a.pbc))
    images = interpolate(a, b, opts.images, opts.interpolation, mic)
    for i, img in enumerate(images[1:-1], 1):
        close = close_pairs_atoms(img, MIN_DISTANCE)
        if close:
            x, y, d = close[0]
            raise NebError(L(f"像 {i} で原子 {x + 1} と {y + 1} の距離が {d:.3f} Å です ({MIN_DISTANCE} Å 未満。ほかに {len(close) - 1} 組)。原子の並びか端の構造を確かめてください",
                             f"in image {i}, atoms {x + 1} and {y + 1} are {d:.3f} Å apart (below {MIN_DISTANCE} Å; {len(close) - 1} more pairs); check the atom order or the end points"))
    notes = []
    fixed = sorted(spec.structure.fixed_atoms)
    if fixed:
        from ase.geometry import find_mic
        dr = b.positions[fixed] - a.positions[fixed]
        if any(a.pbc):
            dr, _ = find_mic(dr, a.cell, a.pbc)
        moved = [(i, float(np.linalg.norm(v))) for i, v in zip(fixed, dr) if np.linalg.norm(v) > 1e-3]
        if moved:
            notes.append(L(f"固定原子のうち {len(moved)} 個が始状態と終状態で動いています (最大 {max(d for _, d in moved):.3f} Å)。像ごとの位置は補間した位置のまま固定されます",
                           f"{len(moved)} fixed atoms move between the end points (max {max(d for _, d in moved):.3f} Å); in each image they stay fixed at the interpolated position"))
    notes.append(L(f"補間: ASE の NEB.interpolate (method = {opts.interpolation}、mic = {mic})", f"interpolation: ASE NEB.interpolate (method = {opts.interpolation}, mic = {mic})"))
    return images, notes


def _xyz(atoms, fmt: str = "xyz") -> str:
    from ase.io import write

    buf = io.StringIO()
    write(buf, atoms, format=fmt)
    return buf.getvalue()


def _images_text(images) -> str:
    from ase.io import write

    buf = io.StringIO()
    write(buf, images, format="extxyz")
    return buf.getvalue()


def _record(spec, opts, images, notes, start_ref, end_ref, extra=None) -> str:
    from adit import __version__

    d = {"generated_by": f"adit {__version__}", "code": spec.method.code, "options": asdict(opts), "n_images_total": len(images),
         "start": start_ref, "end": end_ref, "notes": notes, "images_file": IMAGES_FILE,
         "path_length_ang": _path_length(images)}
    d.update(extra or {})
    return json.dumps(d, ensure_ascii=False, indent=2) + "\n"


def _path_length(images) -> list[float]:
    from ase.geometry import find_mic

    s = [0.0]
    for a, b in zip(images[:-1], images[1:]):
        dr = b.positions - a.positions
        if any(a.pbc):
            dr, _ = find_mic(dr, a.cell, a.pbc)
        s.append(s[-1] + float(np.sqrt((dr ** 2).sum())))
    return s


def _readme_common(opts, images, notes) -> list[str]:
    climb_note = [] if opts.climb else [L(
        "  - climbing image は使いません。最高エネルギーの像は鞍点そのものとは限らないため、像の最大値から得る上昇量は、遷移状態を探して求めた障壁ではありません",
        "  - No climbing image is used. The highest-energy image need not coincide with the saddle point, so the rise to the highest image is not a barrier from a located transition state")]
    return [L("== 反応経路 (NEB) ==", "== Reaction path (NEB) =="),
            L(f"  始状態・終状態を含めて {len(images)} 個の像 (中間 {opts.images} 個)。全像は images.extxyz、条件は neb.json",
              f"  {len(images)} images including the end points ({opts.images} intermediate); all images in images.extxyz, settings in neb.json"),
            *[f"  - {n}" for n in notes], *climb_note]


def _validate_or_raise(s: CalculationSpec, cfg) -> None:
    from adit.project import ProjectError
    from adit.validate import validate

    errs = validate(s, cfg)
    if errs:
        raise ProjectError(L("次の点を直すと生成できます:", "fix the following to generate:") + "\n" + "\n".join(f"  - {e}" for e in errs), errs)


def _vasp(spec, cfg, out: Path, images, opts, notes, rec, overwrite) -> list[Path]:
    from adit.codes.vasp import VaspGenerator
    from adit.project import build_project, write_project

    if len(images) > 100:
        raise NebError(L(f"VASP の NEB の像のディレクトリは 00〜99 の 2 桁なので、始状態と終状態を含めて 100 個までです (いまは {len(images)} 個。中間の像を {100 - 2} 個までにしてください)",
                         f"VASP NEB image directories are two digits (00..99), so at most 100 images including the end points (now {len(images)}; use at most {100 - 2} intermediate images)"))
    s0 = _spec_for(spec, images[0], "geometry_optimization", f"neb {len(images)} images")
    extra = dict(s0.method.extra_incar)
    clash = [k for k in extra if k.strip().upper() in ("IMAGES", "SPRING", "LCLIMB")]
    if clash:
        raise NebError(L(f"追加の INCAR に {clash} があります。NEB の生成が書くので外してください", f"the extra INCAR has {clash}; the NEB generator writes these, so remove them"))
    extra["IMAGES"] = opts.images
    if opts.vasp_spring is not None:
        extra["SPRING"] = float(opts.vasp_spring)
    s0 = s0.model_copy(update={"method": s0.method.model_copy(update={"extra_incar": extra})})
    gen = VaspGenerator()
    texts = {f"{i:02d}/POSCAR": gen.poscar(_spec_for(spec, img)) for i, img in enumerate(images)}
    texts.update({NEB_FILE: rec, IMAGES_FILE: _images_text(images)})
    last = f"{len(images) - 1:02d}"
    lines = _readme_common(opts, images, notes) + [
        L(f"  00 … {last}/POSCAR  各像の構造。VASP は IMAGES = {opts.images} で 01 … {len(images) - 2:02d} を計算します (00 と {last} は計算しません)",
          f"  00 ... {last}/POSCAR  structure of each image; with IMAGES = {opts.images} VASP computes 01 ... {len(images) - 2:02d} (not 00 and {last})"),
        L("  INCAR の IMAGES は、追加の INCAR (extra_incar) の欄として書いてあります。SPRING は指定したときだけ書きます (VASP の既定 -5 は NEB)",
          "  IMAGES is written as an extra INCAR entry; SPRING only when given (the VASP default -5 means NEB)"),
        L("  最適化の方法 (IBRION) の選び方は VASP wiki の Nudged elastic bands の頁にあります。climbing image (LCLIMB) は VTST のタグなので書いていません",
          "  How to choose the optimizer (IBRION) is on the VASP wiki page Nudged elastic bands; the climbing image (LCLIMB) is a VTST tag and is not written"),
        L(f"  端の像のエネルギー: endpoints/initial と endpoints/final を実行し、それぞれの OUTCAR を 00/ と {last}/ に写すと、解析 (adit-analyze) が読めます",
          f"  End-point energies: run endpoints/initial and endpoints/final and copy each OUTCAR into 00/ and {last}/ so that the analysis (adit-analyze) can read them"),
        f"       cp endpoints/initial/OUTCAR 00/ && cp endpoints/final/OUTCAR {last}/"]
    ends = [(out / "endpoints" / tag, _spec_for(spec, img, "single_point", f"neb end point {tag}")) for tag, img in (("initial", images[0]), ("final", images[-1]))]
    batch.check_output(out, overwrite)
    out.mkdir(parents=True, exist_ok=True)
    build_project(s0, cfg, output_dir=out, extra_texts=texts, drop=("POSCAR",), extra_readme=lines)
    for d, s in ends:
        build_project(s, cfg)
    write_project(s0, cfg, out, overwrite=True, extra_texts=texts, drop=("POSCAR",), extra_readme=lines)
    for d, s in ends:
        d.parent.mkdir(parents=True, exist_ok=True)
        write_project(s, cfg, d, overwrite=overwrite)
    return [out] + [d for d, _ in ends]


def _pw_chunks(text: str) -> tuple[list[list[str]], dict[str, list[str]]]:
    lines = text.splitlines()
    namelists, cards, i = [], {}, 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("&"):
            j = i
            while lines[j].strip() != "/":
                j += 1
            namelists.append(lines[i:j + 1])
            i = j + 1
        elif s:
            j = i
            while j < len(lines) and lines[j].strip():
                j += 1
            cards[s.split()[0].upper()] = lines[i:j]
            i = j
        else:
            i += 1
    return namelists, cards


def neb_in(pw_texts: list[str], opts: NebOptions, nstep: int) -> str:
    namelists, cards = _pw_chunks(pw_texts[0])
    drop_ctrl = ("calculation", "tprnfor", "tstress")
    engine = []
    for nl in namelists:
        name = nl[0].strip().upper()
        body = nl[1:-1]
        if name == "&CONTROL":
            body = [l for l in body if l.split("=")[0].strip().lower() not in drop_ctrl]
        if name in ("&CELL", "&FCP", "&RISM") and not body:
            continue
        engine += [nl[0], *body, "/"]
    if "ATOMIC_POSITIONS" not in cards or "ATOMIC_SPECIES" not in cards:
        raise NebError(L("pw.in を neb.in に組み替えられません (ATOMIC_SPECIES か ATOMIC_POSITIONS が見つからない)", "cannot rebuild pw.in as neb.in (ATOMIC_SPECIES or ATOMIC_POSITIONS missing)"))
    engine += cards["ATOMIC_SPECIES"]
    engine.append("BEGIN_POSITIONS")
    last = len(pw_texts) - 1
    for i, t in enumerate(pw_texts):
        _, c = _pw_chunks(t)
        engine.append("FIRST_IMAGE" if i == 0 else ("LAST_IMAGE" if i == last else "INTERMEDIATE_IMAGE"))
        engine += c["ATOMIC_POSITIONS"]
    engine.append("END_POSITIONS")
    for name, card in cards.items():
        if name not in ("ATOMIC_SPECIES", "ATOMIC_POSITIONS"):
            engine += card
    path = ["&PATH", "   string_method = 'neb'", f"   nstep_path = {nstep}", f"   num_of_images = {len(pw_texts)}"]
    if opts.climb:
        path.append("   CI_scheme = 'auto'")
    if opts.qe_opt_scheme:
        path.append(f"   opt_scheme = '{opts.qe_opt_scheme}'")
    path.append("/")
    return "\n".join(["BEGIN", "BEGIN_PATH_INPUT", *path, "END_PATH_INPUT", "BEGIN_ENGINE_INPUT", *engine, "END_ENGINE_INPUT", "END"]) + "\n"


def _espresso(spec, cfg, out: Path, images, opts, notes, rec, overwrite) -> list[Path]:
    from adit.codes.espresso import INPUT_FILE, EspressoGenerator, profile_cmd
    from adit.project import write_project

    if "pw.x" not in profile_cmd(cfg, spec):
        raise NebError(L("neb.x は pw.x と同じ形で呼ぶので、commands.espresso に 'pw.x' の文字が要ります (既定のままで使えます)",
                         "neb.x is called in the same form as pw.x, so commands.espresso must contain 'pw.x' (the default is fine)"))
    s0 = _spec_for(spec, images[0], "single_point", f"neb {len(images)} images")
    _validate_or_raise(s0, cfg)
    gen = EspressoGenerator()
    lib = gen.resolve(s0, cfg)
    nstep = spec.task.max_steps if spec.task.max_steps > 0 else 1
    text = neb_in([gen.pw_in(_spec_for(spec, img, "single_point"), lib) for img in images], opts, nstep)

    def run(cmd: str) -> str:
        new = cmd.replace("pw.x", "neb.x", 1).replace(f"-in {INPUT_FILE}", "-inp neb.in")
        if "-inp neb.in" not in new:
            raise NebError(L(f"実行コマンドを neb.x の形に直せません: {cmd}", f"cannot turn the run command into neb.x form: {cmd}"))
        return new

    lines = _readme_common(opts, images, notes) + [
        L("  neb.in        neb.x の入力 (BEGIN_PATH_INPUT の &PATH と、pw.x と同じ条件の BEGIN_ENGINE_INPUT。全像の座標を書いてあります)",
          "  neb.in        neb.x input (&PATH in BEGIN_PATH_INPUT and a BEGIN_ENGINE_INPUT with the pw.x settings; coordinates of every image are given)"),
        L(f"  nstep_path = {nstep} (計算の設定の最大ステップ数)、num_of_images = {len(images)}"
          + ("、CI_scheme = 'auto'" if opts.climb else "") + (f"、opt_scheme = '{opts.qe_opt_scheme}'" if opts.qe_opt_scheme else "")
          + "。書いていない変数は neb.x の既定です (INPUT_NEB)",
          f"  nstep_path = {nstep} (max steps of the task), num_of_images = {len(images)}"
          + (", CI_scheme = 'auto'" if opts.climb else "") + (f", opt_scheme = '{opts.qe_opt_scheme}'" if opts.qe_opt_scheme else "")
          + "; other variables are neb.x defaults (INPUT_NEB)"),
        L("  結果: adit.dat (像ごとのエネルギー)、adit.int (補間)、adit.path、adit.xyz。解析 (adit-analyze) は adit.dat と adit.int を読みます",
          "  Results: adit.dat (energy per image), adit.int (interpolation), adit.path, adit.xyz; the analysis (adit-analyze) reads adit.dat and adit.int")]
    write_project(s0, cfg, out, overwrite=overwrite, extra_texts={"neb.in": text, NEB_FILE: rec, IMAGES_FILE: _images_text(images)},
                  drop=(INPUT_FILE,), run_transform=run, extra_readme=lines)
    return [out]


def _cp2k(spec, cfg, out: Path, images, opts, notes, rec, overwrite) -> list[Path]:
    from adit.codes.cp2k import INPUT_FILE, Cp2kGenerator
    from adit.project import write_project
    from adit.spec import HARTREE_PER_BOHR_IN_EV_PER_ANG

    if not all(images[0].pbc):
        raise NebError(L("CP2K の BAND は、3 方向とも周期の構造 (分子なら箱に入れて周期にしたもの) だけ生成します (非周期の座標の中心合わせを確かめていないため)",
                         "CP2K BAND is generated only for structures periodic in all directions (put a molecule in a periodic box); centering of non-periodic replicas was not verified"))
    s0 = _spec_for(spec, images[0], "single_point", f"neb {len(images)} images")
    _validate_or_raise(s0, cfg)
    gen = Cp2kGenerator()
    data = gen.resolve(s0, cfg)
    text = gen.cp2k_inp(s0, data)
    if text.count("  RUN_TYPE ENERGY_FORCE\n") != 1 or "&MOTION" in text:
        raise NebError(L("cp2k.inp を BAND の形に直せません", "cannot turn cp2k.inp into a BAND input"))
    text = text.replace("  RUN_TYPE ENERGY_FORCE\n", "  RUN_TYPE BAND\n")
    band = ["&MOTION", "  &BAND"]
    if opts.climb:
        band.append("    BAND_TYPE CI-NEB")
    band.append(f"    NUMBER_OF_REPLICA {len(images)}")
    if opts.cp2k_k_spring is not None:
        band.append(f"    K_SPRING {opts.cp2k_k_spring:.10g}")
    band += ["    &CONVERGENCE_CONTROL", f"      MAX_FORCE {spec.task.force_tolerance_ev_per_ang / HARTREE_PER_BOHR_IN_EV_PER_ANG:.6e}", "    &END CONVERGENCE_CONTROL",
             "    &OPTIMIZE_BAND", "      OPT_TYPE DIIS", "      &DIIS", f"        MAX_STEPS {max(1, spec.task.max_steps)}", "      &END DIIS", "    &END OPTIMIZE_BAND"]
    replicas = {}
    for i, img in enumerate(images):
        name = f"replica_{i:02d}.xyz"
        replicas[name] = _xyz(img)
        band += ["    &REPLICA", f"      COORD_FILE_NAME {name}", "    &END REPLICA"]
    band.append("  &END BAND")
    fixed = gen._fixed_blocks(s0)
    if fixed:
        band += ["  &CONSTRAINT", *fixed, "  &END CONSTRAINT"]
    band.append("&END MOTION")
    text = text + "\n".join(band) + "\n"
    lines = _readme_common(opts, images, notes) + [
        L(f"  cp2k.inp の RUN_TYPE BAND と &MOTION/&BAND。BAND_TYPE {'CI-NEB' if opts.climb else '(書かない。CP2K の既定 IT-NEB)'}、NUMBER_OF_REPLICA {len(images)}、"
          f"MAX_FORCE は計算の設定の力の収束基準、MAX_STEPS は最大ステップ数。ALIGN_FRAMES・ROTATE_FRAMES などは CP2K の既定のままです",
          f"  RUN_TYPE BAND and &MOTION/&BAND in cp2k.inp: BAND_TYPE {'CI-NEB' if opts.climb else '(not written; CP2K default IT-NEB)'}, NUMBER_OF_REPLICA {len(images)}, "
          f"MAX_FORCE from the force criterion, MAX_STEPS from the max steps; ALIGN_FRAMES, ROTATE_FRAMES etc. are CP2K defaults"),
        L("  replica_00.xyz …  各像の座標 (&REPLICA の COORD_FILE_NAME)", "  replica_00.xyz ...  coordinates of each image (COORD_FILE_NAME of &REPLICA)"),
        L("  解析 (adit-analyze) は output.log の像ごとのエネルギーと像間の距離を読みます。CP2K の出力には補間曲線がないため、像の値だけを表示します",
          "  Analysis (adit-analyze) reads the image energies and inter-image distances from output.log. CP2K does not output an interpolated curve, so only the image values are shown")]
    write_project(s0, cfg, out, overwrite=overwrite, extra_texts={INPUT_FILE: text, NEB_FILE: rec, IMAGES_FILE: _images_text(images), **replicas},
                  extra_readme=lines)
    return [out]


def _images_mode(spec, cfg, out: Path, images, opts, notes, rec, overwrite) -> list[Path]:
    from adit import __version__

    items = []
    for i, img in enumerate(images):
        items.append(batch.Item(f"image_{i:02d}", _spec_for(spec, img, "single_point", f"neb image {i}/{len(images) - 1}"), [
            L("== 反応経路の像 ==", "== Reaction-path image =="),
            L(f"  {len(images)} 個の像のうち {i} 番目 (0 が始状態)。補間しただけの構造での一点計算です (NEB の最適化はしていません)。全体は ../README.txt",
              f"  image {i} of {len(images)} (0 = initial state); a single point on the interpolated structure (no NEB optimization). See ../README.txt")]))
    dirs = batch.write_items(items, cfg, out, overwrite=overwrite)
    (out / NEB_FILE).write_text(rec, encoding="utf-8")
    (out / IMAGES_FILE).write_text(_images_text(images), encoding="utf-8")
    batch.write_scan_json(out, "image", [f"{i:02d}" for i in range(len(images))], dirs)
    lines = [L(f"ADIT {__version__} が生成した、補間した反応経路の像ごとの一点計算です ({spec.method.code})", f"Single points on interpolated path images, generated by ADIT {__version__} ({spec.method.code})"), "",
             *_readme_common(opts, images, notes),
             L("  これは NEB の最適化ではありません (像を補間した位置のまま、エネルギーだけを計算します)。", "  This is not an NEB optimization (only energies at the interpolated positions)."),
             L("  経路に沿った長さ [Å] は neb.json の path_length_ang。", "  The path length [Å] is path_length_ang in neb.json."),
             "", L("== 実行したあと ==", "== After running =="),
             L("  adit-analyze <このディレクトリ> --scan   (像ごとのエネルギーの表)", "  adit-analyze <this directory> --scan   (energy per image)")]
    batch.write_top(out, cfg, spec, dirs, lines, "像の計算", "image runs")
    return dirs


def write_neb(spec: CalculationSpec, cfg, out_dir: Path | str, start_ref: str, end_ref: str, opts: NebOptions, *, where: Path | str = ".",
              overwrite: bool = False) -> list[Path]:
    out = Path(out_dir).expanduser()
    images, notes = plan(spec, start_ref, end_ref, opts, where)
    rec = _record(spec, opts, images, notes, start_ref, end_ref)
    if opts.mode == "images":
        return _images_mode(spec, cfg, out, images, opts, notes, rec, overwrite)
    return {"vasp": _vasp, "espresso": _espresso, "cp2k": _cp2k}[spec.method.code](spec, cfg, out, images, opts, notes, rec, overwrite)


__all__ = ["NebError", "NebOptions", "check_endpoints", "interpolate", "plan", "neb_in", "write_neb", "NATIVE_CODES"]
