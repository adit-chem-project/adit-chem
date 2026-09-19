
from __future__ import annotations

from adit.errors import AditValueError
import json
import re
import tempfile
from pathlib import Path

import numpy as np

from adit import batch
from adit.lang import L
from adit.spec import AtomsData, CalculationSpec

CODES = ("dftbplus", "vasp", "espresso")
BACKENDS = ("auto", "phonopy", "ase")
PHONON_FILE = "phonons.json"
COLLECT_SCRIPT = "phonon_collect.py"
BAND_POINTS = 51
DOS_GRID_POINTS = 400
MAX_SUPERCELL_ATOMS = 2000
MAX_DOS_QPOINTS_X_BANDS = 200_000


class PhononError(AditValueError):
    pass


def _thz_per_ev() -> float:
    from ase import units
    return units._e / units._hplanck / 1e12   # 1 eV = 241.8 THz (E = h ν)


def has_phonopy() -> bool:
    try:
        import phonopy  # noqa: F401
        return True
    except ImportError:
        return False


def parse_dim(text: str) -> tuple[int, int, int]:
    parts = [p for p in re.split(r"[x×,\s]+", str(text).strip()) if p]
    try:
        dim = tuple(int(p) for p in parts)
    except ValueError as ex:
        raise PhononError(L(f"3 つの整数を 2x2x2 の形で書いてください: {text!r}", f"write three integers as 2x2x2: {text!r}")) from ex
    if len(dim) != 3 or any(k < 1 for k in dim):
        raise PhononError(L(f"1 以上の整数を 3 つ書いてください (例 2x2x2): {text!r}", f"give three positive integers (e.g. 2x2x2): {text!r}"))
    return dim


def band_paths(cell, npoints: int = BAND_POINTS):
    from ase.cell import Cell
    from ase.dft.kpoints import parse_path_string

    bp = Cell(np.asarray(cell)).bandpath(npoints=0)
    sp = bp.special_points
    paths, labels, conns = [], [], []
    for piece in parse_path_string(bp.path):
        for a, b in zip(piece[:-1], piece[1:]):
            paths.append(np.array([sp[a] + (sp[b] - sp[a]) * t for t in np.linspace(0, 1, npoints)]))
        conns += [True] * (len(piece) - 2) + [False]
        labels += list(piece)
    return paths, labels, conns, bp.path


def _unit_atoms(spec: CalculationSpec):
    a = spec.atoms
    if not all(a.pbc):
        raise PhononError(L("フォノン (超格子) は 3 方向とも周期の構造で作ります", "phonons (supercells) need a structure periodic in all three directions"))
    return a


def _check_dos(unit, backend: str, dos_mesh, dos_width_thz) -> None:
    if dos_width_thz is not None and (not np.isfinite(dos_width_thz) or dos_width_thz <= 0):
        raise PhononError(L(f"DOS を広げる幅は 0 より大きい数にしてください: {dos_width_thz}", f"the DOS width must be a number greater than 0: {dos_width_thz}"))
    if dos_mesh is None:
        if dos_width_thz is not None:
            raise PhononError(L("DOS の幅を使うには、q 点のメッシュ (--phonon-dos-mesh) も要ります", "a DOS width needs a q-point mesh too (--phonon-dos-mesh)"))
        return
    if len(dos_mesh) != 3 or any(int(k) < 1 for k in dos_mesh):
        raise PhononError(L("DOS の q 点のメッシュは 1 以上の整数を 3 つ書いてください", "the DOS q-point mesh needs three positive integers"))
    if backend == "phonopy" and dos_width_thz is not None:
        raise PhononError(L("phonopy の DOS はテトラヘドロン法 (phonopy の既定) で、幅を使いません。--phonon-dos-width を外してください",
                            "the phonopy DOS uses the tetrahedron method (phonopy default) without a width; drop --phonon-dos-width"))
    if backend == "ase" and dos_width_thz is None:
        raise PhononError(L("ASE の DOS には広げる幅が要ります (--phonon-dos-width [THz]。ADIT は値を持ちません)",
                            "the ASE DOS needs a broadening width (--phonon-dos-width in THz; ADIT has no default)"))
    n = int(np.prod([int(k) for k in dos_mesh])) * 3 * len(unit)
    if n > MAX_DOS_QPOINTS_X_BANDS:
        raise PhononError(L(f"DOS の q 点 × バンドの数が {n} になり、上限 {MAX_DOS_QPOINTS_X_BANDS} を超えます (メモリが大きくなりすぎるため。"
                            f"phonopy 4.5.0 で 30x30x30 × 6 バンドが 627 MB でした)。メッシュを小さくしてください",
                            f"{n} DOS q-points x bands exceed the limit {MAX_DOS_QPOINTS_X_BANDS} (memory; 30x30x30 x 6 bands took 627 MB with phonopy 4.5.0); use a smaller mesh"))


def plan(spec: CalculationSpec, dim, *, distance: float | None = None, backend: str = "auto"):
    from ase import Atoms

    if spec.method.code not in CODES:
        raise PhononError(L(f"フォノンの一括生成は {' / '.join(CODES)} だけです (力を出力から読めるコード)", f"phonon generation supports {' / '.join(CODES)} only (codes whose forces can be read)"))
    if backend not in BACKENDS:
        raise PhononError(L(f"作り方は {' / '.join(BACKENDS)} のどれかです", f"backend must be one of {' / '.join(BACKENDS)}"))
    if distance is not None and (not np.isfinite(distance) or distance <= 0):
        raise PhononError(L(f"変位の大きさは 0 より大きい数にしてください: {distance}", f"the displacement must be a number greater than 0: {distance}"))
    if spec.structure.fixed_atoms or spec.structure.fixed_axes:
        raise PhononError(L("固定原子があると、変位の計算で力の一部が 0 にされます。フォノンの生成では固定原子を外してください", "fixed atoms would zero some forces in the displaced runs; remove them for phonons"))
    unit = _unit_atoms(spec)
    if len(unit) * int(np.prod(dim)) > MAX_SUPERCELL_ATOMS:
        raise PhononError(L(f"超格子の原子数が {len(unit) * int(np.prod(dim))} になり、上限 {MAX_SUPERCELL_ATOMS} を超えます", f"the supercell would have {len(unit) * int(np.prod(dim))} atoms, above the limit {MAX_SUPERCELL_ATOMS}"))
    use = backend
    if backend == "auto":
        use = "phonopy" if has_phonopy() else "ase"
    elif backend == "phonopy" and not has_phonopy():
        raise PhononError(L("phonopy が入っていません (pip install phonopy。または --phonon-backend ase)", "phonopy is not installed (pip install phonopy, or --phonon-backend ase)"))
    rec = {"backend": use, "dim": list(dim), "code": spec.method.code, "unit_cell": AtomsData.from_ase(unit).model_dump(mode="json"), "band_points": BAND_POINTS}
    if use == "phonopy":
        import phonopy
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms

        u = PhonopyAtoms(symbols=unit.get_chemical_symbols(), cell=np.asarray(unit.cell), scaled_positions=unit.get_scaled_positions())
        ph = Phonopy(u, supercell_matrix=np.diag(dim))
        if distance is None:
            ph.generate_displacements()
        else:
            ph.generate_displacements(distance=distance)
        scs = [Atoms(symbols=s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True) for s in ph.supercells_with_displacements]
        disp = np.asarray(ph.displacements)
        rec.update(phonopy_version=phonopy.__version__, n_displacements=len(scs),
                   distance_ang=float(np.linalg.norm(disp[0][1:4])) if len(disp) else distance,
                   primitive_matrix=np.asarray(ph.primitive_matrix).tolist() if ph.primitive_matrix is not None else None)
        return scs, rec, ph
    from ase.phonons import Phonons

    delta = 0.01 if distance is None else float(distance)
    ph = Phonons(unit, supercell=tuple(dim), delta=delta)
    big = unit * tuple(dim)
    off = len(unit) * ph.offset
    scs, keys = [], []
    for a in range(len(unit)):
        for i, v in enumerate("xyz"):
            for sign, s in ((-1, "-"), (1, "+")):
                x = big.copy()
                x.positions[off + a, i] += sign * delta
                scs.append(x)
                keys.append(f"{a}{v}{s}")
    rec.update(n_displacements=len(scs), distance_ang=delta, ase_keys=keys, ase_offset=int(ph.offset))
    return scs, rec, None


def write_phonons(spec: CalculationSpec, cfg, out_dir: Path | str, dim, *, distance: float | None = None, backend: str = "auto",
                  dos_mesh=None, dos_width_thz: float | None = None, overwrite: bool = False) -> list[Path]:
    from adit import __version__

    out = Path(out_dir).expanduser()
    scs, rec, ph = plan(spec, dim, distance=distance, backend=backend)
    _check_dos(spec.atoms, rec["backend"], dos_mesh, dos_width_thz)
    rec.update(dos_mesh=[int(k) for k in dos_mesh] if dos_mesh is not None else None, dos_width_thz=dos_width_thz)
    base = spec.model_dump(mode="json")
    base["handoff"] = None
    if base["task"]["type"] != "single_point":
        base["task"] = dict(base["task"], type="single_point")
    items = []
    for k, sc in enumerate(scs, 1):
        data = dict(base, structure=dict(base["structure"], atoms=AtomsData.from_ase(sc).model_dump(mode="json"), velocities=None, fixed_atoms=[], fixed_axes={}))
        data["meta"] = dict(base["meta"], comment=(base["meta"].get("comment", "") + f" phonon displacement {k}/{len(scs)} ({rec['backend']})").strip())
        items.append(batch.Item(f"disp-{k:03d}", CalculationSpec.model_validate(data), [
            L("== フォノンの変位 ==", "== Phonon displacement =="),
            L(f"  {len(scs)} 個の変位のうち {k} 番目 ({rec['backend']}、超格子 {'x'.join(map(str, dim))}、{rec['distance_ang']:g} Å)。力だけを使います。全体は ../README.txt",
              f"  displacement {k} of {len(scs)} ({rec['backend']}, supercell {'x'.join(map(str, dim))}, {rec['distance_ang']:g} Å); only the forces are used. See ../README.txt")]))
    rec["dirs"] = [it.dir for it in items]
    dirs = batch.write_items(items, cfg, out, overwrite=overwrite)
    if ph is not None:
        ph.save(str(out / "phonopy_disp.yaml"))
    batch.write_json(out / PHONON_FILE, {"generated_by": f"adit {__version__}", **rec})
    (out / COLLECT_SCRIPT).write_text(
        "#!/usr/bin/env python3\n" + L("# ADIT が生成。各変位の力を集めて、フォノン分散 (band.yaml) と (メッシュを与えたなら) 状態密度 (total_dos.dat) を書く。ADIT が入った Python で実行する\n",
                                        "# generated by ADIT. Collects the forces of each displacement and writes band.yaml (and total_dos.dat if a mesh was given); run with the Python that has ADIT installed\n")
        + "from pathlib import Path\nfrom adit.phonon_setup import collect\n\nprint(collect(Path(__file__).resolve().parent)['summary'])\n", encoding="utf-8", newline="\n")
    lines = [L(f"ADIT {__version__} が生成した、フォノン (有限変位) の計算です ({spec.method.code}、{rec['backend']})",
               f"Phonon (finite displacement) runs generated by ADIT {__version__} ({spec.method.code}, {rec['backend']})"), "",
             L(f"  超格子 {'x'.join(map(str, dim))} ({len(scs[0])} 原子)、変位 {len(scs)} 個、変位の大きさ {rec['distance_ang']:g} Å",
               f"  supercell {'x'.join(map(str, dim))} ({len(scs[0])} atoms), {len(scs)} displacements of {rec['distance_ang']:g} Å")]
    if spec.kpoints is not None and spec.kpoints.mode == "mesh":
        lines.append(L(f"  k 点は分割数の指定 ({'x'.join(map(str, spec.kpoints.mesh))}) なので、同じ分割数を超格子にも使っています (k 点の間隔は超格子の倍率だけ細かくなります)。",
                       f"  k-points are given as a mesh ({'x'.join(map(str, spec.kpoints.mesh))}), so the same mesh is used for the supercell (the spacing becomes finer by the supercell factor)."))
    elif spec.kpoints is not None and spec.kpoints.mode == "density":
        sc_mesh = spec.kpoints.resolved_mesh(scs[0].cell)
        lines.append(L(f"  k 点は密度の指定なので、超格子では {'x'.join(map(str, sc_mesh))} になります。", f"  k-points are given as a density, so the supercell uses {'x'.join(map(str, sc_mesh))}."))
    if rec["backend"] == "phonopy":
        lines.append(L(f"  phonopy {rec['phonopy_version']} で作りました (phonopy_disp.yaml)。基本セルの選び方と変位の大きさは、指定のないものは phonopy の既定です。",
                       f"  Made with phonopy {rec['phonopy_version']} (phonopy_disp.yaml); the primitive cell and the displacement are phonopy defaults unless given."))
    else:
        lines.append(L("  phonopy が無かったので、ASE の Phonons と同じ変位 (参照セルの各原子を ±x ±y ±z) を作りました。集めるときも ASE で計算します。",
                       "  phonopy was not available, so the ASE Phonons displacements (each atom of the reference cell by ±x ±y ±z) were made; collection also uses ASE."))
    if rec["dos_mesh"]:
        lines.append(L(f"  状態密度: q 点のメッシュ {'x'.join(map(str, rec['dos_mesh']))}"
                       + (f"、広げる幅 {dos_width_thz:g} THz (ASE の DOS)" if rec["backend"] == "ase" else " (phonopy の既定のテトラヘドロン法)"),
                       f"  DOS: q-point mesh {'x'.join(map(str, rec['dos_mesh']))}" + (f", width {dos_width_thz:g} THz (ASE DOS)" if rec["backend"] == "ase" else " (phonopy default tetrahedron method)")))
    else:
        lines.append(L("  状態密度 (total_dos.dat) は書きません (q 点のメッシュ --phonon-dos-mesh を与えたときだけ書きます。ADIT はメッシュの値を持たないため)。",
                       "  total_dos.dat is not written (only when a q-point mesh is given with --phonon-dos-mesh; ADIT has no default mesh)."))
    lines += ["", L("== 実行したあと ==", "== After running =="),
              L(f"  python {COLLECT_SCRIPT}   (ADIT が入った Python で。力を集めて band.yaml を書きます。phonopy で作ったなら phonopy も要ります)",
                f"  python {COLLECT_SCRIPT}   (with the Python that has ADIT installed; writes band.yaml; phonopy is needed if it was used here)"),
              L("  そのあと adit-analyze <このディレクトリ> でフォノン分散の図 (band.yaml と total_dos.dat を読みます)",
                "  then adit-analyze <this directory> plots the dispersion (it reads band.yaml and total_dos.dat)"),
              L("  振動数の単位は THz。負の値は虚振動数 (phonopy の慣習)。安定かどうかは ADIT は判断しません。", "  Frequencies in THz; negative values are imaginary (phonopy convention); ADIT does not judge stability.")]
    batch.write_top(out, cfg, spec, dirs, lines, "変位の計算", "displacement runs")
    return dirs


def _write_band_yaml_ase(path: Path, paths, labels, conns, freqs_thz: np.ndarray, rec_cell) -> None:
    rec = np.asarray(rec_cell)
    qs = np.vstack(paths)
    lines = ["# written by ADIT from ASE Phonons (not by phonopy). frequency in THz, distance in 1/Angstrom without 2*pi",
             "calculator: ase-phonons", f"nqpoint: {len(qs)}", f"npath: {len(paths)}", "segment_nqpoint:"]
    lines += [f"- {len(p)}" for p in paths]
    lines.append("labels:")
    k = 0
    for i in range(len(paths)):
        lines.append(f"- [ '{labels[k]}', '{labels[k + 1]}' ]")
        k += 1 if conns[i] else 2
    lines.append("phonon:")
    dist, prev, start = 0.0, None, 0
    for p in paths:
        for j, q in enumerate(p):
            qc = q @ rec
            if prev is not None and j > 0:
                dist += float(np.linalg.norm(qc - prev))
            prev = qc
            lines += [f"- q-position: [ {q[0]:.7f}, {q[1]:.7f}, {q[2]:.7f} ]", f"  distance: {dist:.7f}", "  band:"]
            for b, f in enumerate(freqs_thz[start + j], 1):
                lines += [f"  - # {b}", f"    frequency: {f:.10f}"]
        start += len(p)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def collect(out_dir: Path | str) -> dict:
    from adit.outputs import read_forces

    out = Path(out_dir).expanduser()
    rec = json.loads((out / PHONON_FILE).read_text(encoding="utf-8"))
    forces, bad = [], []
    for d in rec["dirs"]:
        try:
            forces.append(read_forces(out / d))
        except Exception as ex:
            bad.append(f"{d}: {ex}")
    if bad:
        raise PhononError(L(f"力を読めない変位が {len(bad)} 個あります (全部の変位が走り終わってから集めます):\n  ", f"{len(bad)} displacements have no readable forces (collect after all of them have finished):\n  ")
                          + "\n  ".join(bad[:10]))
    unit = AtomsData.model_validate(rec["unit_cell"]).to_ase()
    npts = int(rec.get("band_points") or BAND_POINTS)
    mesh = rec.get("dos_mesh")
    files = {}
    if rec["backend"] == "phonopy":
        import phonopy

        ph = phonopy.load(str(out / "phonopy_disp.yaml"), produce_fc=False, log_level=0)
        ph.forces = np.array(forces)
        ph.produce_force_constants()
        ph.save(str(out / "phonopy_params.yaml"))
        paths, labels, conns, text = band_paths(ph.primitive.cell, npts)
        ph.run_band_structure(paths, labels=labels, path_connections=conns)
        ph.write_yaml_band_structure(filename=str(out / "band.yaml"))
        files = {"band": str(out / "band.yaml"), "params": str(out / "phonopy_params.yaml")}
        allf = np.concatenate([np.asarray(x).ravel() for x in ph.band_structure.frequencies])
        if mesh:
            ph.run_mesh(mesh)
            ph.run_total_dos()
            ph.write_total_dos(filename=str(out / "total_dos.dat"))
            files["dos"] = str(out / "total_dos.dat")
    else:
        from ase.phonons import Phonons

        conv = _thz_per_ev()
        with tempfile.TemporaryDirectory() as tmp:
            ph = Phonons(unit, supercell=tuple(rec["dim"]), delta=float(rec["distance_ang"]), name=str(Path(tmp) / "phonon"))
            for key, f in zip(rec["ase_keys"], forces):
                ph.cache[key] = {"forces": f}
            ph.read()
            paths, labels, conns, text = band_paths(unit.cell, npts)
            freqs = np.asarray(ph.band_structure(np.vstack(paths), verbose=False)) * conv
            _write_band_yaml_ase(out / "band.yaml", paths, labels, conns, freqs, unit.cell.reciprocal())
            files = {"band": str(out / "band.yaml")}
            allf = freqs.ravel()
            w = rec.get("dos_width_thz")
            if mesh and w:
                dos = ph.get_dos(kpts=tuple(mesh)).sample_grid(npts=DOS_GRID_POINTS, width=w / conv)
                e, g = dos.get_energies() * conv, dos.get_weights() / conv
                with open(out / "total_dos.dat", "w", encoding="utf-8") as fh:
                    fh.write(f"# written by ADIT from ASE Phonons (q mesh {'x'.join(map(str, mesh))}, Gaussian width {w:g} THz). frequency [THz], DOS [1/THz]\n")
                    fh.writelines(f"{x:20.10f}{y:20.10f}\n" for x, y in zip(e, g))
                files["dos"] = str(out / "total_dos.dat")
    summary = L(f"フォノン ({rec['backend']}): 経路 {text}、振動数 {float(np.min(allf)):.4g}〜{float(np.max(allf)):.4g} THz、負の値 {int((allf < 0).sum())} 個。書いたもの: {', '.join(files.values())}",
                f"phonons ({rec['backend']}): path {text}, frequencies {float(np.min(allf)):.4g} to {float(np.max(allf)):.4g} THz, {int((allf < 0).sum())} negative values. Wrote: {', '.join(files.values())}")
    return {"files": files, "path": text, "min_thz": float(np.min(allf)), "max_thz": float(np.max(allf)), "n_negative": int((allf < 0).sum()), "summary": summary}


if __name__ == "__main__":
    import sys
    print(collect(sys.argv[1])["summary"])


__all__ = ["PhononError", "parse_dim", "band_paths", "plan", "write_phonons", "collect", "has_phonopy"]
