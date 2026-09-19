
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from adit.lang import L
from adit.analysis import plotstyle

_IMG = re.compile(r"^\d{2,}$")
_IMAGE_DIR = re.compile(r"^image_(\d+)$")
_ACT = re.compile(r"activation energy \((->|<-)\)\s*=\s*(-?[\d.]+)\s*eV")
_CP2K_TAIL_BYTES = 1 << 20


def find_neb(run_dir: Path) -> dict | None:
    run_dir = Path(run_dir)
    subs = sorted((p for p in run_dir.iterdir() if p.is_dir() and _IMG.match(p.name)), key=lambda p: int(p.name)) if run_dir.is_dir() else []
    if len(subs) >= 3 and [int(p.name) for p in subs] == list(range(len(subs))):
        return {"kind": "vasp", "dirs": subs}
    imgs = sorted((p for p in run_dir.iterdir() if p.is_dir() and _IMAGE_DIR.match(p.name)),
                  key=lambda p: int(_IMAGE_DIR.match(p.name).group(1))) if run_dir.is_dir() else []
    if len(imgs) >= 3 and [int(_IMAGE_DIR.match(p.name).group(1)) for p in imgs] == list(range(len(imgs))):
        return {"kind": "images", "dirs": imgs}
    if (run_dir / "cp2k.inp").is_file() and _cp2k_band_text(run_dir) is not None:
        return {"kind": "cp2k"}
    stems = []
    for dat in sorted(run_dir.glob("*.dat")) if run_dir.is_dir() else []:
        if dat.with_suffix(".int").is_file():
            stems.append(dat.with_suffix(""))
    if len(stems) == 1:
        return {"kind": "espresso", "prefix": stems[0]}
    if len(stems) > 1:
        return {"kind": "espresso", "prefix": None, "candidates": [str(s) for s in stems]}
    return None


def _numbers(path: Path, ncol: int) -> np.ndarray:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        w = line.split()
        if len(w) >= ncol:
            try:
                rows.append([float(x.replace("D", "E")) for x in w[:ncol]])
            except ValueError:
                continue
    return np.array(rows).reshape(-1, ncol)


def _oszicar_e0(p: Path) -> float | None:
    if not p.is_file():
        return None
    m = re.findall(r"\bE0=\s*(-?[\d.]+(?:[Ee][-+]?\d+)?)", p.read_text(encoding="utf-8", errors="replace"))
    return float(m[-1]) if m else None


def neb_from_images(images) -> dict:
    from ase.mep import NEBTools
    from ase.utils.forcecurve import fit_images

    nt = NEBTools(images)
    fit = fit_images(images)
    b_fit, de = nt.get_barrier(fit=True)
    b_raw, _ = nt.get_barrier(fit=False)
    return {"path_A": [float(x) for x in fit.path], "energies_rel_ev": [float(x) for x in fit.energies],
            "fit_path_A": [float(x) for x in fit.fit_path], "fit_energies_rel_ev": [float(x) for x in fit.fit_energies],
            "forward_barrier_ev": float(b_fit), "backward_barrier_ev": float(b_fit - de), "reaction_energy_ev": float(de),
            "forward_barrier_raw_ev": float(b_raw), "backward_barrier_raw_ev": float(b_raw - de),
            "e_first_ev": float(images[0].get_potential_energy()), "interpolated": True,
            "interpolation": L("ASE の NEBTools (像のエネルギーと、経路に沿った力の成分から 3 次で補間)", "ASE NEBTools (cubic fit through image energies and the force component along the path)")}


def _read_vasp(dirs: list[Path]) -> dict:
    from ase.io import read

    out = {"kind": "vasp", "n_images": len(dirs), "image_dirs": [d.name for d in dirs]}
    images, missing = [], []
    for d in dirs:
        oc = d / "OUTCAR"
        if oc.is_file():
            try:
                images.append(read(str(oc), index=-1, format="vasp-out")); continue
            except Exception:
                pass
        images.append(None); missing.append(d.name)
    if not missing:
        r = neb_from_images(images)
        out.update(r, energy_source=L("各像の OUTCAR の最後のステップ (ASE の読み: energy(sigma→0)、力は TOTAL-FORCE)",
                                      "last step of each image's OUTCAR (ASE reader: energy(sigma->0), forces from TOTAL-FORCE)"))
        return out
    e0 = [_oszicar_e0(d / "OSZICAR") for d in dirs]
    lack = [d.name for d, e in zip(dirs, e0) if e is None]
    if lack:
        out["reasons"] = [L(f"エネルギーの無い像があるので計算していません: {', '.join(lack)} (OUTCAR も OSZICAR も読めません。端の像は、始めと終わりの構造の計算の出力を写すと読めます)",
                            f"not computed; images without an energy: {', '.join(lack)} (neither OUTCAR nor OSZICAR is readable; copy the end-point calculations' outputs into the end-point directories)")]
        return out
    e = np.array(e0) - e0[0]
    pos = []
    for d in dirs:
        a = None
        for name in ("CONTCAR", "POSCAR"):
            if (d / name).is_file() and (d / name).stat().st_size > 0:
                try:
                    a = read(str(d / name), format="vasp"); break
                except Exception:
                    pass
        pos.append(a)
    path = list(range(len(dirs)))
    x_label = "image"
    if all(a is not None for a in pos):
        from ase.geometry import find_mic
        s = [0.0]
        for a, b in zip(pos[:-1], pos[1:]):
            dr, _ = find_mic(b.positions - a.positions, a.cell, a.pbc)
            s.append(s[-1] + float(np.sqrt((dr ** 2).sum())))
        path, x_label = s, "path_A"
    de = float(e[-1])
    out.update(path_A=[float(x) for x in path], x_axis=x_label, energies_rel_ev=[float(x) for x in e], e_first_ev=float(e0[0]),
               forward_barrier_raw_ev=float(e.max()), backward_barrier_raw_ev=float(e.max() - de), reaction_energy_ev=de,
               forward_barrier_ev=None, backward_barrier_ev=None, interpolated=False, interpolation=None,
               energy_source=L("各像の OSZICAR の最後の E0 (sigma→0)", "last E0 (sigma->0) in each image's OSZICAR"),
               reasons=[L(f"OUTCAR の無い像があるので ({', '.join(missing)})、力を使った補間をしていません。障壁は像のエネルギーの最大です",
                          f"some images have no OUTCAR ({', '.join(missing)}), so there is no force-based interpolation; barriers are from the highest image")])
    return out


def _read_qe(prefix: Path) -> dict:
    dat = _numbers(prefix.with_suffix(".dat"), 3)
    itp = _numbers(prefix.with_suffix(".int"), 2)
    out = {"kind": "espresso", "prefix": prefix.name, "files": [prefix.with_suffix(".dat").name, prefix.with_suffix(".int").name]}
    if len(dat) < 2:
        out["reasons"] = [L(f"{prefix.name}.dat に像が 2 つ以上ありません", f"{prefix.name}.dat has fewer than two images")]
        return out
    e = dat[:, 1]
    de = float(e[-1] - e[0])
    out.update(n_images=int(len(dat)), x_axis="reaction_coordinate_normalized", path_A=None, reaction_coordinate=dat[:, 0].tolist(),
               energies_rel_ev=(e - e[0]).tolist(), image_error=dat[:, 2].tolist(),
               image_error_unit=L("neb.x の .dat の 3 列目 (手引きでは eV/a0、標準出力の表の見出しでは eV/A)", "3rd column of the neb.x .dat (eV/a0 in the user guide, eV/A in the stdout table header)"),
               forward_barrier_raw_ev=float(e.max() - e[0]), backward_barrier_raw_ev=float(e.max() - e[-1]), reaction_energy_ev=de,
               energy_source=L(f"{prefix.name}.dat (1 つ目の像からの差 [eV])", f"{prefix.name}.dat (relative to the first image [eV])"))
    if len(itp) >= 2:
        fe = itp[:, 1] - itp[0, 1]
        out.update(fit_reaction_coordinate=itp[:, 0].tolist(), fit_energies_rel_ev=fe.tolist(), interpolated=True,
                   interpolation=L(f"neb.x の {prefix.name}.int (像のエネルギーとその微分からの補間)", f"neb.x {prefix.name}.int (interpolation from image energies and their derivatives)"),
                   forward_barrier_ev=float(fe.max()), backward_barrier_ev=float(fe.max() - de))
    else:
        out.update(interpolated=False, interpolation=None, forward_barrier_ev=None, backward_barrier_ev=None)
    act = {}
    reached_max_steps = False
    for p in sorted(prefix.parent.iterdir()):
        if p.is_file() and p.suffix in (".out", ".log", "") and p.stat().st_size < 50 * 2 ** 20:
            for no, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
                m = _ACT.search(line)
                if m:
                    act[m.group(1)] = {"value_ev": float(m.group(2)), "line": f"{p.name}:{no}"}
                if "reached the maximum number of steps" in line.lower():
                    reached_max_steps = True
    if act:
        out["neb_x_activation"] = {"forward": act.get("->"), "backward": act.get("<-")}
    if reached_max_steps:
        out["stopping_reason"] = "maximum_steps"
        out.setdefault("reasons", []).append(L(
            "neb.x は設定した最大反復数に達して終了しました。収束条件を満たして終了した計算ではありません",
            "neb.x stopped after reaching the configured maximum number of iterations; the run did not finish by meeting its convergence criterion"))
    return out


_CP2K_NUM = r"-?\d+\.\d+(?:[EeDd][-+]?\d+)?"


def _cp2k_band_text(run_dir: Path) -> str | None:
    inp, log = run_dir / "cp2k.inp", run_dir / "output.log"
    if not log.is_file():
        return None
    if not re.search(r"^\s*RUN_TYPE\s+BAND\s*$", inp.read_text(encoding="utf-8", errors="replace"), re.M | re.I):
        return None
    with open(log, "rb") as f:
        f.seek(max(0, log.stat().st_size - _CP2K_TAIL_BYTES))
        text = f.read().decode("utf-8", "replace")
    return text if "ENERGIES [au]" in text else None


def _cp2k_values(text: str, tag: str) -> list[float]:
    lines = text.splitlines()
    ks = [k for k, l in enumerate(lines) if l.strip().startswith(tag)]
    if not ks:
        return []
    k = ks[-1]
    vals = [float(x) for x in re.findall(_CP2K_NUM, lines[k].split("=", 1)[1])] if "=" in lines[k] else []
    for l in lines[k + 1:]:
        if not l.strip() or "=" in l or not re.fullmatch(r"[\s\d.eEdD+-]+", l):
            break
        vals += [float(x) for x in re.findall(_CP2K_NUM, l)]
    return vals


def _read_cp2k(run_dir: Path) -> dict:
    text = _cp2k_band_text(run_dir) or ""
    out = {"kind": "cp2k", "files": ["output.log"]}
    e_h = _cp2k_values(text, "ENERGIES [au]")
    if len(e_h) < 2:
        out["reasons"] = [L("output.log の最後の「ENERGIES [au]」に像が 2 つ以上ありません", "the last 'ENERGIES [au]' block of output.log has fewer than two images")]
        return out
    e = (np.array(e_h) - e_h[0]) * 27.211386245988
    de = float(e[-1])
    m = re.findall(r"NUMBER OF NEB REPLICA\s*=\s*(\d+)", text)
    step = re.findall(r"STEP NUMBER\s*=\s*(\d+)", text)
    kind = re.findall(r"BAND TYPE\s*=\s*(\S+)", text)
    out.update(n_images=len(e_h), x_axis="image", path_A=list(range(len(e_h))), energies_rel_ev=[float(x) for x in e],
               e_first_ev=float(e_h[0] * 27.211386245988), energies_hartree=[float(x) for x in e_h],
               forward_barrier_raw_ev=float(e.max()), backward_barrier_raw_ev=float(e.max() - de), reaction_energy_ev=de,
               forward_barrier_ev=None, backward_barrier_ev=None, interpolated=False, interpolation=None,
               replica_distances=_cp2k_values(text, "DISTANCES REP"),
               replica_distances_note=L("像の間の距離 (CP2K の DISTANCES REP)。単位は出力に書かれていないので、そのまま並べています",
                                        "distances between replicas (CP2K DISTANCES REP); the unit is not written in the output, so they are listed as they are"),
               n_replica_in_output=int(m[-1]) if m else None, step_number=int(step[-1]) if step else None,
               band_type=kind[-1] if kind else None,
               energy_source=L("output.log の最後の「ENERGIES [au]」 (Hartree を eV に)", "the last 'ENERGIES [au]' block of output.log (hartree to eV)"))
    return out


def _read_images(dirs: list[Path]) -> dict:
    import json

    from adit.analysis.readers import load_run

    out = {"kind": "images", "n_images": len(dirs), "image_dirs": [d.name for d in dirs]}
    energies, missing, atoms = [], [], []
    for d in dirs:
        try:
            r = load_run(d)
        except Exception as ex:
            energies.append(None); atoms.append(None); missing.append(f"{d.name}: {ex}")
            continue
        energies.append(float(r.energies_ev[-1]) if r.energies_ev else None)
        atoms.append(r.final)
        if not r.energies_ev:
            missing.append(L(f"{d.name}: エネルギーがありません (まだ計算していないかもしれません)", f"{d.name}: no energy (perhaps not run yet)"))
    if missing:
        out["reasons"] = [L(f"エネルギーの無い像があるので計算していません: {'; '.join(missing[:5])}",
                            f"not computed; images without an energy: {'; '.join(missing[:5])}")]
        return out
    rec = dirs[0].parent / "neb.json"
    path = None
    if rec.is_file():
        try:
            js = json.loads(rec.read_text(encoding="utf-8"))
            p = js.get("path_length_ang")
            if isinstance(p, list) and len(p) == len(dirs):
                path = [float(x) for x in p]
        except ValueError:
            pass
    e = np.array(energies, dtype=float) - energies[0]
    de = float(e[-1])
    out.update(path_A=path if path is not None else list(range(len(dirs))), x_axis="path_A" if path is not None else "image",
               energies_rel_ev=[float(x) for x in e], e_first_ev=float(energies[0]),
               forward_barrier_raw_ev=float(e.max()), backward_barrier_raw_ev=float(e.max() - de), reaction_energy_ev=de,
               forward_barrier_ev=None, backward_barrier_ev=None, interpolated=False, interpolation=None,
               energy_source=L("各像のディレクトリの最終エネルギー (解析の読み取り)", "final energy of each image directory (analysis readers)"),
               reasons=[L("像は補間した位置のままの一点計算です (NEB の最適化はしていません)",
                          "the images are single points on the interpolated positions (no NEB optimization)")])
    if path is None:
        out["reasons"].append(L("neb.json の path_length_ang が読めないので、横軸は像の番号です", "path_length_ang in neb.json cannot be read, so the x axis is the image number"))
    if all(a is not None for a in atoms):
        try:
            from ase.calculators.singlepoint import SinglePointCalculator

            from adit.outputs import read_forces
            images = []
            for d, a, en in zip(dirs, atoms, energies):
                img = a.copy()
                img.calc = SinglePointCalculator(img, energy=en, forces=read_forces(d))
                images.append(img)
            fit = neb_from_images(images)
            reasons = out["reasons"]
            out.update(fit)
            out["reasons"] = reasons
        except Exception as ex:
            out["reasons"].append(L(f"力を読めないので補間していません: {ex}", f"no interpolation; the forces cannot be read: {ex}"))
    return out


def analyze_neb(run_dir: Path, out_dir: Path) -> dict | None:
    found = find_neb(Path(run_dir))
    if found is None:
        return None
    if found["kind"] == "vasp":
        t = _read_vasp(found["dirs"])
    elif found["kind"] == "images":
        t = _read_images(found["dirs"])
    elif found["kind"] == "cp2k":
        t = _read_cp2k(Path(run_dir))
    elif found.get("prefix") is None:
        return {"kind": "espresso", "reasons": [L(f"neb.x の .dat と .int の組が複数あり、どれを読むか決められません: {found['candidates']}",
                                                  f"several neb.x .dat/.int pairs; cannot tell which to read: {found['candidates']}")]}
    else:
        t = _read_qe(found["prefix"])
    rec = Path(run_dir) / "neb.json"
    if rec.is_file():
        try:
            climb = bool(json.loads(rec.read_text(encoding="utf-8")).get("options", {}).get("climb"))
            t["climbing_image"] = climb
            if not climb:
                t.setdefault("reasons", []).append(L(
                    "climbing image は使っていません。最高エネルギーの像は鞍点そのものとは限らないため、像の最大値から得る上昇量は、遷移状態を探して求めた障壁ではありません",
                    "no climbing image was used; the highest-energy image need not coincide with the saddle point, so the rise to the highest image is not a barrier from a located transition state"))
        except (OSError, ValueError, TypeError):
            pass
    if "energies_rel_ev" in t:
        _plot(t, Path(out_dir) / "neb.png")
        t["figure"] = str(Path(out_dir) / "neb.png")
    return t


def _plot(t: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 3.4))
    if t["kind"] == "espresso":
        x, xf, xl = t["reaction_coordinate"], t.get("fit_reaction_coordinate"), "reaction coordinate (normalized)"
    else:
        x, xf = t["path_A"], t.get("fit_path_A")
        xl = "path [Å]" if t.get("x_axis", "path_A") == "path_A" else "image"
    ax.plot(x, t["energies_rel_ev"], "o", label="images")
    if t.get("interpolated") and xf is not None:
        ax.plot(xf, t["fit_energies_rel_ev"], "-", color="k", lw=1.2, label="interpolation")
    fb = t.get("forward_barrier_ev") if t.get("forward_barrier_ev") is not None else t["forward_barrier_raw_ev"]
    kind = "interpolated" if t.get("forward_barrier_ev") is not None else "highest image"
    ax.set_title(f"forward {fb:.3f} eV, reaction {t['reaction_energy_ev']:+.3f} eV ({kind})", fontsize=9)
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel(xl); ax.set_ylabel("E - E(first image) [eV]"); plotstyle.grid(ax); ax.legend(fontsize=8)
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)


def summary_lines(t: dict) -> list[str]:
    if "energies_rel_ev" not in t:
        return [L("NEB: ", "NEB: ") + " ".join(t.get("reasons", []))]
    fb, bb = t.get("forward_barrier_ev"), t.get("backward_barrier_ev")
    s = L(f"NEB ({t['kind']}、像 {t['n_images']} 個): 反応エネルギー {t['reaction_energy_ev']:+.4f} eV、"
          f"像の最大からの障壁 前向き {t['forward_barrier_raw_ev']:.4f} eV / 後ろ向き {t['backward_barrier_raw_ev']:.4f} eV",
          f"NEB ({t['kind']}, {t['n_images']} images): reaction energy {t['reaction_energy_ev']:+.4f} eV, "
          f"barrier from the highest image forward {t['forward_barrier_raw_ev']:.4f} eV / backward {t['backward_barrier_raw_ev']:.4f} eV")
    if fb is not None:
        s += L(f"、補間した曲線からの障壁 前向き {fb:.4f} eV / 後ろ向き {bb:.4f} eV ({t['interpolation']})",
               f"; from the interpolated curve forward {fb:.4f} eV / backward {bb:.4f} eV ({t['interpolation']})")
    else:
        s += L("、補間なし", "; no interpolation")
    lines = [s]
    a = t.get("neb_x_activation")
    if a:
        f = f"{a['forward']['value_ev']:.6f} eV ({a['forward']['line']})" if a.get("forward") else "-"
        b = f"{a['backward']['value_ev']:.6f} eV ({a['backward']['line']})" if a.get("backward") else "-"
        lines.append(L(f"  neb.x の標準出力の activation energy: (->) {f}、(<-) {b}", f"  neb.x stdout activation energy: (->) {f}, (<-) {b}"))
    lines += ["  " + r for r in t.get("reasons", [])]
    return lines


__all__ = ["find_neb", "analyze_neb", "neb_from_images", "summary_lines"]
