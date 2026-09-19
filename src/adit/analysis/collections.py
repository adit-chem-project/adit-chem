
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from adit.lang import L
from adit.analysis import plotstyle

PHONON_FILE = "phonons.json"
ELASTIC_FILE = "elastic.json"
CONFORMER_FILE = "conformers.json"
NEB_FILE = "neb.json"
COMPARE_FILE = "compare.json"
EV_KCAL_MOL = 23.060547830619

_FILES = ((PHONON_FILE, "phonons"), (ELASTIC_FILE, "elastic"), (CONFORMER_FILE, "conformers"),
          (NEB_FILE, "neb_images"), (COMPARE_FILE, "compare"))


def detect_collection(run_dir: Path | str) -> str | None:
    run_dir = Path(run_dir)
    for name, kind in _FILES:
        if (run_dir / name).is_file():
            return kind
    if (run_dir / "scan.json").is_file():
        return "scan"
    if (run_dir / "crest_conformers.xyz").is_file():
        return "crest"
    return None


def _load(path: Path) -> dict:
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as ex:
        raise ValueError(L(f"{path.name} を読めません: {ex}", f"cannot read {path.name}: {ex}")) from ex
    if not isinstance(obj, dict):
        raise ValueError(L(f"{path.name} の中身が辞書ではありません", f"{path.name} does not contain an object"))
    return obj


def analyze_phonon_set(run_dir: Path, collect: bool = False) -> dict:
    run_dir = Path(run_dir)
    t = {"record": PHONON_FILE, "reasons": []}
    try:
        rec = _load(run_dir / PHONON_FILE)
    except ValueError as ex:
        t["reasons"].append(str(ex))
        return t
    for k in ("backend", "dim", "code", "n_displacements", "distance_ang", "dos_mesh", "dos_width_thz", "phonopy_version"):
        if k in rec:
            t[k] = rec[k]
    t["n_dirs"] = len(rec.get("dirs") or [])
    have = {k: (run_dir / n).is_file() or (run_dir / "phonopy" / n).is_file()
            for k, n in (("band", "band.yaml"), ("dos", "total_dos.dat"))}
    t["collected"] = bool(have["band"])
    t["has_dos_file"] = bool(have["dos"])
    if not have["band"] and collect:
        try:
            from adit.phonon_setup import collect as phonon_collect
            res = phonon_collect(run_dir)
            t["collect_summary"] = res["summary"]
            t["collected"] = bool((run_dir / "band.yaml").is_file())
            t["has_dos_file"] = (run_dir / "total_dos.dat").is_file()
        except Exception as ex:
            t["reasons"].append(L(f"集計できません ({ex})。各変位の計算が終わっているか確かめてください",
                                  f"cannot collect ({ex}); check that every displacement has finished"))
    elif not have["band"]:
        t["reasons"].append(L("フォノン分散 (band.yaml) がまだありません。各変位の計算を実行したあと、このディレクトリで phonon_collect.py "
                              "を実行してください (adit-analyze に --collect を付けると、ここから同じ集計を呼びます)",
                              "there is no phonon dispersion (band.yaml) yet; after running every displacement, run phonon_collect.py in this "
                              "directory (or add --collect to adit-analyze to call the same collection from here)"))
    if t.get("dos_mesh") is None and not t["has_dos_file"]:
        t["reasons"].append(L("状態密度 (total_dos.dat) は、生成のときに q 点のメッシュを指定した場合だけ作られます",
                              "total_dos.dat is produced only if a q-point mesh was given when the runs were generated"))
    return t


def _num(v, fmt: str = ".4g") -> str:
    try:
        return format(float(v), fmt)
    except (TypeError, ValueError):
        return "?"


def phonon_set_lines(t: dict) -> list[str]:
    s = [L(f"フォノンの計算の組: {t.get('backend', '?')}、超格子 {'x'.join(map(str, t.get('dim') or []))}、変位 {t.get('n_displacements', '?')} 個"
           f" ({_num(t.get('distance_ang'), 'g')} Å)、集計 {'済み' if t.get('collected') else 'まだ'}",
           f"phonon run set: {t.get('backend', '?')}, supercell {'x'.join(map(str, t.get('dim') or []))}, "
           f"{t.get('n_displacements', '?')} displacements of {_num(t.get('distance_ang'), 'g')} Å, "
           f"{'collected' if t.get('collected') else 'not collected yet'}")]
    if t.get("collect_summary"):
        s.append("  " + str(t["collect_summary"]))
    return s + ["  " + r for r in t.get("reasons", [])]


def analyze_elastic(run_dir: Path, out_dir: Path, collect: bool = False) -> dict:
    run_dir, out_dir = Path(run_dir), Path(out_dir)
    t = {"record": ELASTIC_FILE, "reasons": []}
    try:
        rec = _load(run_dir / ELASTIC_FILE)
    except ValueError as ex:
        t["reasons"].append(str(ex))
        return t
    if "c_gpa" not in rec and collect:
        try:
            from adit.elastic_setup import collect as elastic_collect
            elastic_collect(run_dir)
            rec = _load(run_dir / ELASTIC_FILE)
        except Exception as ex:
            t["reasons"].append(L(f"集計できません ({ex})。各歪みの計算が終わっているか確かめてください",
                                  f"cannot collect ({ex}); check that every strain run has finished"))
    for k in ("code", "task", "components", "strains", "convention"):
        if k in rec:
            t[k] = rec[k]
    t["n_runs"] = len(rec.get("runs") or [])
    if "c_gpa" not in rec:
        t["collected"] = False
        t["reasons"].append(L("弾性定数 (C_ij) がまだありません。各歪みの計算を実行したあと、このディレクトリで elastic_collect.py を"
                              "実行してください (adit-analyze に --collect を付けると、ここから同じ集計を呼びます)",
                              "the elastic constants (C_ij) are not there yet; after running every strain, run elastic_collect.py in this "
                              "directory (or add --collect to adit-analyze to call the same collection from here)"))
        return t
    t["collected"] = True
    t["c_gpa"] = rec["c_gpa"]
    t["fits"] = rec.get("fits") or []
    t["stress_voigt_gpa"] = rec.get("stress_voigt_gpa") or {}
    t["max_abs_residual_gpa"] = max((f["max_abs_residual_gpa"] for f in t["fits"]), default=None)
    t["unit"] = "GPa"
    e0 = t["stress_voigt_gpa"].get("e0")
    if e0 is not None and len(e0) == 6:
        t["unstrained_stress_voigt_gpa"] = [float(x) for x in e0]
        t["unstrained_pressure_gpa"] = float(-np.mean([float(x) for x in e0[:3]]))
        t["unstrained_pressure_note"] = L("歪みなし (e0) の計算の圧力 = -(σxx + σyy + σzz) / 3 (ASE の約束: 圧縮で負)。"
                                          "直線の当てはめの切片 (intercept_gpa) は、当てはめから見た δ = 0 の応力です",
                                          "pressure of the unstrained run (e0) = -(sxx + syy + szz) / 3 (ASE convention: compression negative); "
                                          "the fit intercepts (intercept_gpa) are the stresses at δ = 0 as seen by the fits")
    from adit.analysis.elastic_moduli import ElasticModuliError, from_cij, summary_lines

    try:
        moduli = from_cij(t["c_gpa"])
    except ElasticModuliError as ex:
        t["moduli_error"] = str(ex)
    else:
        t["moduli"] = moduli.as_dict()
        t["moduli_summary"] = summary_lines(moduli)
    t["files"] = {"fits": str(_write_elastic_csv(out_dir, t["fits"]))}
    fig = _plot_elastic(rec, out_dir)
    if fig:
        t["figure_stress_strain"] = fig
    return t


def _write_elastic_csv(out_dir: Path, fits: list[dict]) -> Path:
    p = Path(out_dir) / "elastic_fits.csv"
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["i", "j", "c_gpa", "intercept_gpa", "max_abs_residual_gpa", "n_points"])
        for x in fits:
            w.writerow([x["i"], x["j"], f"{x['c_gpa']:.6f}", f"{x['intercept_gpa']:.6f}", f"{x['max_abs_residual_gpa']:.6g}", x["n_points"]])
    return p


def _plot_elastic(rec: dict, out_dir: Path) -> str | None:
    runs = rec.get("runs") or []
    stress = rec.get("stress_voigt_gpa") or {}
    comps = [j for j in (rec.get("components") or []) if any(r["component"] == j for r in runs)]
    if not comps or not stress:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fits = {(f["i"], f["j"]): f for f in (rec.get("fits") or [])}
    fig, axes = plt.subplots(1, len(comps), figsize=(3.2 * len(comps), 3.2), squeeze=False)
    for ax, j in zip(axes[0], comps):
        pts = [(r["strain"], stress[r["dir"]]) for r in runs if r["component"] in (0, j) and r["dir"] in stress]
        pts.sort(key=lambda p: p[0])
        x = np.array([p[0] for p in pts])
        for i in range(1, 7):
            y = np.array([p[1][i - 1] for p in pts])
            line, = ax.plot(x, y, "o", ms=4, label=f"s{i}")
            f = fits.get((i, j))
            if f is not None and len(x) >= 2:
                xx = np.linspace(x.min(), x.max(), 50)
                ax.plot(xx, f["c_gpa"] * xx + f["intercept_gpa"], "-", lw=1.0, color=line.get_color())
        ax.set_xlabel(f"strain {j} (Voigt)"); plotstyle.grid(ax)
        ax.axhline(0, color="gray", lw=0.6); ax.axvline(0, color="gray", lw=0.6)
    axes[0][0].set_ylabel("stress [GPa] (ASE sign)")
    axes[0][-1].legend(fontsize=7, ncol=2)
    p = Path(out_dir) / "elastic_stress_strain.png"
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return str(p)


def elastic_lines(t: dict) -> list[str]:
    if not t.get("collected"):
        return [L("弾性定数の計算の組: まだ集計していません", "elastic run set: not collected yet")] + ["  " + r for r in t.get("reasons", [])]
    c = t["c_gpa"]
    s = [L(f"弾性定数 C_ij [GPa] (行 i = 応力の成分、列 j = 歪みの成分。Voigt 1=xx 2=yy 3=zz 4=yz 5=xz 6=xy。対称化していません。"
           f"当てはめの残差の最大 {_num(t.get('max_abs_residual_gpa'))} GPa):",
           f"elastic constants C_ij [GPa] (row i = stress, column j = strain; Voigt 1=xx ... 6=xy; not symmetrized; "
           f"largest fit residual {_num(t.get('max_abs_residual_gpa'))} GPa):")]
    s += ["  " + " ".join(f"{('-' if v is None else f'{v:9.2f}'):>9}" for v in row) for row in c]
    if t.get("unstrained_pressure_gpa") is not None:
        s.append(L(f"  歪みなしの計算の圧力: {t['unstrained_pressure_gpa']:+.4f} GPa", f"  pressure of the unstrained run: {t['unstrained_pressure_gpa']:+.4f} GPa"))
    return s + ["  " + r for r in t.get("reasons", [])]


def analyze_conformers(run_dir: Path, out_dir: Path) -> dict:
    run_dir, out_dir = Path(run_dir), Path(out_dir)
    t = {"record": CONFORMER_FILE, "reasons": []}
    try:
        rec = _load(run_dir / CONFORMER_FILE)
    except ValueError as ex:
        t["reasons"].append(str(ex))
        return t
    from adit.analysis.readers import load_run

    for k in ("source", "force_field", "requested", "embedded", "kept", "seed", "rmsd_threshold_ang", "rmsd_atoms"):
        if k in rec:
            t[k] = rec[k]
    rows = []
    for c in rec.get("conformers") or []:
        if not c.get("dir"):
            continue
        row = {"rank": c.get("rank"), "dir": c["dir"], "rdkit_id": c.get("rdkit_id"), "ff_energy_kcal_mol": c.get("energy_kcal_mol"),
               "ff_relative_kcal_mol": c.get("relative_kcal_mol"), "ff_converged": c.get("converged"),
               "code": "", "energy_ev": None, "note": ""}
        d = run_dir / c["dir"]
        try:
            r = load_run(d)
            row["code"] = r.code
            if r.energies_ev:
                row["energy_ev"] = float(r.energies_ev[-1])
            else:
                row["note"] = L("エネルギーがありません (まだ計算していないかもしれません)", "no energy (perhaps not run yet)")
        except Exception as ex:
            row["note"] = L(f"読めません: {ex}", f"cannot read: {ex}")
        rows.append(row)
    got = [r for r in rows if r["energy_ev"] is not None]
    if got:
        lo = min(r["energy_ev"] for r in got)
        for r in rows:
            if r["energy_ev"] is not None:
                r["relative_ev"] = r["energy_ev"] - lo
                r["relative_kcal_mol"] = (r["energy_ev"] - lo) * EV_KCAL_MOL
    ff = [r["ff_energy_kcal_mol"] for r in rows if r.get("ff_energy_kcal_mol") is not None]
    if ff:
        lo_ff = min(ff)
        for r in rows:
            if r.get("ff_energy_kcal_mol") is not None:
                r["ff_relative_kcal_mol"] = r["ff_energy_kcal_mol"] - lo_ff
    t["conformers"] = rows
    t["n_with_energy"] = len(got)
    t["relative_note"] = L("相対エネルギーは、それぞれの欄 (力場・計算コード) の中で最も低い配座からの差です。どの配座を使うかは判定しません",
                           "relative energies are measured from the lowest conformer within each column (force field, code); no conformer is chosen for you")
    t["files"] = {"table": str(_write_conformer_csv(out_dir, rows))}
    fig = _plot_conformers(rows, out_dir, str(t.get("force_field") or ""))
    if fig:
        t["figure"] = fig
    return t


_CONF_COLS = ["rank", "dir", "rdkit_id", "ff_energy_kcal_mol", "ff_relative_kcal_mol", "ff_converged", "code", "energy_ev",
              "relative_ev", "relative_kcal_mol", "note"]


def _write_conformer_csv(out_dir: Path, rows: list[dict]) -> Path:
    p = Path(out_dir) / "conformers_energies.csv"
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_CONF_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in _CONF_COLS})
    return p


def _plot_conformers(rows: list[dict], out_dir: Path, force_field: str) -> str | None:
    got = [r for r in rows if r.get("relative_kcal_mol") is not None or r.get("ff_relative_kcal_mol") is not None]
    if len(got) < 2:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(got))
    fig, ax = plt.subplots(figsize=(0.7 * len(got) + 3.0, 3.2))
    ax.bar(x - 0.2, [r.get("ff_relative_kcal_mol") or 0.0 for r in got], width=0.4, label=force_field or "force field")
    ax.bar(x + 0.2, [r.get("relative_kcal_mol") or 0.0 for r in got], width=0.4, label=(got[0].get("code") or "code"))
    ax.set_xticks(x); ax.set_xticklabels([r["dir"] for r in got], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("E - E(lowest) [kcal/mol]"); ax.grid(alpha=0.3, axis="y"); ax.legend(fontsize=8)
    p = Path(out_dir) / "conformers.png"
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return str(p)


def conformer_lines(t: dict) -> list[str]:
    rows = t.get("conformers") or []
    s = [L(f"配座の候補: {len(rows)} 個の計算 (力場 {t.get('force_field', '?')}、元 {t.get('source', '?')}、"
           f"エネルギーを読めたもの {t.get('n_with_energy', 0)} 個)",
           f"conformer candidates: {len(rows)} runs (force field {t.get('force_field', '?')}, from {t.get('source', '?')}, "
           f"{t.get('n_with_energy', 0)} with an energy)")]
    s.append(L("  ディレクトリ  力場 [kcal/mol]  計算コード [eV]  最も低い配座との差 [kcal/mol]",
               "  directory    force field [kcal/mol]  code [eV]  above the lowest [kcal/mol]"))
    for r in rows:
        e = f"{r['energy_ev']:.6f}" if r.get("energy_ev") is not None else "-"
        rel = f"{r['relative_kcal_mol']:.3f}" if r.get("relative_kcal_mol") is not None else "-"
        ff = f"{r['ff_relative_kcal_mol']:.3f}" if r.get("ff_relative_kcal_mol") is not None else "-"
        s.append(f"  {r['dir']}  {ff}  {e}  {rel}" + (f"  ({r['note']})" if r.get("note") else ""))
    return s + ["  " + r for r in t.get("reasons", [])]


__all__ = ["detect_collection", "analyze_phonon_set", "analyze_elastic", "analyze_conformers",
           "phonon_set_lines", "elastic_lines", "conformer_lines"]
