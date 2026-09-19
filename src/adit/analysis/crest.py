"""Read a finished CREST conformer search: relative energies, degeneracies, heavy-atom RMSD and, when a temperature is given, Boltzmann weights."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
from ase import Atoms

from adit.analysis import plotstyle
from adit.analysis.free_energy import KB_EV_K
from adit.analysis.readers import HARTREE_EV
from adit.errors import AditValueError
from adit.lang import L

ENSEMBLE_FILE = "crest_conformers.xyz"
ENERGIES_FILE = "crest.energies"
LOG_FILE = "crest.log"
EV_KCAL_MOL = 23.060547830619
EV_KJ_MOL = 96.485332123
# CREST docs (example_1): "Erel/kcal  Etot  weight/tot  conformer  set  degen  origin"
_TABLE_HEAD = re.compile(r"Erel/kcal\s+Etot\s+weight/tot")


class CrestError(AditValueError):
    pass


def find_crest_dir(run_dir: Path | str) -> Path | None:
    run_dir = Path(run_dir)
    for d in (run_dir, run_dir / "crest"):
        if (d / ENSEMBLE_FILE).is_file():
            return d
    return None


def read_ensemble(path: Path | str) -> list[tuple[float | None, Atoms]]:
    """Plain xyz ensemble; the comment line carries the total energy in Hartree (CREST cregen.f90, format f18.8)."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    out: list[tuple[float | None, Atoms]] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].strip())
        except ValueError as ex:
            raise CrestError(L(f"{Path(path).name} の {i + 1} 行目が原子数ではありません: {lines[i].strip()[:40]!r}",
                               f"line {i + 1} of {Path(path).name} is not an atom count: {lines[i].strip()[:40]!r}")) from ex
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        energy = None
        for tok in comment.split():
            try:
                energy = float(tok)
                break
            except ValueError:
                continue
        block = lines[i + 2:i + 2 + n]
        if len(block) < n:
            raise CrestError(L(f"{Path(path).name} の最後の構造が途中で終わっています", f"the last structure in {Path(path).name} is truncated"))
        syms, pos = [], []
        for raw in block:
            w = raw.split()
            try:
                syms.append(w[0]); pos.append([float(w[1]), float(w[2]), float(w[3])])
            except (IndexError, ValueError) as ex:
                raise CrestError(L(f"{Path(path).name} の座標の行を読めません: {raw.strip()[:60]!r}",
                                   f"cannot read a coordinate line of {Path(path).name}: {raw.strip()[:60]!r}")) from ex
        out.append((energy, Atoms(symbols=syms, positions=pos)))
        i += 2 + n
    if not out:
        raise CrestError(L(f"{Path(path).name} に構造がありません", f"{Path(path).name} contains no structure"))
    return out


def read_energies(path: Path | str) -> dict[int, float]:
    """crest.energies: '<index> <E - E_lowest in kcal/mol>' per conformer (CREST cregen.f90)."""
    out: dict[int, float] = {}
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        w = raw.split()
        if len(w) < 2:
            continue
        try:
            out[int(w[0])] = float(w[1])
        except ValueError:
            continue
    return out


def read_log_degeneracy(path: Path | str) -> dict[int, int]:
    """Degeneracy per conformer from the final ensemble table of the CREST log (column 'degen'); empty when absent."""
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return {}
    heads = [k for k, l in enumerate(lines) if _TABLE_HEAD.search(l)]
    if not heads:
        return {}
    out: dict[int, int] = {}
    for raw in lines[heads[-1] + 1:]:
        w = raw.split()
        if not w:
            break
        try:
            nums = [float(x) for x in w]
        except ValueError:
            break
        # rows that start a conformer group carry: idx Erel Etot weight group_weight set degen [origin]
        if len(nums) >= 7:
            out[int(nums[5])] = int(nums[6])
    return out


def boltzmann_weights(relative_ev, degeneracy, temperature_k: float) -> np.ndarray:
    """w_i = g_i exp(-dE_i / k_B T) / sum_j g_j exp(-dE_j / k_B T), k_B = 8.617333262e-5 eV/K."""
    if temperature_k <= 0:
        raise CrestError(L("温度は正の値 [K] にしてください", "the temperature must be positive [K]"))
    de = np.asarray(relative_ev, dtype=float)
    g = np.asarray(degeneracy, dtype=float)
    x = g * np.exp(-(de - de.min()) / (KB_EV_K * temperature_k))
    return x / x.sum()


def _heavy_rmsd(ref: Atoms, other: Atoms) -> float | None:
    from adit.analysis.geometry_series import kabsch_rmsd

    a = np.array(ref.get_chemical_symbols()) != "H"
    b = np.array(other.get_chemical_symbols()) != "H"
    if a.sum() != b.sum() or not a.any():
        return None
    return kabsch_rmsd(ref.positions[a], other.positions[b])[0]


def analyze_crest(crest_dir: Path | str, out_dir: Path | str, temperature_k: float | None = None) -> dict:
    """Table of the CREST conformers in crest_dir; writes crest_conformers.csv and crest_conformers.png to out_dir."""
    crest_dir, out_dir = Path(crest_dir), Path(out_dir)
    t: dict = {"source": str(crest_dir / ENSEMBLE_FILE), "reasons": [], "temperature_k": temperature_k}
    ens = read_ensemble(crest_dir / ENSEMBLE_FILE)
    file_rel = read_energies(crest_dir / ENERGIES_FILE) if (crest_dir / ENERGIES_FILE).is_file() else {}
    degen = read_log_degeneracy(crest_dir / LOG_FILE)
    if not (crest_dir / ENERGIES_FILE).is_file():
        t["reasons"].append(L(f"{ENERGIES_FILE} が無いので、相対エネルギーは {ENSEMBLE_FILE} のコメント行 (全エネルギー [Eh]) から出しました",
                              f"no {ENERGIES_FILE}; relative energies come from the comment lines of {ENSEMBLE_FILE} (total energy in Eh)"))
    if not degen:
        t["reasons"].append(L(f"縮退度は {LOG_FILE} の最後の表 (列 degen) から読みます。無いので全部 1 としました",
                              f"degeneracies are read from the final table of {LOG_FILE} (column degen); none found, so every conformer counts as 1"))
    abs_eh = [e for e, _ in ens]
    have_abs = all(e is not None for e in abs_eh)
    rows = []
    for k, (e_eh, atoms) in enumerate(ens, start=1):
        row = {"index": k, "n_atoms": len(atoms), "energy_eh": e_eh, "energy_ev": None if e_eh is None else e_eh * HARTREE_EV,
               "crest_energies_kcal_mol": file_rel.get(k), "degeneracy": int(degen.get(k, 1)),
               "rmsd_heavy_A": _heavy_rmsd(ens[0][1], atoms)}
        rows.append(row)
    if have_abs:
        lo = min(abs_eh)
        for r in rows:
            r["relative_ev"] = (r["energy_eh"] - lo) * HARTREE_EV
        t["relative_source"] = L(f"{ENSEMBLE_FILE} のコメント行の全エネルギー", f"total energies on the comment lines of {ENSEMBLE_FILE}")
    elif file_rel:
        for r in rows:
            r["relative_ev"] = None if r["crest_energies_kcal_mol"] is None else r["crest_energies_kcal_mol"] / EV_KCAL_MOL
        t["relative_source"] = ENERGIES_FILE
    else:
        raise CrestError(L(f"{ENSEMBLE_FILE} のコメント行にエネルギーが無く、{ENERGIES_FILE} もありません",
                           f"no energy on the comment lines of {ENSEMBLE_FILE} and no {ENERGIES_FILE}"))
    for r in rows:
        rel = r.get("relative_ev")
        r["relative_kcal_mol"] = None if rel is None else rel * EV_KCAL_MOL
        r["relative_kj_mol"] = None if rel is None else rel * EV_KJ_MOL
    if any(r["rmsd_heavy_A"] is None for r in rows):
        t["reasons"].append(L("重原子の数が配座の間で違うか重原子が無いので、RMSD を出せない配座があります",
                              "the heavy-atom count differs between conformers (or there is none), so some RMSDs are missing"))
    if temperature_k is not None:
        rel = [r["relative_ev"] for r in rows]
        if any(v is None for v in rel):
            t["reasons"].append(L("相対エネルギーの無い配座があるので、重みは出していません", "a conformer has no relative energy, so no weights are given"))
        else:
            w = boltzmann_weights(rel, [r["degeneracy"] for r in rows], temperature_k)
            for r, wi in zip(rows, w):
                r["weight"] = float(wi)
            t["population_lowest"] = float(w[0])
            t["weight_formula"] = "w_i = g_i exp(-dE_i / k_B T) / sum_j g_j exp(-dE_j / k_B T); k_B = 8.617333262e-5 eV/K"
    else:
        t["reasons"].append(L("温度が指定されていないので、Boltzmann の重みは出していません (--conformer-temperature K で指定します)",
                              "no temperature given, so Boltzmann weights are not computed (give --conformer-temperature K)"))
    t["conformers"] = rows
    t["n_conformers"] = len(rows)
    t["rmsd_note"] = L("RMSD は最も低い配座 (1 番) との重原子の Kabsch 重ね合わせ。原子の並べ替え (対称な置換) は考えていません",
                       "RMSD: heavy-atom Kabsch superposition on the lowest conformer (index 1); symmetric atom permutations are not considered")
    t["files"] = {"table": str(_write_csv(out_dir, rows))}
    fig = _plot(rows, out_dir, temperature_k)
    if fig:
        t["figure"] = fig
    return t


_COLS = ["index", "n_atoms", "energy_eh", "relative_kcal_mol", "relative_kj_mol", "relative_ev", "crest_energies_kcal_mol",
         "degeneracy", "rmsd_heavy_A", "weight"]


def _write_csv(out_dir: Path, rows: list[dict]) -> Path:
    p = out_dir / "crest_conformers.csv"
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in _COLS})
    return p


def _plot(rows: list[dict], out_dir: Path, temperature_k: float | None) -> str | None:
    got = [r for r in rows if r.get("relative_kcal_mol") is not None]
    if not got:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(got))
    fig, ax = plt.subplots(figsize=(min(0.35 * len(got) + 3.0, 12.0), 3.2))
    ax.bar(x, [r["relative_kcal_mol"] for r in got], width=0.7, label="E - E(lowest)")
    ax.set_xticks(x); ax.set_xticklabels([str(r["index"]) for r in got], fontsize=7)
    ax.set_xlabel("conformer (CREST index)"); ax.set_ylabel("E - E(lowest) [kcal/mol]"); plotstyle.grid(ax)
    if temperature_k is not None and all("weight" in r for r in got):
        ax2 = ax.twinx()
        ax2.plot(x, [r["weight"] for r in got], "o-", ms=3, lw=0.8, color="tab:red", label=f"weight ({temperature_k:g} K)")
        ax2.set_ylabel("Boltzmann weight"); ax2.set_ylim(0, 1)
        ax2.legend(loc="upper right", fontsize=8)
    p = out_dir / "crest_conformers.png"
    fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    return str(p)


def summary_lines(t: dict) -> list[str]:
    rows = t.get("conformers") or []
    temp = t.get("temperature_k")
    s = [L(f"CREST の配座: {len(rows)} 個 ({t.get('source', '?')})" + (f"、重みの温度 {temp:g} K" if temp is not None else ""),
           f"CREST conformers: {len(rows)} ({t.get('source', '?')})" + (f", weights at {temp:g} K" if temp is not None else ""))]
    s.append(L("  番号  E−E(最低) [kcal/mol]  縮退度  重原子 RMSD [Å]" + ("  重み" if temp is not None else ""),
               "  index  E-E(lowest) [kcal/mol]  degeneracy  heavy-atom RMSD [Å]" + ("  weight" if temp is not None else "")))
    for r in rows:
        rel = f"{r['relative_kcal_mol']:.3f}" if r.get("relative_kcal_mol") is not None else "-"
        rmsd = f"{r['rmsd_heavy_A']:.3f}" if r.get("rmsd_heavy_A") is not None else "-"
        w = f"  {r['weight']:.4f}" if "weight" in r else ""
        s.append(f"  {r['index']}  {rel}  {r['degeneracy']}  {rmsd}{w}")
    s.append("  " + t.get("rmsd_note", ""))
    return s + ["  " + r for r in t.get("reasons", [])]


__all__ = ["CrestError", "find_crest_dir", "read_ensemble", "read_energies", "read_log_degeneracy", "boltzmann_weights",
           "analyze_crest", "summary_lines"]
