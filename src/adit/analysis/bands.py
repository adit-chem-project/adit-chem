
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from adit.lang import L

RY_EV = 13.605693122994


@dataclass
class BandData:
    kpts_frac: np.ndarray
    energies_ev: np.ndarray  # (nk, nbands); the spin-up channel when spin polarised
    cell: np.ndarray  # (3, 3) Å
    labels: list[tuple[int, str]] = field(default_factory=list)
    fermi_ev: float | None = None
    energies_down_ev: np.ndarray | None = None  # (nk, nbands) spin-down channel

    @property
    def n_spins(self) -> int:
        return 2 if self.energies_down_ev is not None else 1

    def all_energies(self) -> np.ndarray:
        # both channels side by side: (nk, nbands_up + nbands_down)
        if self.energies_down_ev is None:
            return self.energies_ev
        return np.concatenate([self.energies_ev, self.energies_down_ev], axis=1)

    def x_axis(self) -> np.ndarray:
        rec = 2 * np.pi * np.linalg.inv(self.cell).T
        kc = self.kpts_frac @ rec
        d = np.linalg.norm(np.diff(kc, axis=0), axis=1)
        return np.concatenate([[0.0], np.cumsum(d)])


def gap_details(energies_ev: np.ndarray, fermi_ev: float | None, kpts_frac=None, labels=None,
                energies_down_ev: np.ndarray | None = None) -> dict | None:
    if fermi_ev is None or np.size(energies_ev) == 0:
        return None
    e = np.asarray(energies_ev, dtype=float)
    if energies_down_ev is not None:
        e = np.concatenate([e, np.asarray(energies_down_ev, dtype=float)], axis=1)
    occupied = e <= fermi_ev
    if not occupied.any() or occupied.all():
        return None
    crossing = occupied.any(axis=0) & ~occupied.all(axis=0)
    if crossing.any():
        return {"gap_ev": 0.0, "metal": True, "n_bands_crossing": int(crossing.sum()),
                "note": L("フェルミ準位を横切るバンドがあるので、ギャップは 0 (金属) です。VBM と CBM は定義しません",
                          "a band crosses the Fermi level, so the gap is 0 (metal); no VBM or CBM is defined")}
    vb = np.where(occupied, e, -np.inf)
    cb = np.where(~occupied, e, np.inf)
    vbm_k = int(np.argmax(vb.max(axis=1)))
    cbm_k = int(np.argmin(cb.min(axis=1)))
    vbm = float(vb[vbm_k].max())
    cbm = float(cb[cbm_k].min())
    direct_per_k = cb.min(axis=1) - vb.max(axis=1)
    direct_k = int(np.argmin(direct_per_k))
    name = dict(labels or {})
    out = {"gap_ev": cbm - vbm, "metal": False, "vbm_ev": vbm, "cbm_ev": cbm,
           "vbm_kpoint_index": vbm_k, "cbm_kpoint_index": cbm_k,
           "vbm_label": name.get(vbm_k, ""), "cbm_label": name.get(cbm_k, ""),
           "direct": vbm_k == cbm_k,
           "direct_gap_ev": float(direct_per_k[direct_k]), "direct_gap_kpoint_index": direct_k,
           "direct_gap_label": name.get(direct_k, ""),
           "note": L("VBM と CBM が同じ k 点にあるかどうかを「直接」と書いています (機械的な事実)。"
                     "k 点は与えた経路の中での番号で、経路の外により低い CBM があるかどうかは分かりません。",
                     "'direct' only records whether the VBM and CBM are at the same k-point. The k-point index refers to "
                     "the given path; whether a lower CBM exists outside the path is unknown.")}
    if energies_down_ev is not None:
        out["spin_channels"] = 2
    if kpts_frac is not None:
        k = np.asarray(kpts_frac, dtype=float)
        out["vbm_kpoint_frac"] = k[vbm_k].tolist()
        out["cbm_kpoint_frac"] = k[cbm_k].tolist()
    return out


def load_bands(run_dir: Path, code: str, fermi_ev: float | None) -> BandData | None:
    bd = run_dir / "bands"
    if not (bd / "kpath.json").is_file():
        return None
    kp = json.loads((bd / "kpath.json").read_text(encoding="utf-8"))
    kpts = np.array(kp["kpts"]); labels = [(int(i), str(l)) for i, l in kp["labels"]]; cell = np.array(kp["cell"])
    if code == "dftbplus":
        up, down = _dftb_band_out(bd / "band.out", len(kpts))
    elif code == "espresso":
        up, down = _qe_bands_out(bd / "output.log")
    elif code == "vasp":
        up, down = _vasp_eigenval(bd / "EIGENVAL")
    else:
        return None
    if not up:
        return None
    n = min(len(up), len(kpts))
    e_down = np.array(down[:n]) if down and len(down) >= n else None
    return BandData(kpts_frac=kpts[:n], energies_ev=np.array(up[:n]), cell=cell, labels=[(i, l) for i, l in labels if i < n],
                    fermi_ev=fermi_ev, energies_down_ev=e_down)


def _dftb_band_out(p: Path, nk: int | None = None):
    if not p.is_file():
        return None, None
    blocks, spins, cur = [], [], None
    for l in p.read_text(encoding="utf-8", errors="replace").splitlines():
        s = l.strip()
        if s.startswith("KPT"):
            cur = []; blocks.append(cur)
            m = re.search(r"SPIN\s+(\d+)", s)
            spins.append(int(m.group(1)) if m else 1)
        elif cur is not None:
            w = l.split()
            if len(w) == 3:
                cur.append(float(w[1]))
    up = [b for s, b in zip(spins, blocks) if b and s == 1]
    down = [b for s, b in zip(spins, blocks) if b and s == 2]
    if not down and nk and len(up) == 2 * nk:
        # spin is the outer loop in band.out: the second half of the blocks is the other channel
        up, down = up[:nk], up[nk:]
    return up, (down or None)


def _qe_bands_out(p: Path):
    if not p.is_file():
        return None, None
    txt = p.read_text(encoding="utf-8", errors="replace")
    i_up, i_down = txt.find("SPIN UP"), txt.find("SPIN DOWN")
    if 0 <= i_up < i_down:
        return _qe_blocks(txt[i_up:i_down]), _qe_blocks(txt[i_down:])
    return _qe_blocks(txt), None


def _qe_blocks(txt: str):
    parts = re.split(r"\n\s+k =.*?bands \(ev\):\s*\n", txt)[1:]
    out = []
    for part in parts:
        vals = []
        for l in part.splitlines():
            if not l.strip():
                if vals:
                    break
                continue
            if l.strip().startswith(("occupation", "highest", "Writing", "the Fermi")):
                break
            try:
                vals += [float(x) for x in re.findall(r"-?\d+\.\d+", l)]
            except ValueError:
                break
        if vals:
            out.append(vals)
    return out


def _vasp_eigenval(p: Path):
    if not p.is_file():
        return None, None
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        _, nk, nb = [int(x) for x in lines[5].split()[:3]]
    except (ValueError, IndexError):
        return None, None
    up, down, i = [], [], 6
    for _ in range(nk):
        while i < len(lines) and not lines[i].strip():
            i += 1
        i += 1
        vals, vals_down = [], []
        for _b in range(nb):
            w = lines[i].split(); i += 1
            vals.append(float(w[1]))
            if len(w) >= 5:  # ISPIN = 2: band, E(up), E(down), occ(up), occ(down)
                vals_down.append(float(w[2]))
        up.append(vals)
        if vals_down:
            down.append(vals_down)
    return up, (down if len(down) == len(up) else None)


def plot_bands(bd: BandData, path_png: Path, window_ev: float = 10.0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = bd.x_axis()
    channels = [(bd.energies_ev, "C0", "-", "spin up")]
    if bd.energies_down_ev is not None:
        channels.append((bd.energies_down_ev, "C1", "--", "spin down"))
    fig, ax = plt.subplots(figsize=(6, 4))
    for energies, color, ls, label in channels:
        e = energies - (bd.fermi_ev or 0.0)
        for b in range(e.shape[1]):
            ax.plot(x, e[:, b], lw=1.0, color=color, ls=ls, label=label if (b == 0 and len(channels) > 1) else None)
    if len(channels) > 1:
        ax.legend(fontsize=8, loc="upper right")
    for i, l in bd.labels:
        ax.axvline(x[i], color="gray", lw=0.6)
    ax.set_xticks([x[i] for i, _ in bd.labels]); ax.set_xticklabels([l.replace("G", "Γ") for _, l in bd.labels])
    if bd.fermi_ev is not None:
        ax.axhline(0, color="gray", ls="--", lw=0.8); ax.set_ylabel("E − E$_F$ [eV]")
        ax.set_ylim(-window_ev, window_ev)
    else:
        ax.set_ylabel("E [eV]")
    ax.set_xlim(x[0], x[-1]); ax.grid(alpha=0.2)
    fig.savefig(path_png, dpi=110, bbox_inches="tight"); plt.close(fig)
