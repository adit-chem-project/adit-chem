"""Read Bader charges from the ACF.dat written by the Henkelman bader program. ADIT does not partition charges itself."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from adit.lang import L


def run_note() -> list[str]:
    return [
        L("bader は ADIT に入っていません。Henkelman グループの bader (GPL) を別に入れてください:",
          "bader is not part of ADIT. Install the Henkelman group's bader (GPL) separately:"),
        "  https://theory.cm.utexas.edu/henkelman/code/bader/",
        "  VASP:  chgsum.pl AECCAR0 AECCAR2 && bader CHGCAR -ref CHGCAR_sum",
        "  cube:  bader density.cube",
        L("できた ACF.dat を adit-analyze --bader ACF.dat に渡してください。",
          "Pass the resulting ACF.dat to adit-analyze --bader ACF.dat."),
    ]


def __getattr__(name: str):
    # RUN_NOTE stays importable; the language is read on each access.
    if name == "RUN_NOTE":
        return run_note()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


@dataclass
class BaderCharges:
    positions_bohr: np.ndarray
    electrons: np.ndarray
    min_dist: np.ndarray
    volume: np.ndarray
    vacuum_charge: float | None = None
    total_electrons: float | None = None

    def net_charge(self, valence_electrons) -> np.ndarray:
        z = np.asarray(valence_electrons, dtype=float)
        if z.shape != self.electrons.shape:
            raise ValueError(L(f"価電子数の個数 {z.size} が原子数 {self.electrons.size} と違います",
                               f"{z.size} valence electron counts do not match {self.electrons.size} atoms"))
        return z - self.electrons


def read_acf(path: Path) -> BaderCharges:
    text = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    rows, vacuum, total = [], None, None
    for line in text:
        s = line.strip()
        if s.upper().startswith("VACUUM CHARGE"):
            vacuum = float(s.split(":")[1])
        elif s.upper().startswith("NUMBER OF ELECTRONS"):
            total = float(s.split(":")[1])
        parts = s.split()
        if len(parts) == 7:
            try:
                values = [float(x) for x in parts]
            except ValueError:
                continue
            if values[0] == int(values[0]) and int(values[0]) == len(rows) + 1:
                rows.append(values[1:])
    if not rows:
        raise ValueError(L(f"{Path(path).name} に Bader の表 (番号 X Y Z 電荷 最小距離 体積) がありません",
                           f"{Path(path).name} has no Bader table (index X Y Z charge min. distance volume)"))
    data = np.array(rows)
    return BaderCharges(data[:, 0:3], data[:, 3], data[:, 4], data[:, 5], vacuum, total)


def summary_lines(bc: BaderCharges, symbols=None, valence_electrons=None) -> list[str]:
    lines = [L(f"Bader 電荷 (外部の bader プログラムの ACF.dat、原子 {bc.electrons.size} 個)",
               f"Bader charges (ACF.dat from the external bader program, {bc.electrons.size} atoms)")]
    net = bc.net_charge(valence_electrons) if valence_electrons is not None else None
    for i, e in enumerate(bc.electrons):
        name = f"{symbols[i]}" if symbols is not None and i < len(symbols) else L("原子", "atom")
        diff_ja = f"、価電子との差 {net[i]:+.4f}" if net is not None else ""
        diff_en = f", difference from the valence count {net[i]:+.4f}" if net is not None else ""
        lines.append(L(f"  {name} {i + 1}: 電子 {e:.4f}{diff_ja} (最近接面まで {bc.min_dist[i]:.3f} Bohr)",
                       f"  {name} {i + 1}: {e:.4f} electrons{diff_en} ({bc.min_dist[i]:.3f} Bohr to the nearest surface)"))
    if bc.total_electrons is not None:
        lines.append(L(f"  合計 {bc.total_electrons:.4f} 電子"
                       + (f"、真空 {bc.vacuum_charge:.4f}" if bc.vacuum_charge is not None else ""),
                       f"  total {bc.total_electrons:.4f} electrons"
                       + (f", vacuum {bc.vacuum_charge:.4f}" if bc.vacuum_charge is not None else "")))
    lines.append(L("  最近接面までの距離が小さいときは、電荷密度の格子が粗い可能性があります。",
                   "  A small distance to the nearest surface may mean the charge-density grid is too coarse."))
    if valence_electrons is None:
        lines.append(L("  価電子数を渡すと、価電子との差 (いわゆる Bader 電荷) も出します。",
                       "  Pass the valence electron counts to also get the difference from them (the usual Bader charge)."))
    return lines
