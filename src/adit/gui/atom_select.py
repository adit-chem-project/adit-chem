"""Turn atoms picked in the 3D view into the text the input fields expect."""

from __future__ import annotations

import re

from adit.lang import L

_INDEX_ONLY = re.compile(r"^\s*index\s+([\d,\-\s]+)$")


def base_indices(selected, n_base: int) -> list[int]:
    """Map displayed-atom indices (0-based, possibly of a repeated cell) to 1-based indices of the base cell, click order kept."""
    out: list[int] = []
    for i in selected:
        k = int(i) % max(1, int(n_base)) + 1
        if k not in out:
            out.append(k)
    return out


def compact_ranges(indices) -> str:
    """1-based indices -> '1-4,7' (sorted, unique)."""
    idx = sorted(set(int(i) for i in indices))
    runs: list[list[int]] = []
    for i in idx:
        if runs and runs[-1][1] == i - 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return ",".join(f"{a}" if a == b else f"{a}-{b}" for a, b in runs)


def _expand(part: str) -> list[int] | None:
    a, _, b = part.partition("-")
    try:
        lo = int(a); hi = int(b) if b else lo
    except ValueError:
        return None
    return list(range(lo, hi + 1)) if 1 <= lo <= hi else None


def merge_fixed_text(existing: str, indices) -> str:
    """Add atoms to the fixed-atoms field ('1-4,7' with optional '7:xy' parts, which are kept as they are)."""
    plain: set[int] = set(int(i) for i in indices)
    axes: list[str] = []
    for part in existing.replace(" ", "").split(","):
        if not part:
            continue
        if ":" in part:
            axes.append(part); continue
        got = _expand(part)
        if got is None:
            axes.append(part)        # unreadable text stays; validation reports it later
        else:
            plain.update(got)
    parts = ([compact_ranges(plain)] if plain else []) + axes
    return ",".join(parts)


def merge_select_text(existing: str, indices) -> str:
    """Add atoms to the selection field ('index 1,2'); other expressions are extended with 'or'."""
    new = compact_ranges(indices)
    text = existing.strip()
    if not text:
        return f"index {new}"
    m = _INDEX_ONLY.match(text)
    if m:
        old: set[int] = set()
        for part in m.group(1).replace(" ", "").split(","):
            got = _expand(part) if part else None
            if got:
                old.update(got)
        return f"index {compact_ranges(old | set(int(i) for i in indices))}"
    return f"{text} or index {new}"


SERIES_FIELD = {2: "distances", 3: "angles", 4: "dihedrals"}


def series_text(existing: str, indices) -> str:
    """Append one atom group ('1,2') to a distance/angle/dihedral series field ('1,2; 3,4')."""
    group = ",".join(str(int(i)) for i in indices)
    text = existing.strip().rstrip(";").strip()
    return f"{text}; {group}" if text else group


def destination_name(dest: str) -> str:
    return {"fixed": L("固定原子", "Fixed atoms"), "select": L("原子の選び方", "Atom selection"),
            "distances": L("距離の時系列", "Distance series"), "angles": L("角度の時系列", "Angle series"),
            "dihedrals": L("二面角の時系列", "Dihedral series")}.get(dest, dest)


__all__ = ["base_indices", "compact_ranges", "merge_fixed_text", "merge_select_text", "series_text", "SERIES_FIELD", "destination_name"]
