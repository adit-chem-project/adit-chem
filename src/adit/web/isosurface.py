"""Isosurface triangles of a volumetric file for the browser, bounded in count."""

from __future__ import annotations

import json
from pathlib import Path

from adit.lang import L
from adit.web.structure3d import scene_from_atoms

MAX_WEB_TRIANGLES = 60_000


def volumetric_files(run_dir: Path | str | None) -> list[str]:
    from adit.analysis.volumetric import find_files

    return [p.name for p in find_files(run_dir)] if run_dir else []


def isosurface_payload(run_dir: Path | str, name: str, level: str | float | None = None, stride: int = 1,
                       max_triangles: int = MAX_WEB_TRIANGLES) -> dict:
    """Stats of the file (always) and, when a level is given, the surfaces at +level / -level plus the atoms as a scene."""
    from adit.analysis.isosurface import IsosurfaceError, IsosurfaceTooLarge, grid_stats, isosurface, isosurfaces, payload
    from adit.analysis.volumetric import find_files, read_grid

    d = Path(run_dir)
    allowed = {p.name: p for p in find_files(d)}
    if name not in allowed:
        raise IsosurfaceError(L(f"このディレクトリにその体積データはありません: {name}", f"no such volumetric file in this directory: {name}"))
    grid = read_grid(allowed[name])
    st = grid_stats(grid)
    out: dict = {"file": name, "stats": st, "surfaces": [], "notes": [], "suggested_stride": None}
    atoms = grid.atoms.copy()
    if not grid.periodic:
        atoms.pbc = False
    out["scene"] = json.loads(scene_from_atoms(atoms).to_json())
    out["center"] = [round(float(x), 4) for x in atoms.get_positions().mean(axis=0)]
    if level is None or (isinstance(level, str) and not level.strip()):
        return out
    try:
        lv = float(str(level).replace(",", ""))
    except ValueError as ex:
        raise IsosurfaceError(L(f"等値を数で入れてください: {level!r}", f"the level must be a number: {level!r}")) from ex
    stride = max(1, int(stride))
    try:
        meshes = isosurfaces(grid, lv, stride, max_triangles) if lv > 0 else [isosurface(grid, lv, stride, max_triangles)]
    except IsosurfaceTooLarge as ex:
        out["notes"].append(str(ex)); out["suggested_stride"] = int(ex.suggested_stride)
        return out
    out["surfaces"] = payload(meshes)
    out["level"] = lv; out["stride"] = stride
    for m in meshes:
        for n in m.notes:
            if n not in out["notes"]:
                out["notes"].append(n)
    if any(m.n_triangles == 0 for m in meshes):
        out["notes"].append(L("三角形が 0 枚の等値があります (その値を通る場所が格子にありません)",
                              "a level produced no triangles (no grid values cross it)"))
    return out


def isosurface_json(run_dir: Path | str, name: str, **kw) -> str:
    return json.dumps(isosurface_payload(run_dir, name, **kw), separators=(",", ":")).replace("<", "\\u003c")


__all__ = ["MAX_WEB_TRIANGLES", "volumetric_files", "isosurface_payload", "isosurface_json"]
