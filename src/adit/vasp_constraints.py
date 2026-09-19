"""Preserve Cartesian constraints when writing VASP lattice-direction flags."""
# Selective dynamics always refers to direct lattice vectors, including when
# positions are written in Cartesian coordinates: https://vasp.at/wiki/POSCAR
# ASE's VASP writer supports FixScaled, but omits FixCartesian:
# https://docs.ase-lib.org/_modules/ase/io/vasp.html

from __future__ import annotations

from adit.errors import AditValueError
import numpy as np
from ase.constraints import FixAtoms, FixCartesian, FixedLine, FixedPlane, FixScaled

from adit.lang import L


class VaspConstraintError(AditValueError):
    pass


def _direct_mask(cell, forbidden, dimension: int) -> np.ndarray:
    # Find lattice vectors spanning the permitted Cartesian displacement space.
    vectors = np.asarray(cell, dtype=float)
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths == 0) or not np.isfinite(vectors).all() or np.linalg.matrix_rank(vectors) != 3:
        raise VaspConstraintError(L("軸固定の変換には独立した 3 本のセルベクトルが必要です",
                                    "converting axis constraints requires three independent cell vectors"))
    permitted = np.all(np.abs((vectors / lengths[:, None]) @ np.asarray(forbidden).T) <= 1e-12, axis=1)
    if int(permitted.sum()) != dimension:
        raise VaspConstraintError(L(
            "このセルでは Cartesian 軸の固定を POSCAR の Selective dynamics で保持できません。"
            "T/F は格子ベクトル方向を指定します。固定方向を変えずに保存するには元の構造・設定を保持してください。",
            "Cartesian axis constraints cannot be preserved by POSCAR Selective dynamics for this cell. "
            "T/F flags refer to lattice-vector directions. Keep the original structure/settings to retain the fixed directions."))
    return ~permitted


def cartesian_to_direct_mask(cell, fixed) -> np.ndarray:
    """Convert a Cartesian fixed-axis mask without changing allowed motion."""
    fixed = np.asarray(fixed, dtype=bool)
    if fixed.all() or not fixed.any():
        return fixed.copy()
    return _direct_mask(cell, np.eye(3)[fixed], int((~fixed).sum()))


def prepare_vasp_constraints(atoms):
    """Return a copy with merged FixScaled masks understood by ASE's writer."""
    out = atoms.copy()
    cartesian = np.zeros((len(atoms), 3), dtype=bool)
    scaled = np.zeros_like(cartesian)
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            cartesian[constraint.index] = True
        elif isinstance(constraint, FixCartesian):
            cartesian[constraint.index] |= constraint.mask
        elif isinstance(constraint, FixScaled):
            scaled[constraint.index] |= constraint.mask
        elif not isinstance(constraint, (FixedLine, FixedPlane)):
            raise VaspConstraintError(L(
                f"POSCAR へ保持できない制約です: {type(constraint).__name__}",
                f"cannot preserve this constraint in POSCAR: {type(constraint).__name__}"))
    for i in range(len(atoms)):
        if scaled[i].all() or cartesian[i].all():
            scaled[i] = True
        else:
            scaled[i] |= cartesian_to_direct_mask(atoms.cell, cartesian[i])
    for constraint in atoms.constraints:
        if isinstance(constraint, (FixedLine, FixedPlane)):
            indices = [i for i in constraint.index if not scaled[i].all()]
            if not indices:
                continue
            direction = constraint.dir
            forbidden = (np.eye(3) - np.outer(direction, direction)
                         if isinstance(constraint, FixedLine) else [direction])
            scaled[indices] |= _direct_mask(atoms.cell, forbidden, 1 if isinstance(constraint, FixedLine) else 2)
    out.set_constraint([FixScaled(i, mask=mask) for i, mask in enumerate(scaled) if mask.any()])
    return out
