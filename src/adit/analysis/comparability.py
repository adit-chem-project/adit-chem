"""Cross-code comparability checks for native output analysis."""
# This module reads existing run directories and checks mechanical facts that
# must line up before comparing MSD or energy across different engines. It does
# not judge whether a setup is chemically meaningful.

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms

from adit.analysis.readers import RunData, load_run
from adit.lang import L
from adit.project import load_project


DEFAULT_CELL_ATOL_ANG = 1e-6
DEFAULT_TIME_ATOL_FS = 1e-6


@dataclass
class AuditIssue:
    severity: str
    category: str
    message: str
    run: str | None = None

    def to_dict(self) -> dict:
        return {"severity": self.severity, "category": self.category, "run": self.run, "message": self.message}


@dataclass
class RunAudit:
    label: str
    path: str
    code: str | None = None
    spec_code: str | None = None
    task: str | None = None
    natoms: int | None = None
    symbols: list[str] = field(default_factory=list)
    composition: dict[str, int] = field(default_factory=dict)
    pbc: tuple[bool, bool, bool] | None = None
    cell_ang: list[list[float]] | None = None
    frame_source: str = ""
    n_frames: int = 0
    frame_dt_fs: float | None = None
    energy_points: int = 0
    final_energy_ev: float | None = None
    temperature_points: int | None = None
    time_points: int | None = None
    time_unit: str = "fs"
    length_unit: str = "angstrom"
    energy_unit: str = "eV per simulated system"
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "path": self.path,
            "code": self.code,
            "spec_code": self.spec_code,
            "task": self.task,
            "natoms": self.natoms,
            "symbols": self.symbols,
            "composition": self.composition,
            "pbc": list(self.pbc) if self.pbc is not None else None,
            "cell_ang": self.cell_ang,
            "frame_source": self.frame_source,
            "n_frames": self.n_frames,
            "frame_dt_fs": self.frame_dt_fs,
            "energy_points": self.energy_points,
            "final_energy_ev": self.final_energy_ev,
            "temperature_points": self.temperature_points,
            "time_points": self.time_points,
            "units": {"time": self.time_unit, "length": self.length_unit, "energy": self.energy_unit},
            "notes": self.notes,
        }


@dataclass
class ComparabilityAudit:
    runs: list[RunAudit]
    issues: list[AuditIssue] = field(default_factory=list)
    reference: str | None = None
    mode: str = "msd"

    @property
    def comparable(self) -> bool:
        return not any(i.severity == "error" for i in self.issues)

    def to_dict(self) -> dict:
        return {
            "comparable": self.comparable,
            "mode": self.mode,
            "reference": self.reference,
            "runs": [r.to_dict() for r in self.runs],
            "issues": [i.to_dict() for i in self.issues],
        }

    def summary_text(self) -> str:
        head = L(
            f"コードをまたぐ比較の点検 ({self.mode}): {len(self.runs)} 計算、結果: {'問題なし' if self.comparable else '問題あり'}",
            f"cross-code check ({self.mode}): {len(self.runs)} runs, result: {'passed' if self.comparable else 'failed'}",
        )
        lines = [head]
        for run in self.runs:
            bits = [
                run.code or "?",
                run.task or "?",
                L(f"{run.natoms if run.natoms is not None else '?'} 原子", f"{run.natoms if run.natoms is not None else '?'} atoms"),
                L(f"{run.n_frames} フレーム", f"{run.n_frames} frames"),
            ]
            if run.frame_dt_fs is not None:
                bits.append(L(f"1 フレーム {run.frame_dt_fs:g} fs", f"{run.frame_dt_fs:g} fs per frame"))
            if run.final_energy_ev is not None:
                bits.append(f"E = {run.final_energy_ev:.8g} eV")
            lines.append(f"  {run.label}: " + " / ".join(bits))
        for issue in self.issues:
            where = f"{issue.run}: " if issue.run else ""
            level = L("エラー", "error") if issue.severity == "error" else L("注意", "warning")
            lines.append(f"{level}: {where}{issue.message}")
        lines.append(L("注: 化学的に比べてよいかは判断していません", "note: chemical comparability is not judged"))
        return "\n".join(lines)


def audit_run_dirs(
    run_dirs: Sequence[Path | str],
    *,
    labels: Sequence[str] | None = None,
    require_md: bool = True,
    cell_atol_ang: float = DEFAULT_CELL_ATOL_ANG,
    time_atol_fs: float = DEFAULT_TIME_ATOL_FS,
) -> ComparabilityAudit:
    """Read existing run directories and audit mechanical comparability."""
    # The returned object is read-only analysis data: no files are written. Values
    # taken from native outputs are normalized by ``load_run``: coordinates in
    # angstrom, times in fs, and energies in eV for the simulated system.
    # ``require_md=True`` audits MSD comparability: trajectories, element-symbol
    # order, PBC, cells, and frame spacing must match, but energy is optional.
    # ``False`` audits energy comparability: composition, PBC, cells, and
    # whole-system energies must be readable; trajectories and atom
    # order are not blocking.
    paths = [Path(p).expanduser() for p in run_dirs]
    names = list(labels) if labels is not None else [p.name or str(p) for p in paths]
    if len(names) != len(paths):
        raise ValueError(L("labels の数が run_dirs と違います", "labels must have the same length as run_dirs"))
    if len(set(names)) != len(names):
        if labels is not None:
            raise ValueError(L("比較する計算のラベルは重複できません", "run labels must be unique"))
        counts: dict[str, int] = {}
        unique = []
        for name in names:
            counts[name] = counts.get(name, 0) + 1
            unique.append(f"{name}#{counts[name]}")
        names = unique

    runs: list[RunAudit] = []
    issues: list[AuditIssue] = []
    atoms_by_label: dict[str, Atoms] = {}

    for path, label in zip(paths, names):
        run_audit, atoms, local_issues = _audit_one(path, label, require_md=require_md, time_atol_fs=time_atol_fs)
        runs.append(run_audit)
        issues.extend(local_issues)
        if atoms is not None:
            atoms_by_label[label] = atoms

    if len(runs) >= 2:
        issues.extend(_cross_check(runs, atoms_by_label, require_md=require_md, cell_atol_ang=cell_atol_ang, time_atol_fs=time_atol_fs))

    return ComparabilityAudit(runs=runs, issues=issues, reference=runs[0].label if runs else None,
                              mode="msd" if require_md else "energy")


def _audit_one(path: Path, label: str, *, require_md: bool, time_atol_fs: float) -> tuple[RunAudit, Atoms | None, list[AuditIssue]]:
    out = RunAudit(label=label, path=str(path))
    issues: list[AuditIssue] = []
    spec = None
    run: RunData | None = None

    try:
        spec = load_project(path)
        out.spec_code = spec.method.code
        out.task = spec.task.type
    except Exception as ex:
        issues.append(AuditIssue("warning", "spec.missing", L(f"spec.json を読めません: {ex}", f"cannot read spec.json: {ex}"), label))

    try:
        run = load_run(path)
        out.code = run.code
        out.frame_source = run.frame_source
        out.frame_dt_fs = run.frame_dt_fs
        out.energy_points = len(run.energies_ev)
        if len(run.energies_ev):
            final_energy = float(run.energies_ev[-1])
            if np.isfinite(final_energy):
                out.final_energy_ev = final_energy
            else:
                issues.append(AuditIssue("error", "energy.nonfinite", L(
                    "最終エネルギーが有限の数ではありません",
                    "the final energy is not finite"), label))
        out.temperature_points = None if run.temperatures_k is None else len(run.temperatures_k)
        out.time_points = None if run.times_fs is None else len(run.times_fs)
        out.notes.extend(run.notes)
    except Exception as ex:
        issues.append(AuditIssue("error", "run.unreadable", L(f"出力を読めません: {ex}", f"cannot read native output: {ex}"), label))

    if out.spec_code and out.code and out.spec_code != out.code:
        issues.append(AuditIssue("error", "code.mismatch", L(f"spec.json は {out.spec_code}、出力は {out.code} です",
                                                             f"spec.json says {out.spec_code}, but native output reads as {out.code}"), label))
    if require_md and out.task and out.task != "molecular_dynamics":
        issues.append(AuditIssue("error", "task.not_md", L(f"MSD 比較用ですが、計算種別は {out.task} です",
                                                           f"MSD comparison was requested, but the task is {out.task}"), label))

    atoms = _first_atoms(run)
    if atoms is None and spec is not None and hasattr(spec, "structure"):
        atoms = spec.structure.atoms.to_ase()
        issues.append(AuditIssue("warning", "trajectory.fallback_spec", L("軌跡が読めないため spec.json の構造で原子列を確認しました",
                                                                           "trajectory was not readable, so atom identity was checked from spec.json"), label))

    if atoms is not None:
        _fill_atoms(out, atoms)
    else:
        issues.append(AuditIssue("error", "atoms.unavailable", L(
            "構造を読めないため、比較する系の原子を確認できません",
            "structure is unavailable, so the atoms in the compared system cannot be checked"), label))
    if run is not None:
        out.n_frames = _safe_len(run.frames)
        if require_md and not run.frames:
            issues.append(AuditIssue("error", "trajectory.missing", L("軌跡または最終構造がありません", "no trajectory or final structure was read"), label))
        if require_md and out.n_frames < 2:
            issues.append(AuditIssue("error", "trajectory.too_short", L(f"MSD には 2 フレーム以上が必要ですが、読めたのは {out.n_frames} フレームです",
                                                                         f"MSD needs at least 2 frames, but only {out.n_frames} were read"), label))
        if require_md and out.frame_dt_fs is None:
            issues.append(AuditIssue("error", "trajectory.dt_missing", L("軌跡フレーム間隔 [fs] が分かりません",
                                                                         "trajectory frame spacing in fs is unknown"), label))
        if not require_md and out.energy_points == 0:
            issues.append(AuditIssue("error", "energy.missing", L("比較用のエネルギーが読めません", "no comparable energy was read"), label))
        if not require_md and out.energy_points and _energy_is_per_atom(run):
            issues.append(AuditIssue("error", "energy.per_atom", L("エネルギーが原子 1 個あたりの値として読まれています。全系の eV と同じ表では比べません",
                                                                  "energy was read as a per-atom value; it is not compared in the same table as whole-system eV"), label))
        _check_time_series(out, run, issues, time_atol_fs)

    return out, atoms, issues


def _first_atoms(run: RunData | None) -> Atoms | None:
    if run is None or not run.frames:
        return None
    try:
        return run.frames[0]
    except Exception:
        return None


def _safe_len(frames: Sequence[Atoms]) -> int:
    try:
        return len(frames)
    except Exception:
        return 0


def _energy_is_per_atom(run: RunData) -> bool:
    text = "\n".join(run.notes).lower()
    return "per atom" in text or "原子 1 個あたり" in text


def _fill_atoms(out: RunAudit, atoms: Atoms) -> None:
    symbols = list(atoms.get_chemical_symbols())
    out.symbols = symbols
    out.natoms = len(symbols)
    out.composition = dict(Counter(symbols))
    out.pbc = tuple(bool(x) for x in atoms.pbc)
    if any(out.pbc):
        out.cell_ang = np.asarray(atoms.cell, dtype=float).tolist()


def _check_time_series(out: RunAudit, run: RunData, issues: list[AuditIssue], time_atol_fs: float) -> None:
    if run.times_fs is not None and out.energy_points and len(run.times_fs) != out.energy_points:
        issues.append(AuditIssue("warning", "times.energy_length", L(f"時刻 {len(run.times_fs)} 点とエネルギー {out.energy_points} 点の数が違います",
                                                                     f"time points ({len(run.times_fs)}) and energy points ({out.energy_points}) differ"), out.label))
    if run.temperatures_k is not None and run.times_fs is not None and len(run.temperatures_k) != len(run.times_fs):
        issues.append(AuditIssue("warning", "times.temperature_length", L(f"時刻 {len(run.times_fs)} 点と温度 {len(run.temperatures_k)} 点の数が違います",
                                                                          f"time points ({len(run.times_fs)}) and temperature points ({len(run.temperatures_k)}) differ"), out.label))
    if run.times_fs is not None and out.frame_dt_fs is not None and out.n_frames == len(run.times_fs) and len(run.times_fs) >= 2:
        diffs = np.diff(np.asarray(run.times_fs, dtype=float))
        if np.isfinite(diffs).all() and float(np.max(np.abs(diffs - out.frame_dt_fs))) > time_atol_fs:
            issues.append(AuditIssue("warning", "times.frame_dt_mismatch", L("時刻列の間隔と軌跡フレーム間隔が一致しません",
                                                                              "time-series spacing does not match trajectory frame spacing"), out.label))


def _cross_check(runs: list[RunAudit], atoms_by_label: dict[str, Atoms], *, require_md: bool, cell_atol_ang: float, time_atol_fs: float) -> list[AuditIssue]:
    issues: list[AuditIssue] = []
    ref = runs[0]
    ref_atoms = atoms_by_label.get(ref.label)
    ref_symbols = ref.symbols
    ref_counts = Counter(ref_symbols)
    ref_pbc = ref.pbc
    ref_cell = np.asarray(ref.cell_ang, dtype=float) if ref.cell_ang is not None else None

    for run in runs[1:]:
        if not require_md and ref.code and run.code and ref.code != run.code:
            issues.append(AuditIssue("error", "energy.definition_unverified", L(
                f"{ref.code} と {run.code} が出力するエネルギーの内訳 (運動・ポテンシャル・電子の寄与など)を照合できないため、直接比較の点検は通過できません。",
                f"The energy contributions reported by {ref.code} and {run.code} (such as kinetic, potential, or electronic terms) have not been matched, so the direct-comparison check cannot pass."), run.label))
        if ref_symbols and run.symbols:
            counts = Counter(run.symbols)
            if counts != ref_counts:
                issues.append(AuditIssue("error", "atoms.composition", L(f"{run.label} の組成 {dict(counts)} が基準 {ref.label} の {dict(ref_counts)} と違います",
                                                                          f"{run.label} composition {dict(counts)} differs from reference {ref.label} {dict(ref_counts)}"), run.label))
            elif require_md and run.symbols != ref_symbols:
                issues.append(AuditIssue("error", "atoms.order", L(f"{run.label} の原子順が基準 {ref.label} と違います",
                                                                    f"{run.label} atom order differs from reference {ref.label}"), run.label))
        if ref_pbc is not None and run.pbc is not None and run.pbc != ref_pbc:
            issues.append(AuditIssue("error", "pbc.mismatch", L(f"{run.label} の周期境界 {run.pbc} が基準 {ref.label} の {ref_pbc} と違います",
                                                                f"{run.label} PBC {run.pbc} differs from reference {ref.label} {ref_pbc}"), run.label))
        if ref_cell is not None or run.cell_ang is not None:
            _check_cell(ref, ref_cell, run, issues, cell_atol_ang)
        if require_md and ref.n_frames and run.n_frames and run.n_frames != ref.n_frames:
            issues.append(AuditIssue("error", "trajectory.frame_count", L(f"{run.label} のフレーム数 {run.n_frames} が基準 {ref.label} の {ref.n_frames} と違います",
                                                                          f"{run.label} frame count {run.n_frames} differs from reference {ref.label} {ref.n_frames}"), run.label))
        if require_md and ref.frame_dt_fs is not None and run.frame_dt_fs is not None and abs(run.frame_dt_fs - ref.frame_dt_fs) > time_atol_fs:
            issues.append(AuditIssue("error", "trajectory.frame_dt", L(f"{run.label} のフレーム間隔 {run.frame_dt_fs:g} fs が基準 {ref.label} の {ref.frame_dt_fs:g} fs と違います",
                                                                       f"{run.label} frame spacing {run.frame_dt_fs:g} fs differs from reference {ref.label} {ref.frame_dt_fs:g} fs"), run.label))
        if ref.energy_points and run.energy_points and run.energy_points != ref.energy_points:
            issues.append(AuditIssue("warning", "energy.point_count", L(f"{run.label} のエネルギー点数 {run.energy_points} が基準 {ref.label} の {ref.energy_points} と違います",
                                                                        f"{run.label} energy point count {run.energy_points} differs from reference {ref.label} {ref.energy_points}"), run.label))
        if require_md and ref_atoms is not None and run.label in atoms_by_label:
            _check_identity_positions(ref, ref_atoms, run, atoms_by_label[run.label], issues)
    return issues


def _check_cell(ref: RunAudit, ref_cell: np.ndarray | None, run: RunAudit, issues: list[AuditIssue], cell_atol_ang: float) -> None:
    if ref_cell is None or run.cell_ang is None:
        issues.append(AuditIssue("error", "cell.missing", L(f"{ref.label} と {run.label} の片方だけセルがあります",
                                                            f"only one of {ref.label} and {run.label} has a cell"), run.label))
        return
    cell = np.asarray(run.cell_ang, dtype=float)
    if ref_cell.shape != (3, 3) or cell.shape != (3, 3) or not np.isfinite(ref_cell).all() or not np.isfinite(cell).all():
        issues.append(AuditIssue("error", "cell.invalid", L("セルが 3x3 の有限値ではありません", "cell is not a finite 3x3 matrix"), run.label))
        return
    delta = float(np.max(np.abs(cell - ref_cell)))
    if delta > cell_atol_ang:
        issues.append(AuditIssue("error", "cell.mismatch", L(f"{run.label} のセルが基準 {ref.label} と違います (最大差 {delta:.4g} Å)",
                                                             f"{run.label} cell differs from reference {ref.label} (max difference {delta:.4g} Å)"), run.label))


def _check_identity_positions(ref: RunAudit, ref_atoms: Atoms, run: RunAudit, atoms: Atoms, issues: list[AuditIssue]) -> None:
    if len(ref_atoms) != len(atoms):
        return
    if ref.symbols != run.symbols:
        return
    dr = np.asarray(atoms.positions, dtype=float) - np.asarray(ref_atoms.positions, dtype=float)
    if not np.isfinite(dr).all():
        issues.append(AuditIssue("error", "atoms.positions_invalid", L("初期座標に有限でない値があります", "initial coordinates contain non-finite values"), run.label))
        return
    if any(ref_atoms.pbc) and any(atoms.pbc):
        from ase.geometry import find_mic
        try:
            _, distances = find_mic(dr, ref_atoms.cell, pbc=ref_atoms.pbc)
        except (ValueError, np.linalg.LinAlgError):
            issues.append(AuditIssue("warning", "atoms.identity_unverifiable", L(
                "周期境界で初期座標を照合できず、同じ元素の原子の対応を確認できません",
                "initial coordinates could not be matched across periodic boundaries, so identities of same-element atoms remain unverified"), run.label))
            return
    else:
        distances = np.linalg.norm(dr, axis=1)
    if distances.size and float(np.max(distances)) > 1e-4:
        issues.append(AuditIssue("error", "atoms.initial_positions", L(
            f"{run.label} の初期座標が基準 {ref.label} と異なります。同じ初期構造からの MSD 比較としては点検を通過できません。",
            f"{run.label} initial coordinates differ from reference {ref.label}; the check cannot pass for an MSD comparison starting from the same structure."), run.label))


__all__ = ["AuditIssue", "RunAudit", "ComparabilityAudit", "audit_run_dirs"]
