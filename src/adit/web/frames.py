"""Frames of a run directory for the browser player, thinned to a bounded count."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from adit.lang import L
from adit.web.structure3d import scene_from_atoms

MAX_WEB_FRAMES = 200


def frames_payload(run_dir: Path | str, *, stride: int = 1, skip: int = 0, max_frames: int = MAX_WEB_FRAMES,
                   memory_mb: float | None = None) -> dict:
    """Scene of the first frame plus positions per frame; energies/temperatures/times aligned to the frames when possible."""
    from adit.analysis.readers import load_run
    from adit.analysis.trajectory import MEMORY_BUDGET_MB, Trajectory, TrajectoryTooLarge, check_budget

    d = Path(run_dir)
    out: dict = {"n_frames": 0, "n_total": 0, "stride": 1, "notes": []}
    try:
        data = load_run(d)
    except Exception as ex:
        out["notes"].append(L(f"再生できません: {ex}", f"cannot play: {ex}"))
        return out
    n_total = len(data.frames)
    out["n_total"] = n_total
    if n_total < 2:
        out["notes"].append(L("フレームが 1 つしかないので軌跡は再生できません", "only one frame, so there is no trajectory to play"))
        return out
    stride = max(1, int(stride)); skip = int(skip) if n_total > skip else 0
    remaining = max(1, n_total - skip)
    if math.ceil(remaining / stride) > max_frames:
        stride = math.ceil(remaining / max_frames)
        out["notes"].append(L(f"ブラウザに送るフレームは最大 {max_frames} 個なので、{stride} フレームに 1 回に間引いています",
                              f"at most {max_frames} frames are sent to the browser, so every {stride}th frame is used"))
    sel = data.frames[skip::stride]
    nat = data.frames.first_natoms() if isinstance(data.frames, Trajectory) else len(data.frames[0])
    try:
        check_budget(sel, nat, memory_mb or MEMORY_BUDGET_MB, what=L("再生", "playback"))
    except TrajectoryTooLarge as ex:
        stride = int(ex.suggested_stride)
        sel = data.frames[skip::stride]
        out["notes"].append(L(f"軌跡が大きいので {stride} フレームに 1 回に間引いて読みました", f"the trajectory is large, so every {stride}th frame was read"))
    frames = list(sel)
    if not frames:
        return out
    idx = np.arange(skip, skip + stride * len(frames), stride)[: len(frames)]
    scene = scene_from_atoms(frames[0])
    center = np.asarray(frames[0].get_positions(), dtype=float).mean(axis=0)
    out.update(scene=json.loads(scene.to_json()), stride=int(stride), skip=int(skip), n_frames=len(frames),
               positions=[np.round(np.asarray(f.get_positions(), dtype=float) - center, 3).tolist() for f in frames],
               frame_index=[int(i) for i in idx])
    times = None
    if data.times_fs and len(data.times_fs) == n_total:
        times = np.asarray(data.times_fs, dtype=float)[idx]
    elif data.frame_dt_fs:
        times = idx * float(data.frame_dt_fs)
    if times is not None:
        out["times_fs"] = [round(float(t), 3) for t in times]
    e = np.asarray(data.energies_ev, dtype=float) if data.energies_ev else np.zeros(0)
    if len(e) == n_total:
        out["energies_ev"] = [float(x) for x in e[idx]]
        if data.temperatures_k and len(data.temperatures_k) == n_total:
            out["temperatures_k"] = [float(x) for x in np.asarray(data.temperatures_k, dtype=float)[idx]]
    elif len(e) >= 2:
        x = np.linspace(0, len(e) - 1, len(frames))
        out["energies_ev"] = [float(v) for v in np.interp(x, np.arange(len(e)), e)]
        out["notes"].append(L(f"エネルギーの点数 ({len(e)}) とフレーム数 ({n_total}) が違うので、縦線の位置は比例で合わせています",
                              f"the energy list ({len(e)} points) and the frames ({n_total}) differ in count, so the marker is placed proportionally"))
    out["is_md"] = bool(data.temperatures_k)
    out["source"] = data.frame_source or ""
    return out


def frames_json(run_dir: Path | str, **kw) -> str:
    return json.dumps(frames_payload(run_dir, **kw), separators=(",", ":")).replace("<", "\\u003c")


__all__ = ["frames_payload", "frames_json", "MAX_WEB_FRAMES"]
