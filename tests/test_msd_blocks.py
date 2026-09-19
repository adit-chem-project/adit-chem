"""Block diffusion estimates retain frames and the main absolute lag window."""

import numpy as np
import pytest

from adit.analysis import compute


def test_blocks_use_every_frame_and_same_absolute_fit_window():
    # Piecewise velocities make omission or shifting of the last frame observable.
    pos = (np.arange(21, dtype=float) ** 2)[:, None, None]
    result = compute.diffusion_blocks(pos, 5.0, 1, fit_fs=(5.0, 10.0))
    assert result["block_frame_counts"] == [5, 4, 4, 4, 4]
    assert result["n_frames_used"] == 21
    assert [(b["frame_start"], b["frame_stop"]) for b in result["blocks"]] == [
        (0, 5), (5, 9), (9, 13), (13, 17), (17, 21)]
    expected = []
    for block in result["blocks"]:
        segment = pos[block["frame_start"]:block["frame_stop"]]
        # Independent direct MSD and two-point slope, without the FFT or fit helper.
        values = [np.mean((segment[lag:] - segment[:-lag]) ** 2) for lag in (1, 2)]
        expected.append((values[1] - values[0]) / 5 / 2 * 0.1)
        assert block["fit_range_fs"] == [5.0, 10.0]
        assert block["fit_points"] == 2
        assert block["duration_fs"] == (len(segment) - 1) * 5
    assert result["d_blocks_cm2_s"] == pytest.approx(expected)
    assert result["d_err_cm2_s"] == pytest.approx(np.std(expected, ddof=1) / np.sqrt(5))
    assert result["estimate_for"] == "mean_of_block_D"


@pytest.mark.parametrize("fit_fs", [None, (10.0, 50.0)])
def test_short_blocks_do_not_rescale_or_fallback_from_main_window(fit_fs):
    pos = (np.arange(21, dtype=float) ** 2)[:, None, None] * np.ones((1, 1, 3))
    result = compute.msd_analysis(pos, None, None, symbols=["H"], dt_fs=5,
                                  remove_drift=False, fit_fs=fit_fs)
    assert result["D_cm2_s"] is not None
    assert result["fit_range_fs"] == [10.0, 50.0]
    error = result["error"]
    assert error["d_err_cm2_s"] is None
    assert error["d_blocks_cm2_s"] == []
    assert error["reason_code"] == "fit_range_not_available_in_all_blocks"
    assert all(b["fit_range_fs"] == [10.0, 50.0] for b in error["blocks"])


def test_block_count_is_reduced_so_every_block_covers_the_same_window():
    pos = np.arange(21, dtype=float)[:, None, None]
    result = compute.diffusion_blocks(pos, 5, 1, fit_fs=(10, 20))
    assert result["n_blocks"] == 4 and result["n_blocks_requested"] == 5
    assert all(b["fit_points"] == 3 and b["fit_range_fs"] == [10.0, 20.0] for b in result["blocks"])
    assert result["fit_range_fs"] == [10.0, 20.0]
    assert result["d_err_cm2_s"] is not None and "4" in result["adjusted"]


@pytest.mark.parametrize("language", ["ja", "en"])
def test_missing_error_reason_is_localized(language, monkeypatch):
    from adit import lang

    monkeypatch.setattr(lang, "LANGUAGE", language)
    result = compute.diffusion_blocks(np.zeros((7, 1, 1)), 5, 1)
    assert result["reason_code"] == "insufficient_blocks"
    assert ("Cannot form" if language == "en" else "ブロックを 2 つ以上") in result["reason"]
