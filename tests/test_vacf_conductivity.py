
import numpy as np
import pytest

from adit.analysis.conductivity import ConductivityError, nernst_einstein, summary_line
from adit.analysis.vacf import CM1_PER_FS1, VacfError, notes, vacf, velocities_from_positions


def test_spectrum_finds_a_known_frequency():
    dt, frames = 2.0, 2000
    nu_cm1 = 1000.0
    omega = 2 * np.pi * nu_cm1 / CM1_PER_FS1      # [1/fs]
    t = np.arange(frames) * dt
    v = np.zeros((frames, 4, 3))
    rng = np.random.default_rng(0)
    phases = rng.uniform(0, 2 * np.pi, size=(4, 3))
    for a in range(4):
        for c in range(3):
            v[:, a, c] = np.cos(omega * t + phases[a, c])
    res = vacf(v, dt)
    peak = res.freq_cm1[int(np.argmax(res.spectrum[1:])) + 1]
    assert peak == pytest.approx(nu_cm1, rel=0.03)
    assert res.vacf[0] == pytest.approx(1.0)
    assert any("cm⁻¹" in n or "cm^-1" in n for n in notes(res))


def test_green_kubo_d_agrees_with_the_msd_slope():
    from adit.analysis.compute import diffusion_fit, fit_range_fs, msd_fft

    rng = np.random.default_rng(3)
    dt, frames, natoms, d_true = 10.0, 1500, 8, 2.0e-4        # Å²/fs
    steps = rng.normal(0.0, np.sqrt(2 * d_true * dt), size=(frames - 1, natoms, 3))
    pos = np.concatenate([np.zeros((1, natoms, 3)), np.cumsum(steps, axis=0)])
    t = np.arange(frames, dtype=float) * dt
    d_msd = diffusion_fit(t, msd_fft(pos), fit_range_fs(t, None), 3)
    v = velocities_from_positions(pos, dt)
    res = vacf(v, dt, source="finite_difference")
    assert res.d_cm2_s == pytest.approx(d_msd, rel=0.35), (res.d_cm2_s, d_msd)
    assert res.source == "finite_difference"


def test_velocity_source_is_recorded():
    v = np.zeros((10, 2, 3)); v[:, :, 0] = 1.0
    assert vacf(v, 1.0).source == "velocities"
    with pytest.raises(VacfError, match="フレーム|frames"):
        vacf(np.zeros((2, 2, 3)), 1.0)
    with pytest.raises(VacfError, match="間隔|time between"):
        vacf(np.zeros((10, 2, 3)), 0.0)


def test_finite_difference_needs_three_frames():
    with pytest.raises(VacfError, match="3 つ以上|three frames"):
        velocities_from_positions(np.zeros((2, 3, 3)), 1.0)


def test_nernst_einstein_matches_a_hand_calculation():
    c = nernst_einstein(d_cm2_s=1.0e-6, n_ions=18, charge=1.0, volume_ang3=855.4, temperature_k=600.0)
    n = 18 / (855.4e-24)
    expected = n * (1.602176634e-19) ** 2 * 1.0e-6 / (1.380649e-23 * 600.0)
    assert c.sigma_s_per_cm == pytest.approx(expected, rel=1e-12)
    assert c.sigma_s_per_cm == pytest.approx(0.0652, rel=0.01)
    assert "Haven" in c.as_dict()["note"] or "相関" in c.as_dict()["note"]
    assert "mS/cm" in summary_line(c)


def test_conductivity_scales_with_charge_squared_and_density():
    base = nernst_einstein(1e-6, 10, 1.0, 1000.0, 300.0)
    doubled_charge = nernst_einstein(1e-6, 10, 2.0, 1000.0, 300.0)
    doubled_count = nernst_einstein(1e-6, 20, 1.0, 1000.0, 300.0)
    assert doubled_charge.sigma_s_per_cm == pytest.approx(4 * base.sigma_s_per_cm)
    assert doubled_count.sigma_s_per_cm == pytest.approx(2 * base.sigma_s_per_cm)


def test_conductivity_refuses_incomplete_input():
    with pytest.raises(ConductivityError, match="電荷|charge"):
        nernst_einstein(1e-6, 10, 0.0, 1000.0, 300.0)
    with pytest.raises(ConductivityError, match="体積|volume"):
        nernst_einstein(1e-6, 10, 1.0, 0.0, 300.0)
    with pytest.raises(ConductivityError, match="拡散係数|diffusion"):
        nernst_einstein(-1.0, 10, 1.0, 1000.0, 300.0)
