from pathlib import Path

import numpy as np
import pytest

from adit.analysis.readers_qc import read_gamess, read_gaussian, read_qchem

DATA = Path(__file__).parent / "data"
HARTREE_EV = 27.211386245988


def test_gaussian_frequencies_ir_and_raman():
    out = read_gaussian(DATA / "gaussian_dvb_raman.out.gz")
    assert len(out.frequencies_cm1) == 54           # 3N-6 (N = 20)
    assert out.frequencies_cm1[:3] == [53.1117, 84.6059, 149.2231]
    assert out.ir_intensities[:3] == [0.0340, 0.0000, 0.3732]
    assert out.raman_activities[:3] == [0.0000, 2.6793, 0.0000]
    assert out.energies_ev[-1] == pytest.approx(-382.308266561 * HARTREE_EV, rel=1e-12)


def test_gaussian_structure_and_charges():
    out = read_gaussian(DATA / "gaussian_dvb_raman.out.gz")
    assert len(out.frames) >= 1
    atoms = out.frames[-1]
    assert len(atoms) == 20
    assert atoms.get_chemical_symbols().count("C") == 10
    kinds = {c["definition"]: c for c in out.charges}
    assert len(kinds["Mulliken"]["values"]) == 20
    assert kinds["Mulliken"]["values"][0] == pytest.approx(-0.004506)
    assert len(kinds["Mulliken (水素を重原子にまとめた)"]["values"]) == 10


def test_gamess_drops_the_six_zero_modes_and_converts_the_ir_unit():
    out = read_gamess(DATA / "gamess_dvb_ir.out.gz")
    assert len(out.frequencies_cm1) == 54
    assert out.frequencies_cm1[:3] == [47.87, 81.18, 152.30]
    assert out.ir_intensities[2] == pytest.approx(0.3414204, abs=0.01)
    assert out.energies_ev[-1] == pytest.approx(-382.0506353313 * HARTREE_EV, rel=1e-12)


def test_gamess_reads_the_bohr_geometry_and_mulliken_charges():
    out = read_gamess(DATA / "gamess_dvb_ir.out.gz")
    assert out.frames and len(out.frames[-1]) == 20
    assert np.abs(out.frames[-1].get_positions()).max() < 10.0
    charges = out.charges[-1]["values"]
    assert len(charges) == 20
    assert charges[0] == pytest.approx(-0.004343)


def test_qchem_frequencies_ir_raman_and_charges():
    out = read_qchem(DATA / "qchem_dvb_raman.out.gz")
    assert len(out.frequencies_cm1) == 54
    assert out.frequencies_cm1[0] == pytest.approx(-106.88)
    assert out.ir_intensities[2] == pytest.approx(0.419)
    assert out.raman_activities[0] == pytest.approx(2.048)
    assert len(out.charges) == 1 and len(out.charges[0]["values"]) == 20


def test_orca_raman_through_load_run(tmp_path):
    import gzip

    from adit.analysis.readers import load_run

    (tmp_path / "output.log").write_bytes(gzip.decompress((DATA / "orca_dvb_raman.out.gz").read_bytes()))
    (tmp_path / "orca.inp").write_text("! B3LYP def2-SVP NumFreq\n", encoding="utf-8")
    data = load_run(tmp_path)
    assert data.raman_activities is not None
    assert len(data.raman_activities) == 54
    assert data.raman_activities[1] == pytest.approx(2.788888)


def test_load_run_uses_the_new_readers(tmp_path):
    import gzip

    from adit.analysis.readers import load_run

    (tmp_path / "gaussian.gjf").write_text("#p b3lyp\n", encoding="utf-8")
    (tmp_path / "output.log").write_bytes(gzip.decompress((DATA / "gaussian_dvb_raman.out.gz").read_bytes()))
    data = load_run(tmp_path)
    assert data.code == "gaussian"
    assert data.frequencies_cm1 and len(data.frequencies_cm1) == 54
    assert data.raman_activities and len(data.raman_activities) == 54
    assert data.frames and len(data.frames[-1]) == 20
    assert data.charges and data.charges[0]["definition"].startswith("Mulliken")


def test_a_directory_without_an_output_falls_back(tmp_path):
    from adit.analysis.readers import load_run

    (tmp_path / "qchem.in").write_text("$molecule\n0 1\nH 0 0 0\n$end\n", encoding="utf-8")
    (tmp_path / "final.xyz").write_text("1\n\nH 0.0 0.0 0.0\n", encoding="utf-8")
    data = load_run(tmp_path)
    assert data.code == "qchem"
    assert not data.frequencies_cm1


def test_gamess_marks_imaginary_frequencies_negative(tmp_path):
    out_file = tmp_path / "gamess.out"
    out_file.write_text("       FREQUENCY:       342.20 I     14.71       10.15\n"
                        "    IR INTENSITY:      0.10000     0.00000     0.00000\n", encoding="utf-8")
    out = read_gamess(out_file)
    assert out.frequencies_cm1 == [-342.20, 14.71, 10.15]
    assert len(out.ir_intensities) == 3
