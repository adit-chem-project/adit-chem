import os
from pathlib import Path

import numpy as np
import pytest

from adit.analysis.readers import load_run
from adit.analysis.readers_openmx import read_openmx

SAMPLE = """\
Atoms.SpeciesAndCoordinates.Unit   Ang # Ang|AU
Atoms.UnitVectors.Unit             Ang # Ang|AU
<Atoms.UnitVectors
  10.0   0.0   0.0
   0.0  10.0   0.0
   0.0   0.0  10.0
Atoms.UnitVectors>

  Utot.         -8.033540955197

***********************************************************
            Eigenvalues (Hartree) for SCF KS-eq.
***********************************************************

   HOMO =  4
   Eigenvalues
          1  -0.70000000000000  -0.70000000000000
          2  -0.40000000000000  -0.40000000000000
***********************************************************

***********************************************************
                   Mulliken populations
***********************************************************

  Total spin moment (muB)   0.000000000

                    Up spin      Down spin     Sum           Diff
      1    C      2.509748760  2.509748760   5.019497520   0.000000000
      2    H      0.372562810  0.372562810   0.745125620   0.000000000

 Sum of MulP: up   =     4.00000 down          =     4.00000

  Decomposed Mulliken populations
      1    C      9.999999999  9.999999999   9.999999999   0.000000000

***********************************************************
       Fractional coordinates of the final structure
***********************************************************

     1      C     0.00000000000000   0.00000000000000   0.00000000000000
     2      H     0.10000000000000   0.20000000000000   0.30000000000000
***********************************************************
"""


def test_reads_energy_structure_and_populations(tmp_path):
    p = tmp_path / "openmx.out"
    p.write_text(SAMPLE, encoding="utf-8")
    out = read_openmx(p)
    assert out.energies_ev[0] == pytest.approx(-8.033540955197 * 27.211386245988)
    assert out.homo_index == 4
    assert out.populations == [5.019497520, 0.745125620]
    assert len(out.frames) == 1
    atoms = out.frames[0]
    assert atoms.get_chemical_symbols() == ["C", "H"]
    assert np.allclose(atoms.get_positions()[1], [1.0, 2.0, 3.0])
    assert len(out.eigenvalues_ev) == 4


def test_without_a_cell_the_structure_is_not_guessed(tmp_path):
    p = tmp_path / "openmx.out"
    p.write_text(SAMPLE.replace("<Atoms.UnitVectors", "#<Atoms.UnitVectors"), encoding="utf-8")
    out = read_openmx(p)
    assert out.frames == [] and out.note_no_cell


def test_load_run_uses_it(tmp_path):
    (tmp_path / "openmx.dat").write_text("System.Name x\n", encoding="utf-8")
    (tmp_path / "output.log").write_text(SAMPLE, encoding="utf-8")
    data = load_run(tmp_path)
    assert data.code == "openmx"
    assert data.energies_ev and data.charges
    assert data.charges[0]["values"][0] == pytest.approx(5.01949752)
    assert any("最高被占準位" in n or "highest occupied" in n for n in data.notes)


@pytest.mark.skipif(not os.environ.get("ADIT_OPENMX_OUT"), reason="needs a real .out in ADIT_OPENMX_OUT")
def test_real_output_file():
    out = read_openmx(Path(os.environ["ADIT_OPENMX_OUT"]))
    assert out.energies_ev
    if out.populations:
        total = sum(out.populations)
        assert total > 0 and abs(total - round(total)) < 0.01
