from pathlib import Path

import pytest

from adit.analysis.bader import RUN_NOTE, read_acf, summary_lines

ACF = Path(__file__).parent / "data" / "bader_ACF.dat"


def test_reads_the_known_charges():
    bc = read_acf(ACF)
    assert bc.electrons.size == 2
    assert bc.electrons[0] == pytest.approx(6.0, abs=1e-3)
    assert bc.electrons[1] == pytest.approx(2.0, abs=2e-3)
    assert bc.total_electrons == pytest.approx(8.0, abs=3e-3)


def test_net_charge_uses_user_supplied_valence():
    bc = read_acf(ACF)
    net = bc.net_charge([7.0, 1.0])
    assert net[0] == pytest.approx(1.0, abs=1e-3)
    assert net[1] == pytest.approx(-1.0, abs=2e-3)


def test_wrong_number_of_valence_electrons_is_an_error():
    with pytest.raises(ValueError, match="価電子数の個数"):
        read_acf(ACF).net_charge([4.0])


def test_summary_and_note():
    lines = summary_lines(read_acf(ACF), symbols=["O", "H"], valence_electrons=[7.0, 1.0])
    assert any("O 1" in line for line in lines)
    assert any("価電子との差" in line for line in lines)
    assert any("bader" in line for line in RUN_NOTE)


def test_summary_without_valence_says_so():
    assert any("価電子数を渡すと" in line for line in summary_lines(read_acf(ACF)))


def test_rejects_a_file_without_the_table(tmp_path):
    bad = tmp_path / "ACF.dat"
    bad.write_text("nothing here\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Bader の表"):
        read_acf(bad)
