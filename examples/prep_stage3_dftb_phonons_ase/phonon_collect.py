#!/usr/bin/env python3
# ADIT が生成。各変位の力を集めて、フォノン分散 (band.yaml) と (メッシュを与えたなら) 状態密度 (total_dos.dat) を書く。ADIT が入った Python で実行する
from pathlib import Path
from adit.phonon_setup import collect

print(collect(Path(__file__).resolve().parent)['summary'])
