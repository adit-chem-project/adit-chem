#!/usr/bin/env python3
# ADIT が生成。各歪みの応力を集めて、弾性定数の表 (elastic_constants.csv) を書く。ADIT が入った Python で実行する
from pathlib import Path
from adit.elastic_setup import collect

print(collect(Path(__file__).resolve().parent)['summary'])
