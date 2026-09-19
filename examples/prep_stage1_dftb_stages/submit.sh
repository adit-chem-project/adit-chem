#!/bin/bash
# ADIT 0.1.0a1 が生成。段階を順に実行する (前の段階が失敗したらそこで止まる)
cd "$(dirname "$0")"
set -e
for d in stage_01_min stage_02_nvt stage_03_nve; do
  echo "== $d"
  bash "$d/submit.sh"
done
