#!/bin/bash
# adit 0.1.0a1 が生成 (2026-09-12T21:20:59+00:00)。実行先: ローカル (プロファイル local)
# 使い方: 実行ファイルが PATH にある状態で、このディレクトリで  bash submit.sh
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
adit_dir=$(pwd)
trap 'grep -m1 -E "adit-openmm: openmm" "$adit_dir/output.log" > "$adit_dir/code_version.txt" 2>/dev/null' EXIT
python3 run_openmm.py > output.log 2>&1
