#!/bin/bash
# adit 0.1.0a1 が生成 (2026-09-12T21:20:51+00:00)。実行先: ローカル (プロファイル local)
# 使い方: 実行ファイルが PATH にある状態で、このディレクトリで  bash submit.sh
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
adit_dir=$(pwd)
trap 'grep -m1 -E "adit-psi4: psi4" "$adit_dir/output.log" > "$adit_dir/code_version.txt" 2>/dev/null' EXIT
psi4 -i input.dat -o output.log -n 1 > stdout.log 2>&1
