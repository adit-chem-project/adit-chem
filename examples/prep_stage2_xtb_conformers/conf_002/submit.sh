#!/bin/bash
# adit 0.1.0a1 が生成 (2026-09-11T19:01:43+00:00)。実行先: ローカル (プロファイル local)
# 使い方: 実行ファイルが PATH にある状態で、このディレクトリで  bash submit.sh
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
adit_dir=$(pwd)
trap 'grep -m1 -E "xtb version" "$adit_dir/output.log" > "$adit_dir/code_version.txt" 2>/dev/null' EXIT
xtb struct.xyz --gfn 2 --chrg 0 --uhf 0 --acc 1 --etemp 300 --input xtb.inp --parallel 1 --opt normal --json > output.log 2>&1
