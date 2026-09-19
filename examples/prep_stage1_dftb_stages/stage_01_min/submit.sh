#!/bin/bash
# adit 0.1.0a1 が生成 (2026-09-11T17:59:19+00:00)。実行先: ローカル (プロファイル local)
# 使い方: 実行ファイルが PATH にある状態で、このディレクトリで  bash submit.sh
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
adit_dir=$(pwd)
trap 'grep -m1 -E "DFTB\+ .*(release|version)" "$adit_dir/output.log" > "$adit_dir/code_version.txt" 2>/dev/null' EXIT
dftb+ > output.log 2>&1
