#!/bin/bash
# adit 0.1.0a1 が生成 (2026-09-11T17:13:47+00:00)。実行先: ローカル (プロファイル local)
# 使い方: 実行ファイルが PATH にある状態で、このディレクトリで  bash submit.sh
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
adit_dir=$(pwd)
trap 'grep -m1 -E "GROMACS version:" "$adit_dir/adit.log" > "$adit_dir/code_version.txt" 2>/dev/null' EXIT
gmx grompp -f grompp.mdp -c conf.gro -p topol.top -o adit.tpr > grompp.log 2>&1 && gmx mdrun -deffnm adit -ntmpi 1 -ntomp 1 > output.log 2>&1
