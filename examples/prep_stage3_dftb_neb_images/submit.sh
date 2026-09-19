#!/bin/bash
# ADIT 0.1.0a1 が生成。像の計算を 1 つずつ順に実行する (互いに独立なので、失敗しても次へ進む)
cd "$(dirname "$0")"
failed=""
for d in image_00 image_01 image_02 image_03 image_04 image_05 image_06; do
  echo "== $d"
  bash "$d/submit.sh" || failed="$failed $d"
done
if [ -n "$failed" ]; then echo "failed:$failed" >&2; exit 1; fi
