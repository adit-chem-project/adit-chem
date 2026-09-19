#!/bin/bash
# ADIT 0.1.0a1 が生成。配座の計算を 1 つずつ順に実行する (互いに独立なので、失敗しても次へ進む)
cd "$(dirname "$0")"
failed=""
for d in conf_001 conf_002 conf_003; do
  echo "== $d"
  bash "$d/submit.sh" || failed="$failed $d"
done
if [ -n "$failed" ]; then echo "failed:$failed" >&2; exit 1; fi
