#!/bin/bash
# ADIT 0.1.0a1 が生成。変位の計算を 1 つずつ順に実行する (互いに独立なので、失敗しても次へ進む)
cd "$(dirname "$0")"
failed=""
for d in disp-001 disp-002 disp-003 disp-004 disp-005 disp-006 disp-007 disp-008 disp-009 disp-010 disp-011 disp-012; do
  echo "== $d"
  bash "$d/submit.sh" || failed="$failed $d"
done
if [ -n "$failed" ]; then echo "failed:$failed" >&2; exit 1; fi
