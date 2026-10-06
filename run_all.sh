#!/bin/sh
# Check the four published codes; exit 0 only if all pass.
set -e
cd "$(dirname "$0")"
for q in 2 3 4 5; do
  python3 check_code.py "$q" "codes/code_q$q.txt"
done
echo "ALL PASS"
