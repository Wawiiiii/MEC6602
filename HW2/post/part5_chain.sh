#!/bin/bash
# Part 5 chain (2026-11-17): 1025 with the fine-grid fixed-CL settings, then 2049 only if 1025 converged.
cd "$(dirname "$0")/.."
python3 run_queue.py part5b
st=$(grep -v '^#' results/Mach_08_CL_04/1025x1025_summary.txt 2>/dev/null | awk 'NF' | tail -1 | awk '{print $9}')
if [ "$st" = "converged" ]; then
    python3 run_queue.py part5c
else
    echo "=== $(date '+%F %T')  1025 status '$st': 2049 NOT started" >> results/matrix_part5b.log
fi
echo "=== $(date '+%F %T')  part5 chain done" >> results/matrix_part5b.log
