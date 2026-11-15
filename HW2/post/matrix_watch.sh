#!/usr/bin/env bash
# Regenerates results/MATRIX_REPORT.md every 10 minutes while run_matrix.sh is running, once more
# at the end, then exits. Run from the HW2 folder:
#   nohup setsid post/matrix_watch.sh > /dev/null 2>&1 &
cd "$(dirname "$0")/.."
while pgrep -f '^bash ./run_matrix.sh' > /dev/null; do
    python3 post/matrix_report.py
    sleep 600
done
python3 post/matrix_report.py
