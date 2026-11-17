#!/usr/bin/env bash
# Hourly mechanical status into results/PROGRESS_NIGHT.md, independent of any Claude session:
# per machine, the queue process, the run it is on and its last iteration line. Stops by itself
# after 36 hours. Run from the HW2 folder:  nohup setsid post/night_status.sh > /dev/null 2>&1 &
cd "$(dirname "$0")/.."
for n in $(seq 36); do
    {
        echo
        echo "## $(date '+%Y-%m-%d %H:%M') automatic status (post/night_status.sh)"
        for hq in "j5015-01 j01 matrix.log" "j5015-18 j18 matrix_j18.log" "j5015-17 j17 matrix_j17.log" "j5015-04 j04 matrix_j04.log"; do
            set -- $hq
            if [ "$1" = "$(hostname -s)" ]; then
                q=$(pgrep -f "^python3 run_queue.py" > /dev/null && echo running || echo stopped)
                inp=$(pgrep -af 'build/(kappa/)?euler results/' | awk '{print $NF}' | head -1)
            else
                q=$(timeout 20 ssh -o BatchMode=yes -o ConnectTimeout=6 $1 'pgrep -f "^python3 run_queue.py" > /dev/null && echo running || echo stopped' 2>/dev/null)
                inp=$(timeout 20 ssh -o BatchMode=yes -o ConnectTimeout=6 $1 "pgrep -af 'build/(kappa/)?euler results/'" 2>/dev/null | awk '{print $NF}' | head -1)
            fi
            last=""
            [ -n "$inp" ] && last=$(grep -E '^ *[0-9]+ ' "${inp/_input/_terminal}" 2>/dev/null | tail -1 | awk '{print "it " $1 ", drop " $3 ", CFL " $9 ", lin its " $10}')
            echo "- $1 (queue $2: ${q:-unreachable}): ${inp:-no solver running} ${last:+| $last}"
            tail -n 1 "results/$3" 2>/dev/null | grep -q 'queue done' && echo "  queue $2 finished"
        done
        new=$(cat results/ALERTS*.txt 2>/dev/null | awk -v t="$(date -d '-1 hour' '+%Y-%m-%d %H:%M:%S')" '($1" "$2) >= t')
        [ -n "$new" ] && echo "- alerts in the last hour:" && echo "$new" | sed 's/^/    /'
    } >> results/PROGRESS_NIGHT.md
    python3 post/matrix_report.py 2> /dev/null
    sleep 3600
done
