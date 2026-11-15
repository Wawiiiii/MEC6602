#!/usr/bin/env bash
# Grid-study matrix: every grid 9x9 -> 2049x2049 for each (Mach, alpha) case below, one run after
# the other (each run uses all cores). Same settings as the M 0.5 alpha 0 study: euler_inputs.txt
# is the template, only mesh, mach, alpha, fmg_levels, limiter_k and the output paths change.
# FMG start (4 levels) on 2049 only. Farfield vortex off (default).
# Limiter: off (limiter_k = 1e6) on the smooth subsonic cases (M 0.5), on (default K = 5) on the
# transonic and supersonic cases (shocks).
# Run from the HW2 folder:  nohup setsid ./run_matrix.sh > results/matrix.log 2>&1 &
# The input file of each run is kept next to its outputs (results/<case>/<grid>_input.txt).

set -u
cd "$(dirname "$0")"

EXE=./build/euler
TEMPLATE=euler_inputs.txt
GRIDS="9x9 17x17 33x33 65x65 129x129 257x257 513x513 1025x1025 2049x2049"

# folder name, Mach, alpha, limiter (on/off)
CASES="Mach_05_alpha_0_nolimiter:0.5:0.0:off
Mach_05_alpha_125_nolimiter:0.5:1.25:off
Mach_08_alpha_0:0.8:0.0:on
Mach_08_alpha_125:0.8:1.25:on
Mach_125_alpha_0:1.25:0.0:on
Mach_125_alpha_125:1.25:1.25:on"

for c in $CASES; do
    IFS=: read -r dir mach alpha lim <<< "$c"
    out=results/$dir
    mkdir -p "$out"
    for g in $GRIDS; do
        fmg=1
        [ "$g" = 2049x2049 ] && fmg=4
        inp=$out/${g}_input.txt
        sed -e "s|^mesh *=.*|mesh  = grids/$g.x|" \
            -e "s|^mach *=.*|mach  = $mach|" \
            -e "s|^alpha *=.*|alpha = $alpha|" \
            -e "s|^fmg_levels *=.*|fmg_levels = $fmg|" \
            -e "s|^surface *=.*|surface = $out/s$g.dat|" \
            -e "s|^field *=.*|field   = $out/f$g.dat|" \
            -e "s|^summary *=.*|summary = $out/${g}_summary.txt|" \
            "$TEMPLATE" > "$inp"
        [ "$lim" = off ] && printf '\nlimiter_k = 1e6   # no limiting (smooth subsonic case)\n' >> "$inp"
        echo "=== $(date '+%F %T')  $dir  $g  (M $mach, alpha $alpha, fmg $fmg, limiter $lim)"
        "$EXE" "$inp" > "$out/${g}_terminal.txt" 2>&1
        echo "    exit $?  $(grep -v '^#' "$out/${g}_summary.txt" 2>/dev/null | tail -1)"
    done
done
echo "=== $(date '+%F %T')  matrix done"
