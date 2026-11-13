# HW2 results

Solver output (history, surface Cp, Tecplot field) and plots of the runs used in the report.
Each sub-folder is one case, e.g. `M05_a125_513/`:

- `history.dat` : iteration, ||R_rho||, drop, CL, CD, CM, seconds
- `surface.dat` : x/c, y/c, Cp, p on the wall faces (CL, CD, CM in the header)
- `field.dat`   : Tecplot ordered zone, node values X Y Rho U V P Mach Cp Entropy
- `plots/`      : PNGs from `post/plot_results.py`

Settings unless stated otherwise: second order (Roe + MUSCL / Hemker-Koren, eps = (5h)^3),
extrapolated wall pressure, chord 1, moment about (0.25, 0), converged to machine zero.
