# HW2 – results so far (2026-11-13)

All runs: NACA0012, single-block O-grids of `grids/` (farfield ~150 chords), second-order Roe + MUSCL /
Hemker-Koren (limiter eps = (5h)^3), extrapolated wall pressure (DEVIATIONS.md A20), forces from the
wall-face pressure, chord 1, moment about (0.25, 0), CM nose-up positive. "converged" = density
residual dropped by 1e-11 to 1e-12 (machine zero). Times: one core of an AMD Ryzen 7 9800X3D (serial).

Raw files: `results/runs_bdf1/` (s* = surface Cp + forces in the header, h* = history), multigrid runs
`results/mg_*`, the 513x513 M 0.5 case with field + plots in `results/M05_a125_513/`.
Plots: `python post/plot_results.py --history H --surface S --field F --out DIR`.
Orders: `python post/order_study.py s65.dat s129.dat s257.dat s513.dat s1025.dat --cells 64 128 256 512 1024`.

---

## 1. Grid convergence study, M = 0.5, alpha = 1.25 deg (Part 3, Vassberg & Jameson case)

Implicit BDF1 + JFNK-GMRES (SGS-preconditioned), cfl0 = 1 -> 1e8.

| grid (nodes) | CL | CD | CM | iterations | time |
|---|---|---|---|---|---|
| 65x65 | 0.17580183 | 0.00734443 | -0.00227543 | 18 | 0.3 s |
| 129x129 | 0.17762119 | 0.00166278 | -0.00208267 | 22 | 2.3 s |
| 257x257 | 0.17905556 | 0.00028266 | -0.00221650 | 31 | 17 s |
| 513x513 | 0.17930803 | 0.00005915 | -0.00224442 | 38 | 174 s |
| 1025x1025 | 0.17935083 | 0.00002699 | -0.00224975 | 53 | 1164 s |
| 2049x2049 | (see §6) | | | | |

Observed order p (Newton iteration on G(p), post.pdf slides 26-28) and extrapolated values:

| grids (cells) | p CL | CL_inf | p CD (3 grids) | CD_inf | p CM | CM_inf | p CD with CD_exact = 0 |
|---|---|---|---|---|---|---|---|
| 64-128-256 | 0.34 | (pre-asymptotic) | 2.04 | | n/a (not monotone) | | 2.56 |
| 128-256-512 | **2.51** | 0.179362 | 2.63 | 1.6e-5 | **2.26** | -0.0022518 | 2.26 |
| 256-512-1024 | **2.56** | 0.179360 | **2.80** | 2.2e-5 | **2.39** | -0.0022510 | 1.13 |

- CL, CM: second order (p 2.3-2.6) on the fine grids. CL_inf = 0.17936, CM_inf = -0.002251.
- CD: drops at second order, but tends to a small non-zero limit ~2e-5 (so p with "exact = 0" falls
  to 1.1 on 512->1024). Mostly a lift-related farfield effect (no vortex correction, post.pdf slide 22),
  but alpha = 0 also shows a much smaller floor (~2.4e-6) once the 2049 grid is included:

M = 0.5, alpha = 0, BDF1 + GMRES (raw files `results/Mach_05_alpha_0/`, Linux i9-10850K, 20 threads;
2049 with an FMG start, 4 levels). CL and CM are at machine zero (|CL| < 3e-11, |CM| < 2e-12) on every grid.

| grid (nodes) | CD | iterations | time | p_bar = log2(CD_m / CD_f) |
|---|---|---|---|---|
| 9x9 | 5.1541044652e-02 | 10 | 0.02 s | |
| 17x17 | 3.8128139493e-02 | 12 | 0.02 s | 0.44 |
| 33x33 | 1.9426246278e-02 | 14 | 0.08 s | 0.97 |
| 65x65 | 7.4979461589e-03 | 17 | 0.38 s | 1.37 |
| 129x129 | 1.6390067028e-03 | 20 | 2.0 s | 2.19 |
| 257x257 | 2.6045381790e-04 | 28 | 15 s | 2.65 |
| 513x513 | 3.8539518911e-05 | 36 | 174 s | 2.76 |
| 1025x1025 | 6.3291277464e-06 | 45 | 1749 s | 2.61 |
| 2049x2049 | 2.8252782947e-06 | 21 (+ FMG) | 13106 s | 1.16 |

| triplet (cells) | R | p | CD* |
|---|---|---|---|
| 64-128-256 | 0.2353 | 2.09 | -1.64e-4 |
| 128-256-512 | 0.1610 | 2.64 | -4.04e-6 |
| 256-512-1024 | 0.1452 | 2.78 | 8.60e-7 |
| 512-1024-2048 | 0.1088 | 3.20 | **2.40e-6** |

  CD* from the three finest grids = 2.40e-6 (not 0, the exact value for subsonic inviscid flow).
  p = 3.2 is above the formal order 2: the 2048 grid already feels the floor, so the finest triplet is
  not cleanly asymptotic. Since alpha = 0 carries no circulation, this floor is not the missing
  point vortex; its source is still open (candidates: finite farfield, limiter near the LE/TE,
  wall-pressure extrapolation).

- Comparison with V&J 2010 continuum values (Table 10; C_ref = 1, X_ref = 0.25 as here):

| solver | farfield vortex | CL* | CD* | CM* |
|---|---|---|---|---|
| this code (256-512-1024) | no | 0.179360 | ~2e-5 | -0.0022510 |
| CFL3D | no | 0.179784 | 1.22e-5 | -0.0022706 |
| OVERFLOW | no | 0.179777 | 1.00e-5 | -0.0022626 |
| CFL3D | yes | 0.180352 | -1.3e-7 | -0.0022774 |
| FLO82 | yes | 0.180346 | 5.0e-8 | -0.0022688 |

  The codes without the point vortex all land near CL 0.1798 with a small positive drag floor; with the
  vortex CL rises by ~0.0006 and CD -> 0 (V&J conclusions). This code sits with the no-vortex group,
  which supports the farfield explanation of the CD floor (B11). V&J also find FLO82 pre-asymptotic below
  512x512, like our 64-128-256 triplet. Our CL, CM orders (2.3-2.6) match FLO82's (2.1); CFL3D and
  OVERFLOW are first order for CL and CM.
- Same grids with the slide-14 wall (P1 = P2, "mirror"): CL p = 1.49, CM p ~ 1.0 -> reason for A20.
- 513x513 solution: stagnation Cp = 1.0639 (isentropic value at M 0.5: 1.0638), suction peak
  Cp = -0.721 at x/c ~ 0.06, entropy error |s/s_inf - 1| < 8e-3 (max at the leading edge).
  Plots: `results/M05_a125_513/plots/` (Mach, pressure, Cp contours, surface Cp, entropy, history).

## 2. Methods compared (same converged CL / CD / CM to 1e-9 or better)

| grid | method | iterations / cycles | time |
|---|---|---|---|
| 33x33 | explicit RK5 (slide 18, CFL 1.15) | 11505 | 5 s |
| 33x33 | BDF1 / BDF2 / ROW + GMRES | 15 / 17 / 14 | 0.05-0.10 s |
| 65x65 | explicit RK5 | 34899 | 60 s |
| 65x65 | BDF1 + GMRES | 18 | 0.3 s |
| 65x65 | multigrid W, RK5 smoother | 174 cycles | 1.2 s |
| 65x65 | multigrid W, SGS smoother | 42 cycles | 0.44 s |
| 129x129 | multigrid W, SGS smoother | 44 cycles | 2.3 s |
| 257x257 | BDF1 + GMRES | 31 | 17 s |
| 257x257 | BDF1 + GMRES with FMG start (33 -> 65 -> 129) | 19 | 12.7 s |
| 257x257 | multigrid W, SGS smoother | 48 cycles | 10.3 s |
| 513x513 | BDF1 + GMRES | 38 | 174 s |
| 513x513 | multigrid W, SGS smoother | 173 cycles | 122 s |

Not working yet: multigrid with the RK5 smoother on 257x257 (W-cycle stalls at drop 1.3e-2) and the
V-cycle with RK5 on 33/65 (stalls on a wrong state, CL = -1.0); first-order coarse operators with a
second-order fine operator. W-cycle + SGS smoother works.

## 3. Transonic M = 0.8, alpha = 1.25 deg (BDF1 + GMRES)

| grid | CL | CD | CM | iterations | time |
|---|---|---|---|---|---|
| 65x65 | 0.36379976 | 0.02769958 | -0.04182657 | 27 | 0.5 s |
| 129x129 | 0.36027449 | 0.02380642 | -0.04070550 | 45 | 5.0 s |
| 257x257 | 0.35718619 | 0.02288512 | -0.03969695 | 87 | 45 s |

Multigrid W + SGS: 257x257 converges in 67 cycles (13 s), identical values; on 65x65 and 129x129 it
stalls (drop 4e-2 and 5e-3 after 500 cycles). CL still changes with the grid (0.364 -> 0.360 ->
0.357): finer grids needed; shock quality (Cp oscillations, DEVIATIONS B2) not looked at yet.

## 4. Supersonic M = 1.5, alpha = 0 (BDF1 + GMRES)

| grid | CD | iterations | note |
|---|---|---|---|
| 33x33 | 0.10344 | 22 | |
| 65x65 | 0.10192 (mirror wall) | 22 | |
| 257x257 | 0.09600 | 80 | needs cfl0 = 1 (cfl0 = 10 fails, A21) |

## 5. Tests (all passing)

test_roe 47/47, test_bc 29/29, test_explicit 12/12, test_output 30/30, test_jacobian 19/19,
test_linsolve 38/38, test_implicit 7/7, test_multigrid 37/37 (agent build).

## 5b. OpenMP

All timings above are serial (1 core). With OpenMP (case key `threads`, default all cores) the
results are bit-identical; speed-ups measured on a busy machine with 8 threads: residual 4.3x,
explicit RK 4.3x, implicit BDF1 + GMRES 1.4x (serial SGS preconditioner). See PROGRESS.md.

## 6. Pending (were running when this file was written)

- 2049x2049 M 0.5 alpha 1.25, BDF1 + GMRES (h2049 / s2049 in runs_tmp when finished).
- 1025x1025 M 0.5 alpha 0 (a0_1025).
- Multigrid W + SGS at M 0.5 on 1025x1025.
- An agent was working on the multigrid stalls (DEVIATIONS B18 / B19) when this was pushed.
Note: `grids/2049x2049.x` (82 MB) is not in git: copy it from Moodle into `grids/`.
Note: `results/M05_a125_513/field.dat` (40 MB) is not in git; regenerate it with
`build/Release/euler.exe results/M05_a125_513/case.txt` (~3 min). The plots are in git.
If they are not in `results/` they did not finish before the push; rerun with
`build/Release/euler.exe <case>` (see PROGRESS.md for the settings of the fine grids).
