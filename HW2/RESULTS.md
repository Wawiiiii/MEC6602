# HW2 – results and what to use in the report (updated 2026-11-13, night)

NACA0012, single-block O-grids of `grids/` (farfield ~150 chords). Second-order Roe + MUSCL /
Hemker-Koren (limiter eps = (5h)^3), extrapolated wall pressure (DEVIATIONS A20), forces from the
wall-face pressure, chord 1, moment about (0.25, 0), CM nose-up positive. Run outputs are in
`results/` on the Windows machine only (not in git); every number needed is in this file.

---

## 0. What to use in the report (decision guide)

**One main method for all the analysis (Parts 3, 4, 5):**

| choice | use | why |
|---|---|---|
| spatial scheme | Roe + MUSCL (Hemker-Koren), 2nd order | the homework scheme; order ~2 shown in §2 |
| solver | **implicit BDF1 + JFNK-GMRES**, SGS-preconditioned (`scheme = bdf1`) | fastest and most robust on every case tried (§3) |
| start | **full multigrid start** `fmg_levels = 4` | 257x257 converges in 16-28 iterations instead of 31-87 |
| farfield (subsonic, lifting) | **`farfield_vortex = on`** for the Part 3 grid study | removes the CD floor; CD -> 0 at order ~2.7 (§2). Slide 22 says "not to be done", so present it as an addition and show the off result once |
| wall | `wall_pressure = extrapolated` (default) | slide-14 mirror gives CM order ~1 (A20) |

Settings per grid (case-file keys): `gmres_restart 30 / gmres_max_iter 60 / tol 1e-12` up to 257;
`50 / 150 / 1e-11` on 513; `50 / 200 / 1e-11` on 1025; `30 / 300 / 1e-10` on 2049.
Everything else at the defaults of `cases/example.txt`.

**Mention in one line each (verification, not results):**
- BDF2 and ROW: implemented and selectable; give the same CL/CD/CM to all 8 printed digits.
  Not worth it for steady flow: BDF2 = same time as BDF1, ROW ~2x slower (two linear solves per step).
- Multigrid (Part 2 acceleration): FAS W-cycle, SGS smoother (`scheme = mg`), converges on every
  case, same forces. Show ONE comparison table (§3); it only beats BDF1 on 513 subsonic.
- Explicit RK5: works but ~10^4 iterations (65x65: 34899 it / 60 s vs 18 it / 0.3 s implicit).
  Good to show the speed-up factor once.

**Part 3 (grid study):** M 0.5, alpha 1.25 (Vassberg & Jameson case), BDF1, vortex correction on,
grids 129-257-513-1025 (+2049 if it finishes), orders of CL, CD (exact 0), CM by Newton iteration
(post.pdf slides 26-28). Show alpha = 0 CD as the "no lift" check. Show the vortex-off CD floor in
one sentence / one plot (it explains why the correction matters).

**Part 4 (flow analysis):** M 0.5, M 0.8 and M 1.5 at alpha 1.25 on 513x513 (BDF1). Features to show
are in §4 (stagnation Cp vs exact, suction peak, shock position, bow-shock standoff, post-shock
entropy vs Rankine-Hugoniot). M 1.25 also works; M 1.5 has the cleaner (closer) bow shock.

**Part 5 (fixed CL, M 0.8, CL 0.4):** not done yet.

**Open:** reference chord 1 vs the grid's 1.00893 (sharp-TE extension) – check what Vassberg &
Jameson normalise with (the PDF could not be downloaded here; put it in `docs/`).

---

## 1. Run status (2026-11-13)

- BDF1 / BDF2 / ROW + GMRES with FMG start, 8 flow cases (M 0.5, 0.8, 1.25, 1.5 x alpha 0, 1.25):
  grids 9-257 all converged (machine zero); 513: 12/24 converged, then **4 stalled** (GMRES
  stagnation at high CFL, fixed by CFL backtracking, DEVIATIONS A24). Rerun with the fix:
  **M 0.8 alpha 1.25 BDF1 converged in 33 iterations (4 rejected steps), BDF2 in 36 (5)**, both to
  drop 1e-11, CL 0.354757 = the ROW value; then 1025 and four 2049 checks (BDF1, alpha 1.25,
  each Mach). **Open (B20): M 1.25 alpha 1.25 BDF1 on 513 oscillates** (residual jumps when the
  CFL passes ~100; not a GMRES failure) – being diagnosed; it creeps down (drop 4.8e-4 at
  iteration 42). **M 1.25 alpha 0 BDF2 on 513: converged in 42 iterations with the backtracking**
  (stalled before). M 1.5 alpha 0 BDF1 on 513: 18 iterations.
- Multigrid W + SGS (CFL 10), 8 flow cases, grids 65-513: 29/32 converged (65-257 all, 513 at
  M 0.5, M 0.8, M 1.5 alpha 0 in 349 cycles). **513 M 1.25 (both alpha): stalls at drop ~2e-3**
  (B21, being diagnosed; forces converged to 5 digits). 513 M 1.5 alpha 1.25: slow, drop 5e-8
  after 800 cycles, still dropping. => multigrid is fine as the Part 2 comparison on subsonic /
  transonic cases; do not use it for supersonic fine grids.
- 513 M 1.25 alpha 1.25 BDF1 (B20) did converge in the end: 168 iterations (BDF2 156) vs 17 on 257.
- Tests on the final code (CFL backtracking, vortex option, multigrid fixes): all 10 programs pass
  (`python tests/run_all_tests.py build/fix/Release`).
- Every converged BDF2 / ROW / multigrid run: same CL/CD/CM as BDF1 to all printed digits.
- 9x9: the FMG start crashed (grid cannot be coarsened) – fixed, start skipped on 9x9.

## 2. Grid convergence, M = 0.5, alpha = 1.25 (Part 3)

BDF1 + GMRES, FMG start.

| grid | CL (vortex on) | CD (vortex on) | CM (vortex on) | CD (vortex off) | CD alpha = 0 |
|---|---|---|---|---|---|
| 129 | 0.17856257 | 0.00164336 | -0.00209337 | 0.00166278 | 0.00163901 |
| 257 | 0.18001125 | 0.00026200 | -0.00222798 | 0.00028266 | 0.00026045 |
| 513 | 0.18026629 | 0.00003841 | -0.00225614 | 0.00005915 | 0.00003854 |
| 1025 | (running) | | | 0.00002699 | 0.00000633 |

Observed orders (vortex on, 129-257-513): **CL 2.51 (CL_inf ~ 0.18032), CM 2.26, CD (exact 0) 2.65 /
2.77** (129->257, 257->513). Without the correction CD tends to a floor ~2e-5 (CD order with exact 0
falls to 1.1 on 512->1024): the lift-induced error of a 150-chord farfield without the vortex,
~ CL^2 / (4 pi R) = 1.7e-5. With the correction CD at alpha 1.25 equals the alpha = 0 CD (no lift
error left). The correction raises CL by 0.5 %.

Vortex off (for reference, earlier study): CL p 2.51 / 2.56 (CL_inf 0.17936), CM p 2.26 / 2.39.
Mirror wall (slide 14): CL p 1.49, CM p ~1.0 (reason for the extrapolated wall, A20).

### 2a. alpha = 0 drag, 9x9 to 2049x2049 (vortex off)

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

Comparison with the V&J 2010 continuum values (Table 10; C_ref = 1, X_ref = 0.25 as here), vortex-off alpha 1.25 study:

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

### 2b. Limiter off on the subsonic cases, and the full matrix (started 2026-11-14)

**Why the alpha = 0 CD order of 3.2 is not an order.** The ratio of successive CD changes
(CD_coarse - CD_fine) keeps growing with the limiter on: 4.25, 6.21, 6.89, 9.19 (128 -> 2048). In
the asymptotic range it would settle at a constant (4 for order 2), so no triplet is asymptotic: two
behaviours are mixed, drag falling faster than order 2 on the medium grids (p 2.6-2.8) and a part
that stops decreasing below ~1e-6 (with p 2.78 the 2049 CD would be 1.6e-6, it is 2.8e-6; p_bar
1024 -> 2048 against 0 is 1.16). A small 2049 step makes the fit give a large p, so CD* ~ CD_2049.
CD* is also very sensitive: +-5e-7 on CD_2049 gives p 3.0-3.4 and CD* 1.8e-6 - 3.0e-6. Report
"CD stalls at O(1e-6)", not a three-digit CD*.

**Decision (2026-11-14): limiter off (`limiter_k = 1e6`) on the smooth subsonic cases (M 0.5),
limiter on (K = 5) on the transonic and supersonic cases (shocks).** See DEVIATIONS A26.

Limiter off vs on, M 0.5 alpha 0, CD (exact 0):

| grid | limiter on (K 5) | limiter off | p_bar off |
|---|---|---|---|
| 65x65 | 7.4979e-03 | 2.0062e-03 | 2.27 |
| 129x129 | 1.6390e-03 | 3.2226e-04 | 2.64 |
| 257x257 | 2.6045e-04 | 4.2852e-05 | 2.91 |
| 513x513 | 3.8540e-05 | 4.8055e-06 | 3.16 |
| 1025x1025 | 6.3291e-06 | (running) | |
| 2049x2049 | 2.8253e-06 | (to run) | |

Without the limiter the drag is 4-8x smaller on the same grid: most of the spurious entropy of the
limited runs came from the limiter clipping at the leading-edge suction peak (where the entropy
error peaks). 1025 and 2049 will say whether the ~1e-6 floor also goes away (limiter) or stays
(grid / farfield / TE, B9, B11).

**Fine grids, limiter off (2026-11-14 evening), M 0.5 alpha 0, CD:** 513: 4.8055e-06,
1025: 5.6185e-07 (p_bar 3.10; triplet 256-512-1024 p 3.16, CD* 2.9e-8), **2049: 9.0807e-07 (higher
than 1025; converged to 3.3e-13 in 14 iterations after FMG, so not an iterative error)**. So the
limiter is not the whole story: without it CD drops at order ~3 down to 1025, then a non-vanishing
part of ~1e-6 shows up on 2049 (with the limiter it was ~2.8e-6, the other solver below 2.7e-6).
A rise from 1025 to 2049 means at least two error terms of opposite sign (e.g. a positive part that
does not go to 0 plus a negative O(h^2) part), so no three-grid fit is valid on 512-1024-2048.
Remaining candidates for the ~1e-6 part: finite farfield (~150 c; at alpha 0 the far field is the
thickness doublet, an O(1e-6) drag error is plausible), trailing edge (B9). Next test: cut the
1025 grid at a smaller outer radius and see whether CD changes like 1/R^2.

**kappa study (DEVIATIONS A27), M 0.5 alpha 0, limiter off, unlimited kappa = 0:** CD 513:
1.4525e-05, 1025: 2.9421e-06; p_bar 256->512 2.53, 512->1024 **2.30** (going down towards 2),
vs 3.16 / 3.10 for kappa = 1/3, and 3-5x more drag. Supports the explanation: the drag (spurious
entropy) follows the dissipative error, O(h^3) for kappa = 1/3, O(h^2) for kappa = 0.

**M 0.5 alpha 1.25, limiter off:** CL 0.17935244 / 0.17937717 / 0.17938028, CD 2.6205e-05 /
2.1401e-05 / 2.1694e-05, CM -2.2515423e-03 / -2.2552086e-03 / -2.2556132e-03 (513 / 1025 / 2049).
Triplet 512-1024-2048: CL p 2.99, CL* 0.179381; CM p 3.18, CM* -0.0022557; CD flat at ~2.1e-5 =
the lift-induced farfield floor without the vortex (B11, CL^2/(4 pi R) = 1.7e-5).

Reference: another solver's multigrid results on the same case (values pasted 2026-11-14, source
not recorded here), CD:

| cells | 8 | 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|---|---|---|---|---|
| CD | 2.9240e-02 | 1.6886e-02 | 6.7307e-03 | 2.1339e-03 | 5.4662e-04 | 1.2945e-04 | 3.0458e-05 | 7.4713e-06 | 2.7136e-06 |
| p_bar | | 0.79 | 1.33 | 1.66 | 1.96 | 2.08 | 2.09 | 2.03 | **1.46** |

Its triplets: p 2.08 / 2.11 / 2.27 and CD* -3.4e-7 / 5.2e-7 / 1.47e-6 (128-256-512, 256-512-1024,
512-1024-2048). Clean order 2 up to 1024, then the same flattening on 2048 as ours, at nearly the
same 2049 value (2.71e-6 vs our 2.83e-6). Two different solvers stalling at ~2-3e-6 on 2049 points to
the grids / farfield rather than to our code.

**Queue from 2026-11-15 00:20 (user's order), limiter ON with both slide 16 limiters**
(Hemker-Koren kappa^ 1/3 = default, folders `Mach_<M>_alpha_<a>`; Van Albada kappa^ 0, folders
`..._vanalbada`, DEVIATIONS A28; user: Van Albada is the limiter to use): M 0.8 alpha 0 HK 2049
(running), M 0.5 alpha 0 VA, M 0.5 alpha 1.25 HK, VA, M 0.8 alpha 0 VA, M 0.8 alpha 1.25 HK, VA,
M 1.25 alpha 0 HK, VA, M 1.25 alpha 1.25 HK, VA. M 0.5 alpha 0 HK = `Mach_05_alpha_0` (done).

**Lab machines (2026-11-15 00:46, user's permission for the night, one case per machine, nice 19;
same shared disk, so every result lands in this `results/`; same binaries, so the numbers do not
depend on the machine, only the times do):** j5015-18 (Xeon w7-2595X) = `Mach_08_alpha_125` (HK),
j5015-17 (i7-14700K) = `Mach_08_alpha_125_vanalbada`, j5015-04 (i9-14900K, 31 GB free) =
`Mach_125_alpha_0` (HK) 9 -> 1025 only. Everything else on j5015-01 (i9-10850K, queue j01): M 0.8
alpha 0 HK 2049, M 0.5 alpha 0 VA, M 0.5 alpha 1.25 HK and VA, M 0.8 alpha 0 VA, M 1.25 alpha 0 VA,
M 1.25 alpha 0 HK 2049. Logs: `results/matrix.log` (j01), `results/matrix_j18.log` etc.; alerts
`results/ALERTS.txt`, `results/ALERTS_j18.txt` etc.

**Settings from 2026-11-15 00:59 (user):** every run starts with FMG (4 levels) from 257 up and GMRES restart 200 / 800 iterations from 1025 up; retries = GMRES 200/800 on every grid, then + cfl_max 1e3 (the low CFL caps 50 / 10 are no longer used). Queues restarted with it; running solvers were taken over by the new queue processes.

**Day split, 2026-11-15 10:05 (user: tol 1e-10 on 2049, all lab machines, add M 1.25 alpha 1.25, done by Mon 09:00):** queues j01d, j18d + j18e, j17d + j17e (two per machine: a 2049 run is mostly serial), j04d (grids up to 1025); contents in `run_queue.py` (QUEUES), schedule in `results/PROGRESS_NIGHT.md`. M 1.25 from 1025 up starts with cfl0 10, gmres_fail 1e-2 and the Cauchy stop (20 it, 2e-7).
Later the same day (to keep Mon 09:00): M 1.25 on 1025 / 2049 changed to Cauchy 10 it / 1e-6 with max_iter 150 / 60 (12:30); the second case of each lab queue was moved to its own parallel queue (j18e/f/g, j17e/f at 12:35 and 14:20); at 15:30 the cases still waiting behind slow runs went to new queues j18h (M 0.5 a1.25 HK 2049), j18i (M 1.25 a1.25 HK 9-1025), j17h (M 1.25 a1.25 VA 9-1025), j17i (M 0.5 a1.25 HK 9-1025).

**Matrix (running):** every grid 9 -> 2049, BDF1 + GMRES, `euler_inputs.txt` settings, FMG 4
levels on 2049 only, farfield vortex off. Started with `run_matrix.sh` (M 0.5 alpha 0 limiter off,
9 -> 1025), continued from 2026-11-14 ~15:30 by **`run_queue.py`** (skips grids already converged,
detects stalls, retries), in this order:
`Mach_05_alpha_0_nolimiter_kappa0` (kappa study, 9 -> 1025, DEVIATIONS A27), then 2049 of
`Mach_05_alpha_0_nolimiter`, `Mach_05_alpha_125_nolimiter`, `Mach_08_alpha_0`, `Mach_08_alpha_125`,
`Mach_125_alpha_0`, `Mach_125_alpha_125` (folders in `results/`; each grid keeps its `_input.txt`
and `_terminal.txt`). The limiter-on M 0.5 alpha 1.25 runs done before the switch (9-257) are in
`results/Mach_05_alpha_125_limiter/`.

Stall handling (`run_queue.py`): a run is a stall if its best residual drop has not improved 10x
in the last 50 iterations; stalled, max_iter or diverged runs are retried with attempt 2
(same settings + FMG start from 257 up), attempt 3 (FMG + cfl_max 50, max_iter 1000, window 80),
attempt 4 (FMG + cfl_max 10, max_iter 3000, window 300); every attempt keeps its
`_input_aN.txt` / `_terminal_aN.txt`, and the events are listed in `results/ALERTS.txt`.
First stall (2026-11-14 22:12): M 0.8 alpha 0 513, cold start: the CFL reached ~290 at it 36, the
residual jumped 10x and oscillated at drop ~1e-2 (CL ~1e-4 at alpha 0, the shock moving) -> stall
at it 93. A CFL cap of 50 with FMG was tried first and stuck at drop ~5e-4 (no Newton
convergence at CFL 50; kept as `513x513_*_old_cfl50_killed.txt`), so the retry order was changed
to FMG first: **FMG + normal ramp converged in 17 iterations (drop 9.5e-13), CD 8.3664502e-3 =
the previous campaign's 0.00836645.** **Live tables and orders of every case:
`results/MATRIX_REPORT.md`** (rewritten after every run and every 10 min; by hand:
`python3 post/matrix_report.py`), with a "retry" column for runs that needed attempt 2 or 3.
Copy the final numbers here when it is done.

## 3. Speed (Part 2) – one table is enough

Total over the 8 flow cases, FMG start included (loaded machine, mixed thread counts: ratios to
~30 %, rerun on an idle machine for the report):

| grid | BDF1 | BDF2 | ROW | multigrid W + SGS |
|---|---|---|---|---|
| 65 | 2.5 s (15 it) | 2.6 s (17) | 4.8 s (14) | 13.5 s (65 cycles) |
| 129 | 13.6 s (16) | 14.8 s (19) | 27.3 s (16) | 57.5 s (83) |
| 257 | 109 s (20) | 110 s (22) | 212 s (19) | 272 s (159) |
| 513 (3 subsonic cases) | 453 s (20) | 457 s (23) | 1040 s (21) | 394 s (122) |

Per flow case on 257x257 (iterations or cycles / time, FMG start included) – the ranking holds for
every case: BDF1 ~ BDF2 fastest, ROW ~2x, multigrid 1.3x (subsonic) to 4.5x (supersonic) slower:

| case (257x257) | BDF1 | BDF2 | ROW | multigrid W + SGS |
|---|---|---|---|---|
| M 0.5, a 0 | 20 / 15.4 s | 22 / 14.5 s | 20 / 31.0 s | 104 / 22.8 s |
| M 0.5, a 1.25 | 21 / 17.1 s | 23 / 15.2 s | 20 / 30.8 s | 117 / 23.5 s |
| M 0.8, a 0 | 22 / 15.4 s | 26 / 16.4 s | 22 / 31.9 s | 97 / 20.1 s |
| M 0.8, a 1.25 | 28 / 22.6 s | 33 / 24.2 s | 26 / 42.9 s | 193 / 40.5 s |
| M 1.25, a 0 | 16 / 10.1 s | 19 / 11.1 s | 16 / 21.7 s | 214 / 47.3 s |
| M 1.25, a 1.25 | 17 / 11.4 s | 18 / 11.0 s | 16 / 20.8 s | 214 / 46.6 s |
| M 1.5, a 0 | 16 / 8.3 s | 18 / 9.5 s | 15 / 16.0 s | 175 / 37.0 s |
| M 1.5, a 1.25 | 16 / 8.9 s | 17 / 8.0 s | 15 / 16.9 s | 155 / 33.6 s |

Without FMG start (serial, M 0.5 alpha 1.25): BDF1 18 / 22 / 31 / 38 / 53 iterations on 65-1025.
Explicit RK5: 33x33 11505 it / 5 s, 65x65 34899 it / 60 s.

## 4. Flow features (Part 4), BDF1

| case | grid | CL | CD | CM | features |
|---|---|---|---|---|---|
| M 0.5, a 1.25 | 513 | 0.17931 | 0.000059 | -0.002244 | stagnation Cp 1.0639 (isentropic 1.0641); suction peak -0.721 at x/c 0.066; max Mach 0.671; entropy error < 9e-4 |
| M 0.5, a 0 | 513 | 0 | 0.000039 | 0 | suction -0.485 at x/c 0.12 |
| M 0.8, a 1.25 | 257 | 0.35719 | 0.02289 | -0.03970 | upper shock x/c ~ 0.64 (Cp jump 0.77), weak lower shock x/c ~ 0.33; max Mach 1.38; 513: CL 0.35476 |
| M 0.8, a 0 | 257 | 0 | 0.00849 | 0 | shocks at x/c ~ 0.51 both sides; max Mach 1.25 |
| M 1.25, a 0 | 257 | 0 | 0.09620 | 0 | bow shock standoff 0.32 c; post-shock entropy 0.0050 (normal shock 0.0052); stagnation Cp 1.4288 (pitot 1.4233) |
| M 1.5, a 0 | 257 | 0 | 0.09599 | 0 | standoff 0.089 c; post-shock entropy 0.0300 (0.0295); stagnation Cp 1.5396 (pitot 1.5322) |
| M 1.5, a 1.25 | 257 | 0.06844 | 0.09739 | -0.01415 | |

Entropy on the wall node row of the Tecplot file shows -0.14 at the TE in supersonic cases: an
output artefact (node = average with the wall ghost cells), the interior stays within -6e-4.
Plots: `python post/plot_results.py --history H --surface S --field F --out DIR`.

## 5. Tests

test_roe 47, test_bc 53, test_explicit 12, test_output 30, test_jacobian 19, test_linsolve 38,
test_implicit 7, test_multigrid 41: all passing in the multigrid agent's build (before the CFL
backtracking and vortex changes; rerun pending). One command: `python tests/run_all_tests.py [build/Release]`.

## 6. Fixes of this night (details in DEVIATIONS)

- B18 multigrid RK smoother: period-2 instability at the TE, fixed by correction relaxation 0.5.
- B19 multigrid SGS smoother (and implicit `linear_solver = sgs`): unstable above CFL ~20, defaults now 10.
- A24 GMRES stagnation at high CFL on 513 (M 0.8, M 1.25): CFL backtracking.
- FMG start on 9x9: skipped (cannot coarsen) instead of crashing.
- Farfield vortex correction (slide 22) added as option `farfield_vortex = on`.

## 7. All forces, every flow case and grid (BDF1, vortex correction off)

"same forces" = the other methods that converged to the same CL/CD/CM to all 8 digits (mg =
multigrid W + SGS, run on 65-513). "-" = not finished yet (513 alpha 1.25 at M 0.8 / 1.25 are the
stalled runs being redone with the CFL backtracking; 1025 and 2049 to come). Regenerate with the
campaign summaries: `tables.py` (session scratch script) writes `results/forces_tables.md`.
Coarse grids (9-33) are far from the asymptotic range (CD even rises 9 -> 17 at M 1.5); use 129+
for the study.

### M = 0.5

| grid | CL (a 0) | CD (a 0) | CL (a 1.25) | CD (a 1.25) | CM (a 1.25) | BDF1 it / time a 1.25 | same forces |
|---|---|---|---|---|---|---|---|
| 9 | 0.00000000 | 0.05154104 | 0.18181681 | 0.05498436 | -0.01283231 | 11 / 0.0 s | bdf2, row |
| 17 | 0.00000000 | 0.03812814 | 0.16964818 | 0.03965423 | -0.00542223 | 12 / 0.0 s | bdf2, row |
| 33 | 0.00000000 | 0.01942625 | 0.16805300 | 0.02094337 | -0.00337958 | 12 / 0.1 s | bdf2, row |
| 65 | 0.00000000 | 0.00749795 | 0.17580183 | 0.00734443 | -0.00227543 | 14 / 0.3 s | bdf2, mg, row |
| 129 | 0.00000000 | 0.00163901 | 0.17762119 | 0.00166278 | -0.00208267 | 15 / 2.0 s | bdf2, mg, row |
| 257 | 0.00000000 | 0.00026045 | 0.17905556 | 0.00028266 | -0.00221650 | 21 / 17.1 s | bdf2, mg, row |
| 513 | 0.00000000 | 0.00003854 | 0.17930803 | 0.00005915 | -0.00224442 | 19 / 155.2 s | bdf2, mg, row |
| 1025 | 0.00000000 | 0.00000633 | 0.17935083 | 0.00002699 | -0.00224975 | 53 / 1164 s (no FMG, serial) | - |

### M = 0.8

| grid | CL (a 0) | CD (a 0) | CL (a 1.25) | CD (a 1.25) | CM (a 1.25) | BDF1 it / time a 1.25 | same forces |
|---|---|---|---|---|---|---|---|
| 9 | 0.00000000 | 0.06582499 | 0.24424681 | 0.07215817 | -0.04122004 | 13 / 0.0 s | bdf2, row |
| 17 | 0.00000000 | 0.04509130 | 0.33255842 | 0.05769034 | -0.03100277 | 16 / 0.0 s | bdf2, row |
| 33 | 0.00000000 | 0.02152716 | 0.34890871 | 0.03646500 | -0.03788656 | 15 / 0.2 s | bdf2, row |
| 65 | 0.00000000 | 0.01291188 | 0.36379976 | 0.02769958 | -0.04182657 | 17 / 0.5 s | bdf2, mg, row |
| 129 | 0.00000000 | 0.00925552 | 0.36027449 | 0.02380642 | -0.04070550 | 19 / 2.6 s | bdf2, mg, row |
| 257 | 0.00000000 | 0.00848843 | 0.35718619 | 0.02288512 | -0.03969695 | 28 / 22.6 s | bdf2, mg, row |
| 513 | 0.00000000 | 0.00836645 | - (ROW: 0.354756) | - | - | - | - |

### M = 1.25

| grid | CL (a 0) | CD (a 0) | CL (a 1.25) | CD (a 1.25) | CM (a 1.25) | BDF1 it / time a 1.25 | same forces |
|---|---|---|---|---|---|---|---|
| 9 | 0.00000000 | 0.13056260 | 0.08494550 | 0.13240216 | -0.01978328 | 17 / 0.0 s | bdf2, row |
| 17 | 0.00000000 | 0.11668287 | 0.08337106 | 0.11852599 | -0.01856079 | 17 / 0.0 s | bdf2, row |
| 33 | 0.00000000 | 0.10576515 | 0.08766236 | 0.10762657 | -0.01916629 | 14 / 0.1 s | bdf2, row |
| 65 | 0.00000000 | 0.09912637 | 0.08890169 | 0.10090501 | -0.01906713 | 15 / 0.3 s | bdf2, mg, row |
| 129 | 0.00000000 | 0.09669628 | 0.08923824 | 0.09847520 | -0.01897296 | 15 / 1.3 s | bdf2, mg, row |
| 257 | 0.00000000 | 0.09619828 | 0.08935520 | 0.09798076 | -0.01894692 | 17 / 11.4 s | bdf2, mg, row |
| 513 | 0.00000000 | 0.09614396 | - | - | - | - | - |

### M = 1.5

| grid | CL (a 0) | CD (a 0) | CL (a 1.25) | CD (a 1.25) | CM (a 1.25) | BDF1 it / time a 1.25 | same forces |
|---|---|---|---|---|---|---|---|
| 9 | 0.00000000 | 0.11410997 | 0.05992336 | 0.11556284 | -0.01504732 | 15 / 0.0 s | bdf2, row |
| 17 | 0.00000000 | 0.13101900 | 0.06446640 | 0.13228511 | -0.01467186 | 18 / 0.0 s | bdf2, row |
| 33 | 0.00000000 | 0.10343552 | 0.06699156 | 0.10485700 | -0.01431915 | 16 / 0.1 s | bdf2, row |
| 65 | 0.00000000 | 0.09850019 | 0.06800790 | 0.09990277 | -0.01428968 | 15 / 0.2 s | bdf2, mg, row |
| 129 | 0.00000000 | 0.09632959 | 0.06840162 | 0.09772119 | -0.01418344 | 17 / 1.2 s | bdf2, mg, row |
| 257 | 0.00000000 | 0.09599382 | 0.06844216 | 0.09738849 | -0.01414976 | 16 / 8.9 s | bdf2, mg, row |
| 513 | - | - | - | - | - | - | - |

Note: `grids/2049x2049.x` (82 MB) is not in git: copy it from Moodle into `grids/`.
