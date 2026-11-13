# HW2 – 2D Euler solver: progress notes

Last updated: 2026-11-13 evening, steps 6 (outputs) and 7 (implicit BDF1/BDF2/ROW + GMRES) done,
step 8 (FAS multigrid V/W, RK or SGS smoother, FMG) written and partly working; grid study on
64..1024 shows order >= 2 for CL, CM (CD: see RESULTS.md); all numbers in RESULTS.md.

This file sums up everything decided so far and where the code stands, so the work
can be picked up again without the original conversation.
Every deviation / interpretation of the slides and every known risk is listed in
`DEVIATIONS.md` (keep it up to date at each step; it feeds the report).

---

## 1. The homework (HW2/docs/MEC6602E_HW2.pdf)

2D steady Euler, finite volume, NACA0012 on a single-block O-grid (plot3D files,
9x9 to 2049x2049 in `HW2/grids/`).

Outputs required for every run:
- CL, CD, CM (moment about the quarter chord)
- Tecplot field file
- surface Cp file

Parts:
1. Second-order finite-volume scheme with explicit Runge-Kutta (multistage).
2. An acceleration technique.
3. Grid convergence study up to 2049², compared with Vassberg & Jameson 2010
   (FLO82 tables, post.pdf p.29). Order of accuracy p by Newton iteration
   (post.pdf p.26-28).
4. Subsonic, transonic and supersonic cases.
5. Fixed CL: M = 0.8, CL = 0.4, secant iteration on alpha (post.pdf p.9).

Convergence monitoring: L2 norm of the density residual + Cauchy criterion on the forces.

---

## 2. Decisions taken (do not re-open)

- **Follow the slides exactly, no deviation.** Grading says any deviation from the
  class algorithms must be justified; we want none. When the slides give options,
  pick among those; if something is not in the slides, flag it, do not import it.
- **Second order only.** The spatial scheme must be second order.
  **Hard requirement (user, 2026-11-13): the grid convergence study (GCI / order p by Newton
  iteration, post.pdf slides 26-28) must show p ~ 2 for the forces, like the Vassberg & Jameson
  reference.** Status: CD order 65->129 ~ 1.75 with RK (DEVIATIONS A14 table); to be measured on
  65/129/257 with the implicit solver, and fixed if < 2 (suspects: wall pressure B6, limiter B1).
- **Spatial scheme:** Roe flux, second order with MUSCL + limiter
  (Roe_scheme_accel.pdf slides 5-16).
  - Limiter: **Hemker-Koren, kappa^ = 1/3**, smooth eps-form of slide 16 (decided 2026-11-12).
    eps = **(K h)^3, h = sqrt(Omega), K = 5** (Venkatakrishnan) of the cell the state comes from,
    NOT the slide's suggested Omega^2, which makes the residual stall at 1e-3 on 65x65 and 129x129
    (changed 2026-11-12, numbers in DEVIATIONS.md A14). Not switched off at extrema (A13).
  - If Roe ever becomes too troublesome, the user accepts switching to JST (central + scalar
    dissipation, basic_algo slides 8-10). Not needed so far.
  - Reconstructed variables: **primitive (rho, u, v, p)** (decided 2026-11-12; slides only say "U",
    to be justified in one line in the report).
  - Harten entropy fix delta = c~/10 on all three eigenvalues (slide 11).
  - Slide 16 lists a_R = D+U_{I+1}, b_R = D-U_{I+1}. For Hemker-Koren (not symmetric in a, b) that
    is inconsistent with the slide 15 psi(r) form; the code uses a = U_{I+1} - U_I (across the
    face), b = U_{I+2} - U_{I+1} for the right state = slide 15 exactly (mirror of the left state).
    Checked by test_roe B2/B3 and independently by a second agent. Mention in the report.
- **Time / solution methods, in this order:**
  1. explicit multistage RK with local time stepping (required by Part 1),
  2. Newton-Krylov GMRES (JFNK, Roe_scheme_accel slides 39-40),
  3. multigrid (FAS, slides 41-53),
  4. multizone (slides 57-61, OpenMP) + fixed-CL.
- **Why Roe + RK + implicit together:** Roe is the space discretisation (how the
  residual R(W) is computed). RK and implicit/JFNK are ways to march to R(W) = 0.
  The homework demands explicit RK; the implicit solver is the acceleration.
  Using a first-order Jacobian/preconditioner in the implicit operator does NOT
  change the steady-state accuracy (slide 27): the converged solution satisfies
  the second-order R(W) = 0.
- **Boundary conditions** (basic_algo slides 7, 13-16), decided 2026-11-12:
  - Physical boundary faces (wall j = 0, farfield j = nj-1) use the flux of the averaged state
    Fc(1/2 (W_ghost + W_interior)), conservative average (slide 7 method #2, "to be used at physical
    boundaries"). The i-cut is connectivity: Roe + MUSCL like an interior face.
  - Wall (slide 14): ghost = interior rho, p, velocity mirrored. Slide gives rho0 = rho2, p0 = p2 for
    the second ghost and no velocity: the second ghost is set = the first one.
  - Farfield (slides 15-16), normal pointing out of the domain, rho0 c0 from the interior cell.
    Inflow/outflow from the sign of the interior V.n; supersonic/subsonic from |V.n| vs c0 (normal Mach,
    the characteristic count; the slides do not say which Mach number). Ghost = 2 Wb - Wd in
    conservative variables (slide 16 "extrapolated to the ghost centre"), so the face average is exactly
    Wb; second farfield ghost = first one (never read by the residual).
  - Non-dimensional free stream: rho = 1, p = 1/gamma (c = 1), u = M cos(alpha), v = M sin(alpha).
  - Not used (not in the slides): farfield vortex correction (post.pdf slide 22 says "not to be done
    in homework"), extrapolated wall pressure.
- **Prior solvers are references only, not ported.** The user's earlier solver
  `Desktop/Euler_Solver` (NOT "Euler Solver V2", which is wrong) may be used as a reference; things
  can be ported only if they match the slides / homework guidelines.
  - `/apps/partage-elns/1_RawiH/1_Euler_Solver_V2` (C++ JST + RK5 + IRS, FAS multigrid,
    multizone OpenMP). Known issues there: RK coefficients mismatched with JST, no ghost
    cells, IRS applied to dt*R, array-of-structs data, exceptions in hot loops.
  - `/apps/partage-elns/1_RawiH/GPU_JAX_Solver/FICAS` (JAX Cartesian IBM). Useful ideas:
    implicit SGS smoother on a first-order Roe Jacobian at CFL ~100 (explicit + FAS took
    1630 cycles at M0.5 vs 250 with implicit). Its numbers are NOT a validation reference
    (first-order immersed-boundary wall).
  - Non-slide tricks from those codes (farfield vortex, extrapolated wall pressure,
    Swanson-Turkel, smooth HLLEC, ...) are not used.
- Language/tools: C++17, header + cpp per module, CMake. No `-ffast-math` (it would break the
  NaN checks on unfilled ghost cells).
  - Steps 1-2 were done on Linux (GCC 11.5, CMake 3.26, 20 threads, 62 GB RAM).
  - From step 3 on: Windows 11, Visual Studio 2022 (MSVC) + CMake 4.2. Anaconda Python has numpy
    (`/c/Users/rawih/Anaconda3/python.exe`); msys2 g++ also available. No PDF renderer: slides read
    with `pdftotext`.
  - Only grids 9x9 .. 257x257 are on this machine; 513, 1025, 2049 must be downloaded (Moodle).

---

## 3. Working rules

- Go **step by step, very slowly**. One step per round: Claude writes the step,
  builds it, runs its check, then says "you can go and check". Next step only on "go".
- Comment well: what the code does and which slide it comes from.
- Comment style:
  - `//` for one-line comments and comments at the end of a code line,
  - multi-line: `/*` on its own line, text indented with NO leading stars, `*/` under `/*`,
  - one variable per line, with its comment on the same line,
  - CMake keeps `#`.
- If the user has edited a file, keep their edits.
- Build and run tests (C++ and Python) without asking; test everything thoroughly.
- A second agent may be used freely, e.g. for independent cross-checks.

---

## 4. Code layout and conventions

```
HW2/
  CMakeLists.txt        library euler2d (src/*.cpp) + one test executable per step
  include/mesh.h        Mesh struct, read_plot3d(), compute_metrics()
  src/mesh.cpp
  include/output.h      write_tecplot_mesh(), write_tecplot_field(), write_surface_cp()
  src/output.cpp
  include/forces.h      ForceRef, Forces, wall_pressure(), integrate_forces(), compute_forces()
  src/forces.cpp
  include/blockmat.h    BlockMatrix: 4x4 blocks D, W, E, S, N per physical cell (implicit operator)
  include/jacobian.h    flux_jacobian(), roe_abs_matrix(), roe_face_jacobians(), assemble_first_order_jacobian()
  src/jacobian.cpp
  include/linsolve.h    block_matvec(), factor_diagonal(), block_sgs(), gmres() (+ dot, norm2, axpy)
  src/linsolve.cpp
  include/implicit.h    TimeScheme {BDF1, BDF2, ROW}, LinearSolver {GMRES, SGS}, ImplicitSettings, run_implicit()
  src/implicit.cpp
  include/multigrid.h   MultigridSettings, MgSmoother {RK, SGS}, coarsen_mesh(), restrictions, prolongation,
  src/multigrid.cpp     run_multigrid() (V / W cycles), coarse_hierarchy(), inject_piecewise_constant() (FMG)
  RESULTS.md            all results so far (grid study, methods, M 0.8, M 1.5) -- start here for the report
  results/              raw outputs (runs_bdf1/, mg_*, M05_a125_513/ with field + plots, case files)
  post/plot_results.py  Python plots: residual / force history, surface Cp, Mach / Cp / entropy contours, mesh
  post/order_study.py   observed order p (Newton, post.pdf slides 26-28) and extrapolated CL, CM from surface files
  include/flow.h        NEQ, GAMMA, Flow struct (W, Q, R), prim <-> cons
  src/flow.cpp          allocate_flow(), compute_primitives()
  include/roe.h         physical_flux(), harten(), roe_flux()
  src/roe.cpp
  include/muscl.h       hk_slope(), muscl_left(), muscl_right() (inline)
  include/residual.h    Order {First, Second}, i_face_flux(), j_face_flux(), compute_residual()
  src/residual.cpp
  include/bc.h          Freestream, make_freestream(), init_freestream(), wall_ghost(),
  src/bc.cpp            farfield_state(), JminBc {Wall, Farfield}, apply_bc()
  include/solver.h      RkScheme, rk_scheme() (slide 18 table), local_time_step(),
  src/solver.cpp        density_residual_norm(), ExplicitSettings, run_explicit()
  apps/euler_main.cpp   the solver program: euler <case.txt>  or  euler <mesh.x> <mach> <alpha_deg> [--key value]
  cases/example.txt     example case file (key = value, all settings documented)
  tests/test_mesh.cpp   step 1 check
  tests/test_metrics.cpp step 2 check
  tests/test_roe.cpp    step 3 check
  tests/test_bc.cpp     step 4 check
  tests/test_explicit.cpp step 5 check
  tests/test_output.cpp step 6 check
  tests/test_jacobian.cpp step 7a check
  tests/test_linsolve.cpp step 7b check
  tests/test_implicit.cpp step 7c check
  output/               Tecplot files written by the tests
  grids/*.x             plot3D meshes; grids/tecplot/*.dat = Tecplot copies (X, Y, Area)
  docs/                 slides, homework statement, prof's code_example.txt
```

Grid facts (checked):
- File order: i fastest, then j; all x then all y.
- j = 0 is the wall, j = nj-1 the farfield (~148-149 chords away).
- i = 0 and i = ni-1 are the same nodes (cut at the trailing edge, x = 1.00893).
- Wall order in i: TE -> lower surface -> LE -> upper surface -> TE (clockwise).
  With j going outward, the (i, j) cell corners turn counter-clockwise -> positive areas.

Storage (mesh.h):
- `NG = 2` ghost layers on every side (Euler2D slide 18).
- One flat `std::vector<double>` per quantity, i contiguous.
- `node(i, j) = (j + NG) * ni_tot + (i + NG)`, `cell(i, j) = (j + NG) * nic_tot + (i + NG)`;
  i, j can be negative (ghosts). Ghost slots start as NaN so unfilled use shows up.
- Sizes: `ni, nj` nodes; `nic = ni-1, njc = nj-1` cells; `*_tot = n + 2*NG`.

Metrics (Euler2D slide 21):
- `vol[cell(i,j)]` = cell area = 0.5 * (a x b), a, b the diagonals; sign checked > 0.
- i-face (i, j): edge node(i,j) -> node(i,j+1), left cell (i-1, j), right cell (i, j),
  data `nix, niy, si` stored at `node(i, j)`.
- j-face (i, j): edge node(i,j) -> node(i+1,j), left cell (i, j-1), right cell (i, j),
  data `njx, njy, sj` stored at `node(i, j)`.
- Normals are unit vectors from left to right (+i or +j). At the wall (j = 0) the
  normal points into the fluid.
- Residual assembly convention (basic_algo slide 6): one flux per face,
  `R[left] += F*dS`, `R[right] -= F*dS`. Boundary faces feed ghost slots, which are ignored.
- The i-cut needs no special face: face i = 0 and face i = ni-1 are the same physical
  face, each seen from one side, with the ghost cells holding copies across the cut.

---

## 5. Status

### Step 1 – mesh reader + data structure: DONE
- `read_plot3d()` reads the file, checks the block count, sizes, early end of file,
  extra values, and that the i-cut closes (within 1e-10).
- `test_mesh` on 9x9: nodes 9x9 (13x13 with ghosts), cells 8x8 (12x12);
  index checks 0 / 168 / 0 / 143; farfield radius 148.067-149.043.
  2049² reads in 1.45 s, 69 MB.

### Step 2 – metrics + slide 22 checks: DONE
- `compute_metrics()`, `write_tecplot_mesh()`, `test_metrics`.
- Cartesian 3x3 test: areas 1, lengths 1, i normals (1,0), j normals (0,1), residual 0.
- Airfoil grids 9x9 .. 2049²: all areas > 0; sum of areas = shoelace ring area
  (rel. diff 8e-15); | |n| - 1 | <= 2.2e-16; constant flux residual / perimeter <= 9e-17;
  airfoil area 0.0822 (matches NACA0012).
- Tecplot files with areas in `HW2/output/mesh_*.dat`.

### Step 3 – second-order Roe flux + MUSCL / Hemker-Koren + residual: WRITTEN, awaiting check
- `roe_flux()` (slides 7-12), `muscl_left/right()` (slides 13-16), `i_face_flux()`,
  `j_face_flux()`, `compute_residual()` (one flux per face, + left / - right).
- `compute_metrics()` now also fills the ghost-cell areas (needed by the limiter eps):
  copies across the i-cut, mirror at wall and farfield. Steps 1-2 checks unchanged.
- `test_roe` on 9x9 .. 257x257: 47/47 checks pass.
  - Roe: consistency 0, conservation 7e-15, supersonic upwinding (property 1) 7e-14,
    rotation 2e-14, RH shock exact without the fix (V~ - c~ = 0), slip line exact without the fix.
  - MUSCL: linear data exact, mirror symmetry, slide 16 = slide 15 form, step kept, face values
    second order (L1 2.24-2.53).
  - Residual: uniform flow 9e-16 (slide 12), i-cut fluxes identical, conservation 4e-16.
  - Residual vs. exact divergence of a smooth field: 1st order L1 order 0.98, 2nd order L1 order
    2.18 / 1.94. L2 order of the 2nd-order scheme drops to 1.70 on 257: limiter at smooth extrema
    (lines of cells with O(h) error -> O(h^1.5) in L2). With the limiter off (eps -> inf) L2 goes
    back to ~2. Not a bug; keep in mind for the grid convergence analysis.
  - Smooth HK form: on monotone data the face value never leaves its two cells (B5, 47th check);
    at a sharp extremum it can overshoot (0.5, 1, 0 -> 1.125). A TVD switch was tried and rejected
    (6x larger smooth error), see DEVIATIONS.md A13.
  - One 2nd-order residual on 257x257: ~4.8 ms (73 ns/cell, serial).
- Independent check (second agent, Python 40-digit eigen-decomposition |A| = T|L|T^-1 with the
  Harten fix): roe_flux matches to ~4e-14 over 4014 cases incl. transonic / sonic expansions.

### Step 4 – boundary conditions with ghost cells: WRITTEN, awaiting check
- `apply_bc()` fills every ghost the residual reads: i-cut copies, wall mirror (slide 14),
  farfield characteristic state (slides 15-16) extrapolated to the ghost. `j_face_flux()` now uses
  the averaged-state flux on the wall and farfield faces (slide 7 #2).
- `test_bc` on 9x9 .. 257x257: 53/53 checks pass.
  - Wall ghost: rho, p kept, V.n reversed, V.t kept; wall flux = (0, p n, 0); mass/energy through
    the wall on the grids <= 6e-16.
  - Farfield: free stream kept in all regimes; subsonic in/outflow satisfy the 1D characteristic
    relations they come from (independent check of the slide formulas) to 4e-16; supersonic copies.
  - Slide 20 sanity check (all BCs farfield + free stream): residual <= 8e-16 ("nothing moves").
  - No NaN reaches the residual from a NaN-initialised ghost layer; i-cut fluxes identical.
  - Free stream + real wall: only row j = 0 moves at first order, rows 0-1 at second order.
  - Symmetric field on the (exactly symmetric) grids -> mirror-symmetric residual to 3e-16.
- `test_roe` still passes (47/47) with the new boundary faces.

### Step 5 – explicit Euler + multistage RK with local time stepping: WRITTEN, awaiting check
- `rk_scheme()`: slide 18 table (1st- and 2nd-order columns, 3-5 stages) + explicit Euler (slide 17).
  Default: 5 stages, 2nd-order column, CFL 1.15.
- `local_time_step()`: dt = sigma Omega / (Lambda_I + Lambda_J), normals and lengths averaged over the
  two faces of each direction (slides 11, 17). dt frozen over the stages.
- `run_explicit()`: apply_bc before every residual (slide 9 step 4d), stops on a drop of the density
  residual norm (post.pdf slide 4, ||R_rho/Omega||_2), writes a history file, NaN -> "diverged".
- The first runs exposed a stall of the 2nd-order residual (1e-3) on 65 and 129 -> limiter eps changed
  to Venkatakrishnan's (5h)^3 (DEVIATIONS.md A14, full comparison table there).
- `test_explicit` on 33x33 and 65x65: 12/12 checks pass (table coefficients, dt by hand, one Euler
  iteration by hand, slide 20 over 200 RK5 iterations 3e-15, alpha = 0 mirror-symmetric to 6e-15 and
  CL 2e-15, machine-zero convergence 1st and 2nd order). test_roe 47/47, test_bc 53/53 still pass.
- Convergence, 2nd order RK5 CFL 1.15 (drop 1e-12 unless noted):
  | case | 33x33 | 65x65 | 129x129 (drop 1e-10) |
  |---|---|---|---|
  | M 0.5 alpha 1.25, 1st order | 1101 it, 0.3 s | 3505 it, 3.7 s | 7806 it, 35 s |
  | M 0.5 alpha 1.25, 2nd order | 11506 it, 5 s | 34938 it, 60 s | ~76000 it, 516 s |
  | M 0.8 alpha 1.25 (drop 1e-10) | | 20108 it, 34 s | |
  | M 1.2 / 1.5 alpha 0 (drop 1e-10) | | 4506 it, 8 s / 1610 it, 3 s | |
- Preliminary forces (wall-face pressure, chord 1; real outputs are step 6), M 0.5 alpha 1.25, 2nd
  order: CL 0.1652 / 0.1710 / 0.1753, CD 0.0188 / 0.0081 / 0.0024 on 33 / 65 / 129 (CD should -> 0,
  order ~1.75). Total enthalpy error 2.2 % / 1.5 % (33 / 65).
- Explicit RK is slow on fine grids (iterations ~ x2 per refinement): this is what the implicit solver
  (step 7) is for.
- (Numbers above are with the slide-14 mirror wall; the default wall BC changed in step 7, A20.)

### Step 6 – outputs: DONE (2026-11-13), awaiting check
- Forces (post.pdf 16-18): CL, CD, CM about (xref, yref) = quarter chord, chord 1 by default; wall
  pressure = averaged-state pressure of the wall face = the momentum flux of the residual (A16).
- Surface Cp file (x/c, y/c, Cp, p per wall face; comment line with CL CD CM), Tecplot field file
  (node values = average of the 4 surrounding cells, ghosts and corners included; X Y Rho U V P Mach
  Cp Entropy), history file (iteration, ||R_rho||, drop, CL, CD, CM, seconds).
- Cauchy criterion on CL, CD, CM (post.pdf slide 5): case keys cauchy_window, cauchy_tol.
- Corner ghost cells are now filled by apply_bc (B10 fixed).
- `test_output`: 30/30 on 9, 33, 65, 257 (closed-contour identities, linear pressure: F exact, moment
  O(h^2) with err * nic^2 = 3.2 on every grid, forces = wall momentum fluxes, field = free stream with
  all-farfield BCs, no seam at the cut, alpha = 0 -> CL = CM = 0, Cauchy stop with CL within 5e-8).
- Python: `post/plot_results.py --history h.dat --surface s.dat --field f.dat --out plots` (overlays:
  repeat the options); `post/order_study.py s65.dat s129.dat ... --cells 64 128 ...`.

### Step 7 – implicit solver (Part 2): DONE (2026-11-13), awaiting check
- 7a `jacobian.cpp`: first-order Roe Jacobian, frozen |A~| (slide 26), Harten fix included, wall and
  farfield faces by finite differences of the exact boundary flux lumped on the diagonal (slide 30),
  plus the north block for the extrapolated wall pressure. test_jacobian 19/19 (M v vs FD 1e-10).
- 7b `linsolve.cpp`: block SGS (forward + backward sweeps, explicit 4x4 inverses) and restarted
  right-preconditioned GMRES (MGS + Givens). test_linsolve 38/38; SGS(2) cuts GMRES iterations ~8x.
- 7c `implicit.cpp`: BDF1 (slide 25), BDF2 and ROS2 Rosenbrock-Wanner (user request, A17); linear
  solver GMRES (JFNK, J v by finite difference of the full 2nd-order residual, slide 40) or SGS only
  (slide 33); CFL ramp sigma0 ||R0||/||R|| (slide 40) clamped to [cfl0 = 1, 1e8] (A21); update limiter
  (A18). test_implicit 7/7: BDF1/BDF2/ROW = RK5 solution to 4e-13, alpha = 0 symmetric, M 0.8 and 1.5
  converge, SGS-only converges at CFL <= 50.
- Performance (BDF1 + GMRES, M 0.5, drop 1e-12, serial): 65x65 18 it 0.3 s; 129 22 it 2.3 s;
  257 31 it 17 s; 513 38 it 174 s (vs explicit RK5: 65x65 34938 it 60 s).
- **Wall BC changed to the extrapolated ghost pressure (A20)** because with slide 14's P1 = P2 the CM
  converged at order ~1 and CL at ~1.5 (user requirement: order 2).

### Grid convergence study (Part 3), M = 0.5, alpha = 1.25, BDF1 + GMRES, final defaults
| grid (cells) | CL | CD | CM | implicit it | time |
|---|---|---|---|---|---|
| 64 | 0.17580183 | 0.00734443 | -0.00227543 | 18 | 0.3 s |
| 128 | 0.17762119 | 0.00166278 | -0.00208267 | 22 | 2.3 s |
| 256 | 0.17905556 | 0.00028266 | -0.00221650 | 31 | 17 s |
| 512 | 0.17930803 | 0.00005915 | -0.00224442 | 38 | 174 s |
| 1024 | 0.17935083 | 0.00002699 | -0.00224975 | 53 | 1164 s |
| 2048 | (running) | | | | |
Observed order (post/order_study.py, Newton on G(p)):
| grids | p CL | CL_inf | p CD (3 grids) | CD_inf | p CM | CM_inf | p CD = log2(CDm/CDf) |
|---|---|---|---|---|---|---|---|
| 128-256-512 | 2.51 | 0.179362 | 2.63 | 1.6e-5 | 2.26 | -0.0022518 | 2.26 |
| 256-512-1024 | 2.56 | 0.179360 | 2.80 | 2.2e-5 | 2.39 | -0.0022510 | 1.13 |
64-128-256 is pre-asymptotic for CL and CM (CM not monotone there).
Drag: with the exact value 0 assumed (slide 28), p drops to 1.13 on 512->1024 because CD tends to a
small non-zero limit ~2e-5 (three-grid estimate, p ~ 2.8). It is lift-related: at alpha = 0 (no lift)
CD = 1.64e-3 / 2.60e-4 / 3.85e-5 on 128 / 256 / 512 -> p = 2.65, 2.76, no floor (1024 running).
Explanation: the farfield is ~150 chords away with no vortex correction (post.pdf slide 22 says not
to do it in the homework), so the circulation of the lifting airfoil is not represented at the
boundary: an O(CL^2 / R) error that does not go away with mesh refinement. To mention in the report.
For comparison with the slide-14 mirror wall (same grids 128-256-512): CL p = 1.49, CM p ~ 1.0.

### Step 8 – FAS multigrid (slides 41-53): WRITTEN 2026-11-13, partly working
- `multigrid.cpp`: coarse grids = every other node (65 -> 4 levels, 2049 -> 9 levels down to 9x9);
  restriction volume-weighted (W) and summed (R); forcing Q = I(R + Q) - R_coarse(I W) (slides 46, 51);
  bilinear 9-3-3-1 prolongation of the correction (slide 50); V or W cycle (slides 52-53);
  pre / post / coarsest smoothing counts; FMG start inside the cycle (mg_fmg) and an FMG start for
  any other scheme (fmg_levels, slide 42); correction limiter (max_change).
  Smoother: RK (slide 18 multistage) or SGS (one implicit Euler step + block SGS on the level's
  first-order Jacobian). Coarse levels first order (mg_coarse_order, default 1). Case: scheme = mg.
- test_multigrid (written by an agent, built in `agent_mg/`, now also a target of the main CMake):
  37/37 (transfer operators exact, converged fine solution unchanged by a cycle to 2e-15 (slide 46),
  all converged variants = single-grid solution to 1e-9).
- Works: W-cycle + SGS smoother at M 0.5 on 65-257 (42-48 cycles; 257: 10 s vs BDF1 17 s);
  W + SGS at M 0.8 on 257 (67 cycles, 13 s). W + RK on 65 (174 cycles, 1.2 s vs RK5 60 s).
- Not working yet: RK smoother W-cycle on 257 (stalls at 1.3e-2); V-cycle + RK on 33/65 (stalls on a
  wrong state, CL -1.0; fine = 2nd order, coarse = 1st order mismatch: converges with 1st order on all
  levels); W + SGS at M 0.8 on 65 and 129 (stalls 4e-2 / 5e-3). To investigate: 2nd-order coarse levels
  with RK, fewer levels for transonic, lower fine CFL in the cycle.
- All results: `RESULTS.md`.

### OpenMP (slides 57-70), 2026-11-13
- Library linked with OpenMP (MSVC OpenMP 2.0, `find_package(OpenMP)`); case key `threads` (0 = all).
- Residual restructured: all face fluxes first (stored in f.Fi / f.Fj, parallel over rows), then each
  cell = sum of its 4 faces (parallel). No races, and the result is bit-identical for any number of
  threads (checked: 1 vs 8 threads give the same residual history and forces to the last digit).
- Also parallel: primitives, time step, RK updates, apply_bc columns, Jacobian assembly (i-faces by
  rows, j-faces by columns, boundaries by columns), block mat-vec, diagonal inversion, implicit copies.
  Residual norm: row sums in parallel, rows added in order (deterministic). Serial: block SGS sweeps
  (Gauss-Seidel is sequential), GMRES dot products (kept deterministic).
- Speed-up (machine busy with other runs, 8 threads): one 2nd-order residual on 257x257 5.2 -> 1.2 ms;
  explicit RK 129x129 3000 it 21.3 -> 4.9 s (4.3x); implicit BDF1 + GMRES 257x257 20.3 -> 14.6 s (1.4x,
  the serial SGS preconditioner now dominates). All test suites pass with OpenMP.
- Possible next step: multicolour (red-black / diagonal) ordering to parallelise the SGS sweeps.

### Build and run (Windows, from `HW2/`)
```bash
cmake -S . -B build && cmake --build build --config Release
./build/Release/test_mesh.exe    grids/9x9.x
./build/Release/test_metrics.exe grids/9x9.x output/mesh_9x9.dat
./build/Release/test_roe.exe     grids/9x9.x grids/17x17.x grids/33x33.x grids/65x65.x grids/129x129.x grids/257x257.x
./build/Release/test_bc.exe      grids/9x9.x grids/17x17.x grids/33x33.x grids/65x65.x grids/129x129.x grids/257x257.x
./build/Release/test_explicit.exe grids/33x33.x grids/65x65.x
./build/Release/test_output.exe  grids/9x9.x grids/33x33.x grids/65x65.x grids/257x257.x
./build/Release/test_jacobian.exe grids/33x33.x
./build/Release/test_linsolve.exe
./build/Release/test_implicit.exe grids/33x33.x
./build/Release/euler.exe cases/example.txt         # solver, all settings in a case file (keys documented there)
./build/Release/euler.exe grids/65x65.x 0.5 1.25 --scheme bdf1   # or any key on the command line as --key value
python post/plot_results.py --history history.dat --surface surface_cp.dat --field field.dat --out plots
python post/order_study.py s65.dat s129.dat s257.dat s513.dat --cells 64 128 256 512
```
Fine grids: use gmres_restart 50, gmres_max_iter 150 (513) to 300 (2049). 2049x2049 needs ~9 GB.
`euler.exe` cannot be rebuilt while it runs; a second build folder (`build/dev`) is handy for that.

---

## 6. Next steps (roadmap)

3. ~~Second-order Roe flux with MUSCL~~ written (see §5), awaiting the user's check.
4. ~~Boundary conditions with ghost cells~~ written (see §5), awaiting the user's check.
5. ~~Explicit Euler, then multistage RK~~ written (see §5), awaiting the user's check.
6. ~~Outputs~~ done (see §5). **Open item:** reference chord 1.0 or 1.00893 for the Vassberg &
   Jameson comparison (case key `chord`).
7. ~~Implicit solver~~ done (see §5): BDF1 / BDF2 / ROW + GMRES or SGS. Original plan below.
   Remaining: Part 3 report numbers on 1025 / 2049; Part 4 (sub/trans/supersonic analysis on a fine
   grid); Part 5 fixed CL (M 0.8, CL 0.4, secant on alpha, post.pdf slide 9) -> step 8.
   Original plan:
   **Fully implicit solver (Part 2), JFNK GMRES** (Roe_scheme_accel slides 25-40). Order agreed
   2026-11-12: RK (5) -> outputs (6) -> implicit (7). The converged RK solution is the reference the
   implicit solver must reproduce (same R = 0, slide 27). Sub-steps, each tested:
   7a. First-order Roe Jacobian blocks dR/dW (slides 26-29), checked against finite differences
       of the first-order residual. Open: how |A| = T|Lambda|T^-1 is built; BCs in the operator
       or neglected (slide 30).
   7b. Implicit Euler (slide 25) solved with block SGS / LU-SGS sweeps (slides 33-35, the
       professor's recommendation, slide 32), 2-3 forward/backward sweeps, under-relaxation.
   7c. GMRES (slide 39) with the Jacobian-free product Jv ~ (R(W + eps v) - R(W)) / eps (slide 40),
       restart 20-50, preconditioner = 7a/7b, CFL ramp sigma^n = sigma^0 ||R^0|| / ||R^n||.
8. **Fixed CL** (Part 5, post.pdf slide 9): secant iteration on alpha inside the iterations.
9. Optional: **multigrid** (FAS, slides 41-53), **OpenMP** (slides 57-70; MSVC = OpenMP 2.0).

---

## 7. To resume with Claude

Point Claude at this file and the memory notes, e.g.:
"Read HW2/PROGRESS.md, HW2/DEVIATIONS.md and HW2/RESULTS.md. Steps 1-7 are done, multigrid (8) is
written and partly working. Next: fix the multigrid stalls (RK smoother, transonic coarse grids),
2049 result, Part 4 (M 0.5 / 0.8 / 1.5 analysis, shock quality), Part 5 fixed CL, chord question.
Keep the same rules: slides unless a documented reason, order 2 required, you code and test,
comment style as described."
