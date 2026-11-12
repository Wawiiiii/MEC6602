# HW2 – 2D Euler solver: progress notes

Last updated: 2026-11-12, step 5 (explicit Euler / RK5) written and tested, waiting for the user's check.

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
  include/output.h      write_tecplot_mesh()
  src/output.cpp
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

### Build and run (Windows, from `HW2/`)
```bash
cmake -S . -B build && cmake --build build --config Release
./build/Release/test_mesh.exe    grids/9x9.x
./build/Release/test_metrics.exe grids/9x9.x output/mesh_9x9.dat
./build/Release/test_roe.exe     grids/9x9.x grids/17x17.x grids/33x33.x grids/65x65.x grids/129x129.x grids/257x257.x
./build/Release/test_bc.exe      grids/9x9.x grids/17x17.x grids/33x33.x grids/65x65.x grids/129x129.x grids/257x257.x
./build/Release/test_explicit.exe grids/33x33.x grids/65x65.x
./build/Release/euler.exe cases/example.txt         # solver, all settings in a case file (keys documented there)
./build/Release/euler.exe grids/65x65.x 0.5 1.25   # or on the command line: --order --stages --cfl --max_iter --tol --history --print_every
```

---

## 6. Next steps (roadmap)

3. ~~Second-order Roe flux with MUSCL~~ written (see §5), awaiting the user's check.
4. ~~Boundary conditions with ghost cells~~ written (see §5), awaiting the user's check.
5. ~~Explicit Euler, then multistage RK~~ written (see §5), awaiting the user's check.
6. **Outputs:** CL, CD, CM (quarter chord), surface Cp, Tecplot field (cell -> node
   interpolation, Euler2D slide 20), L2 density residual + Cauchy force criterion.
   **Open item:** reference chord 1.0 or 1.00893 for Vassberg & Jameson comparison.
7. **Fully implicit solver (Part 2), JFNK GMRES** (Roe_scheme_accel slides 25-40). Order agreed
   2026-11-12: RK (5) -> outputs (6) -> implicit (7). The converged RK solution is the reference the
   implicit solver must reproduce (same R = 0, slide 27). Sub-steps, each tested:
   7a. First-order Roe Jacobian blocks dR/dW (slides 26-29), checked against finite differences
       of the first-order residual. Open: how |A| = T|Lambda|T^-1 is built; BCs in the operator
       or neglected (slide 30).
   7b. Implicit Euler (slide 25) solved with block SGS / LU-SGS sweeps (slides 33-35, the
       professor's recommendation, slide 32), 2-3 forward/backward sweeps, under-relaxation.
   7c. GMRES (slide 39) with the Jacobian-free product Jv ~ (R(W + eps v) - R(W)) / eps (slide 40),
       restart 20-50, preconditioner = 7a/7b, CFL ramp sigma^n = sigma^0 ||R^0|| / ||R^n||.
8. **Multigrid** (FAS, slides 41-53). How it combines with JFNK still to decide.
9. **Multizone + OpenMP** (slides 57-61) and **fixed-CL** (post.pdf p.9).

---

## 7. To resume with Claude

Point Claude at this file and the memory notes, e.g.:
"Read HW2/PROGRESS.md. We are at step 3 (second-order Roe + MUSCL). Keep the same
rules: slides only, step by step, you code, comment style as described, tell me when
to check."
