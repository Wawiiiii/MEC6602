# HW2 – 2D Euler solver: progress notes

Last updated: 2026-11-11, after step 2 (metrics) was finished and verified.

This file sums up everything decided so far and where the code stands, so the work
can be picked up again without the original conversation.

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
- **Prior solvers are references only, not ported:**
  - `/apps/partage-elns/1_RawiH/1_Euler_Solver_V2` (C++ JST + RK5 + IRS, FAS multigrid,
    multizone OpenMP). Known issues there: RK coefficients mismatched with JST, no ghost
    cells, IRS applied to dt*R, array-of-structs data, exceptions in hot loops.
  - `/apps/partage-elns/1_RawiH/GPU_JAX_Solver/FICAS` (JAX Cartesian IBM). Useful ideas:
    implicit SGS smoother on a first-order Roe Jacobian at CFL ~100 (explicit + FAS took
    1630 cycles at M0.5 vs 250 with implicit). Its numbers are NOT a validation reference
    (first-order immersed-boundary wall).
  - Non-slide tricks from those codes (farfield vortex, extrapolated wall pressure,
    Swanson-Turkel, smooth HLLEC, ...) are not used.
- Language/tools: C++17, header + cpp per module, CMake. GCC 11.5, CMake 3.26.
  Machine: 20 threads, 62 GB RAM, no GPU. No `-ffast-math` (it would break the NaN
  checks on unfilled ghost cells).

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

---

## 4. Code layout and conventions

```
HW2/
  CMakeLists.txt        library euler2d (src/*.cpp) + one test executable per step
  include/mesh.h        Mesh struct, read_plot3d(), compute_metrics()
  src/mesh.cpp
  include/output.h      write_tecplot_mesh()
  src/output.cpp
  tests/test_mesh.cpp   step 1 check
  tests/test_metrics.cpp step 2 check
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

### Build and run (from `13_MEC6602/`)
```bash
cmake -S HW2 -B HW2/build && cmake --build HW2/build -j
./HW2/build/test_mesh    HW2/grids/9x9.x
./HW2/build/test_metrics HW2/grids/9x9.x HW2/output/mesh_9x9.dat
```

---

## 6. Next steps (roadmap)

3. **Second-order Roe flux with MUSCL** (Roe_scheme_accel slides 5-16):
   Roe averages and |dF1|, |dF2,3,4|, |dF5| (slides 7-9), Harten entropy fix with
   delta = c~/10 (slide 11), per-face algorithm (slide 12), MUSCL + limiter with
   eps = Omega^2 (slides 13-16).
   **Decision needed first:** limiter van Albada (kappa = 0) or Hemker-Koren (kappa = 1/3),
   and which variables are reconstructed (as the slides say).
4. **Boundary conditions with ghost cells** (basic_algo slides 13-16): slip wall
   (reflect velocity), farfield supersonic + subsonic characteristic, i-cut copy.
   Sanity check from basic_algo slide 20: all BCs farfield + freestream init -> nothing moves.
5. **Explicit Euler, then multistage RK** with local time stepping
   (basic_algo slides 17-19; upwind column of the table).
6. **Outputs:** CL, CD, CM (quarter chord), surface Cp, Tecplot field (cell -> node
   interpolation, Euler2D slide 20), L2 density residual + Cauchy force criterion.
   **Open item:** reference chord 1.0 or 1.00893 for Vassberg & Jameson comparison.
7. **JFNK GMRES:** Jv ~ (R(W + eps v) - R(W)) / eps, restart 20-50, preconditioner =
   first-order Roe Jacobian blocks + SGS, CFL ramp sigma^n = sigma^0 ||R^0|| / ||R^n||
   (Roe_scheme_accel slides 25-27, 39-40).
8. **Multigrid** (FAS, slides 41-53). How it combines with JFNK still to decide.
9. **Multizone + OpenMP** (slides 57-61) and **fixed-CL** (post.pdf p.9).

---

## 7. To resume with Claude

Point Claude at this file and the memory notes, e.g.:
"Read HW2/PROGRESS.md. We are at step 3 (second-order Roe + MUSCL). Keep the same
rules: slides only, step by step, you code, comment style as described, tell me when
to check."
