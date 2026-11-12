# HW2 – deviations, interpretations and possible failure modes

Companion to `PROGRESS.md`. Two lists:

- **Part A**: every place where the code is not a word-for-word reading of the slides (an
  interpretation, a choice the slides leave open, or a small addition). Each one needs a line
  of justification in the report (the grading asks for it).
- **Part B**: things that could go wrong or degrade the results, with how to test them later.

Status: `open` = not tested yet, `checked` = a test exists and passes, `n/a yet` = the feature
it depends on is not written yet.

Last updated: 2026-11-12, after step 5.

---

## Part A – deviations and interpretations

| # | Where | What the slides say | What the code does | Why | Status |
|---|---|---|---|---|---|
| A1 | `muscl.h` | Roe slide 16: right state uses a_R = D+U_{I+1}, b_R = D-U_{I+1} | Right state uses a = U_{I+1} - U_I, b = U_{I+2} - U_{I+1} | The slide 16 labels, taken literally with Hemker-Koren, disagree with slide 15's own psi(r) form (up to 36 % of the jump) and drop the smooth-data accuracy from 3rd to 2nd order. The code is exactly slide 15, and the right state is the mirror of the left. Van Albada is symmetric in a, b, so the labels only matter for HK. | checked (test_roe B2/B3, second agent) |
| A2 | `residual.cpp` | Slides 13-16: "U" | MUSCL on the primitive variables (rho, u, v, p) | The slides do not say which variables. Primitive variables keep the reconstructed rho and p closer to physical values next to shocks (decided 2026-11-12). | decision |
| A3 | `mesh.cpp` | Slide 16: eps from the cell volume | Ghost cells get an area too (copy across the cut, mirror at wall/farfield) | The left/right state of a boundary face comes from a ghost cell, which needs an eps. The mirror ghost has the area of its mirror cell. | checked (test_metrics, test_roe) |
| A4 | `residual.cpp` | — | Boundary faces (wall, farfield) use Fc(1/2 (W_g + W_i)), not Roe | This is slide 7 method #2 ("to be used at physical boundaries, except the connectivity") and slide 13. Listed because the Roe slides themselves say nothing about boundaries. | checked (test_bc E4/E5) |
| A5 | `residual.cpp` | Slide 7: W01 = 1/2 (W0 + W1) | The average is taken in conservative variables | W is the conservative vector on the slides. Side effect: at the wall p01 = p_i + (gamma-1)/2 rho_i V_n,i^2 (see B6). | decision |
| A6 | `bc.cpp` | Slide 14: rho0 = rho2, p0 = p2 for the 2nd ghost; no velocity given | 2nd wall ghost = 1st wall ghost (velocity reflected too) | The slide is incomplete for the velocity. The 2nd wall ghost is not read by the current residual (boundary faces use A4, face j = 1 reads ghost 1 only). | decision |
| A7 | `bc.cpp` | Slides 15-16: "supersonic / subsonic inflow / outflow" | Inflow/outflow from the sign of V.n of the interior cell, supersonic from \|V.n\| >= c0 (normal Mach) | The slides do not say which velocity or Mach number. The normal Mach number is what sets how many characteristics enter. | decision |
| A8 | `bc.cpp` | Slide 16: "extrapolated from the interior and the BC to the ghost cell centres" | W_g = 2 Wb - Wd in conservative variables; 2nd farfield ghost = 1st | Linear extrapolation from the cell centre (d) through the face (b). In conservative variables the face average is exactly Wb. The 2nd farfield ghost is never read. | checked (test_bc E5) |
| A9 | `bc.h` | post.pdf 13-15: codes are non-dimensional (no specific scaling) | rho_inf = 1, p_inf = 1/gamma (c_inf = 1), \|V_inf\| = M | Any consistent scaling works. Cp and forces must use q_inf = 1/2 rho_inf V_inf^2 = M^2/2. | decision |
| A10 | `roe.cpp` | Slide 11: delta "a small fraction of the local speed of sound (like c~/10)" | delta = c~/10 on all three eigenvalues, including \|V~\| | Exactly the slide. Listed because applying it to \|V~\| smears contacts / slip lines (B4). | decision |
| A11 | `bc.cpp` | basic_algo slide 14 title: "No-slip wall" | Slip wall (V.n = 0, tangential velocity kept) | The formulas on the slide are the slip wall; no-slip has no meaning for Euler. | decision |
| A12 | — | post.pdf slide 22: vortex correction "not to be done in homework" | Not used | Following the slide. Consequence: see B11. | decision |
| A13 | `muscl.h` | Slide 16 smooth limiter formula (slide 14: "no local extrema created") | Used as written, no extra switch-off when a b < 0 | Tried 2026-11-12: the TVD switch (psi = 0 for r <= 0) made the residual error on a smooth field 6x larger (257: 6.3e-4 vs 9.6e-5) and the 2nd-order L2 order 1.5 instead of 1.7. Reason: eps = Omega^2 ~ h^4 ~ a^2 at a smooth extremum, so eps is what keeps smooth extrema from being clipped while discontinuities (a^2 >> eps) are fully limited; the same idea as van Albada / Venkatakrishnan's eps ~ dx^3. Kept the slide form. On monotone data the face value never leaves its two neighbour cells (test_roe B5); at a sharp extremum it can overshoot by up to 0.127 x the local jump (B2). | checked (test_roe B5) |

| A14 | `muscl.h` `limiter_eps()` | Slide 16: "eps is a small number. We can use the square of the cell volume" | eps = (K h)^3, h = sqrt(Omega), K = 5 (Venkatakrishnan, J. Comput. Phys. 118, 1995) | With eps = Omega^2 the 2nd-order residual **stalls** at a drop of ~1e-3 on 65x65 and 129x129 (limit cycle of the limiter in the wall row near the leading edge; converges on 33x33). This limiter-induced stall is exactly what Venkatakrishnan's eps was designed to remove. Measured at M = 0.5, alpha = 1.25 (table below): with (5h)^3 every grid converges to machine zero, the spurious drag is lower, and the step 3 smooth-field test improves (257: L2 order 1.86 instead of 1.69, error 7.0e-5 instead of 1.0e-4). Limiting at shocks is kept (there a^2 = O(1) >> eps). | checked (test_explicit G6, test_roe C4) |
| A15 | `solver.cpp` | post.pdf slide 4: "L2-norm of the density residual" (no normalisation given) | \|\|R_rho\|\| = sqrt(mean over cells of (R_rho / Omega)^2), i.e. d(rho)/dt | R / Omega is the rate of change of rho in the cell; without the division the huge farfield cells would weigh more than the airfoil cells. Only used for monitoring / stopping (drop relative to the first iteration). | decision |

Limiter eps comparison (2026-11-12), NACA0012 M = 0.5 alpha = 1.25, 2nd order, RK5 CFL 1.15,
preliminary CD from the wall-face pressure (should be 0 in subsonic flow, so CD = error):

| Variant | CD 33x33 | CD 65x65 | CD 129x129 | order 65->129 | residual | CL 129 |
|---|---|---|---|---|---|---|
| slide: eps = Omega^2, averaged-state wall flux | 0.0302 | 0.0113 | 0.0031 | (1.9) | **stalls** at 1e-3 on 65, 129 | 0.1747 |
| eps = Omega | 0.0220 | 0.0088 | 0.0024 | 1.84 | machine zero | 0.1753 |
| **eps = (5h)^3 (chosen)** | 0.0188 | 0.0081 | 0.0024 | 1.75 | machine zero | 0.1753 |
| eps = Omega^2 + Roe flux at the wall (mirror states) | 0.0529 | 0.0163 | 0.0044 | 1.91 | machine zero | 0.1753 |
| no limiter (reference only, not usable with shocks) | 0.0122 | 0.0037 | 0.0010 | 1.87 | machine zero | 0.1768 |
| first order | 0.0836 | 0.0561 | 0.0327 | 0.78 | machine zero | 0.1651 |

Also tried on 33/65 and rejected: pressure-only wall flux with the wall pressure extrapolated
linearly from the first two cells (CL +0.004, CD unchanged), Harten fix on the acoustic waves only
(no effect on CD). Conclusions: (1) every 2nd-order variant converges at order ~1.8-1.9, the scheme
is second order; (2) the limiter is the main source of the spurious drag on these coarse grids
(wall spacing 0.049 chord on 65x65 vs leading-edge radius 0.016); (3) keeping the slide 7 averaged-
state wall flux is better than the Roe wall flux.

Review 2026-11-12 (user: "correct them if they are false and large deviations from the state of
the art"): every item above was checked against the slides, Blazek and usual practice.
- A1 is already the consistent reading (the slide labels are the error, not the code).
- A2-A9, A11, A12 are either exactly the slide or the usual choice where the slides are silent.
- A10 (fix on |V~| too) is the slide and Blazek; it costs some dissipation on entropy/shear waves.
  Kept; measured through B4 (subsonic drag) before any change.
- A13: the only candidate for a fix. The fix was tried and was worse (numbers in A13): reverted.
Nothing changed in the code from this review.

Still to decide (later steps):
- Reference chord for CL, CD, CM: 1.0 or the grid's 1.00893 (the TE of the grids is at x = 1.00893).
- Implicit solver: how |A| is built in the first-order Jacobian (slide 26: T|Lambda|T^-1), and
  whether the boundary conditions enter the implicit operator (slide 30) or are neglected.

---

## Part B – possible failure modes and how to test them

| # | Risk | Why it could happen | How to test later | Status |
|---|---|---|---|---|
| B1 | Second-order solution converges at less than order 2 | Seen in the step 3 test: L2 order 1.70 on 257 (limiter at smooth extrema, L1 still 2). The wall pressure at the cell centre (B6) and the TE kink (B9) can also cost accuracy. | V&J grid study (Part 3): order p of CL, CD, CM by Newton iteration (post.pdf 26-28), subsonic drag = error (post.pdf 28). | open |
| B2 | Oscillations or negative p / rho near shocks | The smooth HK formula overshoots at a sharp extremum (up to 0.127 x the local jump; data 0.5, 1, 0 -> 1.125), see A13. No positivity check on the reconstructed states, so sqrt of a negative number in the Roe average gives NaN. | Transonic M = 0.8 alpha = 1.25 and supersonic M = 1.5 runs: watch the Cp at the shock, min p and min rho of the face states, NaN check each iteration. If oscillations show up, fallback = the TVD switch of A13 (cost known: 6x larger smooth error). | open |
| B3 | Residual stalls before machine zero | A limiter can create a limit cycle. | **Happened** with eps = Omega^2 on 65 and 129 (drop 1e-3); fixed by A14. Now M = 0.5 converges to 1e-12 on 33 and 65 (test_explicit G6), 1e-10 on 129. To recheck on 257 and in the transonic case at every new feature (implicit solver too). | checked (A14) |
| B4 | Extra spurious entropy / drag | Harten fix on \|V~\| (A10) adds dissipation on the entropy and shear waves; delta = c~/10 is a tuning value. | Subsonic case: CD should go to 0 with refinement; total-pressure loss along the wall; compare delta = c~/10 with a smaller value (test only). | open |
| B5 | No limiting in the outer field | eps = (5h)^3 (A14) grows with the cell size: on the outer cells (h ~ 3-30) eps >> (jump)^2 and the slope becomes the unlimited kappa = 1/3 slope (b + 2a)/3, while near the airfoil eps is small (limiting only at discontinuities, see A13). Fine for smooth farfield flow, but a shock that reaches the outer cells (supersonic cases) would not be limited. | Supersonic M = 1.5: check the bow shock and the outer cells for oscillations. | open |
| B6 | Wall pressure only first-order accurate pointwise | The wall flux uses the pressure of the averaged state, i.e. the first-cell value (plus (gamma-1)/2 rho V_n^2, A5), not the pressure extrapolated to the wall. | Order of CL/CM in the grid study; compare the surface Cp from the cell values with a linear extrapolation to the wall (post-processing only). | open |
| B7 | Wrong farfield regime, or switching near \|V.n\| = c | At M around 1 the sides of the farfield circle have \|V.n\| close to c0; the regime can change from one iteration to the next. | M = 1.0 - 1.2 runs: count the regime changes per iteration; check the residual convergence. | open |
| B8 | Bad farfield ghost (negative rho or p) | W_g = 2 Wb - Wd can be unphysical when a strong wave reaches the farfield (early iterations). Ghost 1 is read by the MUSCL stencil of the last interior face. | NaN / positivity check of the ghost values in the first iterations of a cold start at M = 0.8 and 1.5. | open |
| B9 | Loss of accuracy at the trailing edge | The grids have a constant i-spacing jump of 2.76 in the first wall cell at the TE (does not go away with refinement); the MUSCL works in index space. The sharp TE sets the Kutta condition numerically. | Cp near the TE vs V&J; local error of the step 3 C4 test with / without the TE cells (done: no visible effect). | open |
| B10 | Corner ghost cells are NaN | The 4 corner ghost blocks (ghost in i and in j) are never filled. Any later stencil that reads them gives NaN. | The cell-to-node interpolation for Tecplot (step 6) at node (0,0) and (ni-1,0) reads a corner cell; the implicit Jacobian may too. Must be handled explicitly. | n/a yet |
| B11 | CL differs from V&J by a farfield effect | No vortex correction (A12); the farfield is at about 150 chords. | Compare CL with the V&J tables (post.pdf 29); estimate the 1/R effect if the difference is systematic. | open |
| B12 | Explicit RK unstable at the table CFL | The table CFL (slide 18) is for 1D model problems; the boundary flux (A4) and the farfield extrapolation change the stability locally. | RK5 at the table CFL (2.5 first order, 1.15 second order) is stable and converges at M = 0.5, 0.8, 1.2, 1.5 on 65x65. Largest stable CFL not searched yet. | partly checked |
| B13 | Implicit solver does not converge or stalls | First-order approximate Jacobian + boundary conditions neglected in the operator (slide 30) reduce stability; SGS does not converge the linear system (slide 33). | CFL ramp history, linear residual per GMRES iteration, comparison of the converged solution with the RK one (must be the same R = 0 solution). | n/a yet |
| B14 | Memory / time on 2049x2049 | 4.2 M cells. Now: W, Q, R = 12 arrays = 0.4 GB + mesh 0.3 GB; a stored first-order block Jacobian (5 blocks of 4x4 per cell) = 2.7 GB; GMRES with restart 30 = 4 GB. Machine: 31 GB, 16 threads. One 2nd-order residual ~ 0.3 s serial. | Run the 2049 grid once each solver works; time per iteration and peak memory. | n/a yet |
| B15 | OpenMP limits on MSVC | MSVC supports OpenMP 2.0 only (signed loop index, no collapse, limited reductions). | When the multizone / OpenMP step starts; results must be bit-identical with 1 thread. | n/a yet |
| B16 | Tests do not cover what they seem to | The step 3 order test (C4) fills the ghosts by analytic extrapolation, not with the real BCs, so it says nothing about the boundary accuracy. | The real accuracy check is B1 (grid study). | open |
| B17 | Ghost-cell mirror metrics used where real geometry matters | Ghost areas are copies (A3); only eps uses them now. An implicit Jacobian or a time step evaluated in a ghost cell would get a mirror area, not a real one. | Check what each new feature reads from ghost cells. | open |
