/*
   Step 8 check: geometric FAS multigrid (Roe_scheme_accel slides 41-53), multigrid.h / multigrid.cpp.

   M1 coarse grids (slide 41): coarse node (I, J) = fine node (2I, 2J), valid metrics on every level,
      number of levels, an 8 x 8-cell grid is the end of the hierarchy
   M2 restriction of the solution (slide 49): volume-weighted average, keeps a uniform field and
      conserves sum Omega W
   M3 restriction of the residual (slide 49): plain sum of the 4 fine cells, conserves sum R
   M4 prolongation of the correction (slide 50): the 9/16, 3/16, 3/16, 1/16 weights reproduce constants
      and fields linear in the cell indices, periodic across the TE cut, constant across wall / farfield
   M5 full multigrid injection (slide 42): fine cells = their coarse cell
   M6 forcing function (slide 46): on a converged fine solution (R_h = 0) a cycle changes nothing,
      RK and SGS smoothers
   M7 convergence (slides 44, 52-53): same CL, CD, CM as the single-grid explicit solver, far less work;
      W / V cycle, RK / SGS smoother, FMG start (slide 42), on 33x33 and 65x65
   M8 transonic M 0.8: W-cycle converges, same forces as the implicit solver

   Usage: test_multigrid <33x33.x> <65x65.x>
          test_multigrid <33x33.x> <65x65.x> vdebug   (V-cycle diagnostic table only, no checks)
*/

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "implicit.h"
#include "mesh.h"
#include "multigrid.h"
#include "residual.h"
#include "solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using euler::NEQ;

static int g_pass = 0; // number of passed checks
static int g_fail = 0; // number of failed checks

// Print one check line and count it.
static void check(bool ok, const char* fmt, ...)
{
    std::printf(ok ? "  [PASS] " : "  [FAIL] ");
    va_list args;
    va_start(args, fmt);
    std::vprintf(fmt, args);
    va_end(args);
    std::printf("\n");
    if (ok) {
        ++g_pass;
    } else {
        ++g_fail;
    }
}

// Read a plot3D grid and compute its metrics.
static euler::Mesh load_mesh(const std::string& file)
{
    euler::Mesh m = euler::read_plot3d(file);
    euler::compute_metrics(m);
    return m;
}

// Cell array of the mesh, ghost slots NaN like the solver arrays.
static std::vector<double> cell_array(const euler::Mesh& m)
{
    return std::vector<double>(m.vol.size(), std::numeric_limits<double>::quiet_NaN());
}

// Largest |a - b| over the physical cells and the 4 conservative variables.
static double max_diff(const euler::Mesh& m, const euler::Flow& a, const euler::Flow& b)
{
    double d = 0.0; // largest difference
    for (int k = 0; k < NEQ; ++k) {
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                d = std::max(d, std::abs(a.W[k][c] - b.W[k][c]));
            }
        }
    }
    return d;
}

// Wall-clock seconds since t0.
static double since(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

/*
   Multigrid settings of the tests: RK5 smoother of slide 18, 2nd-order column on the fine level,
   1st-order column on the coarse levels (run_multigrid takes fine_rk / coarse_rk as given, the
   default RkScheme would be explicit Euler).
*/
static euler::MultigridSettings mg_settings()
{
    euler::MultigridSettings s;
    s.fine_rk = euler::rk_scheme(euler::Order::Second, 5);
    s.coarse_rk = euler::rk_scheme(euler::Order::First, 5);
    s.print_every = 0;
    return s;
}

// ===================================================================== M1 coarse grids

static void test_coarsening(const euler::Mesh& m33, const euler::Mesh& m65)
{
    // Slide 41: the coarse grid is made of every other node of the fine grid, bit for bit.
    {
        const euler::Mesh c = euler::coarsen_mesh(m65);
        const bool sizes = c.ni == 33 && c.nj == 33 && c.nic == 32 && c.njc == 32; // expected dimensions
        double dmax = 0.0;                                                          // largest node difference
        for (int J = 0; J < c.nj; ++J) {
            for (int I = 0; I < c.ni; ++I) {
                dmax = std::max(dmax, std::abs(c.x[c.node(I, J)] - m65.x[m65.node(2 * I, 2 * J)]));
                dmax = std::max(dmax, std::abs(c.y[c.node(I, J)] - m65.y[m65.node(2 * I, 2 * J)]));
            }
        }
        check(sizes && dmax == 0.0, "M1 65x65 -> %dx%d nodes, coarse node (I,J) = fine node (2I,2J): max diff %.1e",
              c.ni, c.nj, dmax);
    }

    // Metrics of every level of the 65x65 hierarchy: areas and face lengths > 0 and finite.
    {
        const std::vector<euler::Mesh> h = euler::coarse_hierarchy(m65, 4); // coarsest first
        const bool order = h.size() == 3 && h[0].ni == 9 && h[1].ni == 17 && h[2].ni == 33;
        bool ok = order; // all metrics valid
        double fine_area = 0.0; // total area of the 65x65 grid
        for (int j = 0; j < m65.njc; ++j) {
            for (int i = 0; i < m65.nic; ++i) {
                fine_area += m65.vol[m65.cell(i, j)];
            }
        }
        for (const euler::Mesh& c : h) {
            double vmin = 1e300; // smallest cell area
            double area = 0.0;   // total area
            bool faces = true;   // all face lengths > 0 and normals finite
            for (int j = 0; j < c.njc; ++j) {
                for (int i = 0; i < c.nic; ++i) {
                    const double v = c.vol[c.cell(i, j)];
                    vmin = std::min(vmin, std::isfinite(v) ? v : -1.0);
                    area += v;
                }
            }
            for (int j = 0; j < c.njc; ++j) {
                for (int i = 0; i < c.ni; ++i) {
                    const int f = c.node(i, j);
                    faces = faces && c.si[f] > 0.0 && std::isfinite(c.nix[f]) && std::isfinite(c.niy[f]);
                }
            }
            for (int j = 0; j < c.nj; ++j) {
                for (int i = 0; i < c.nic; ++i) {
                    const int f = c.node(i, j);
                    faces = faces && c.sj[f] > 0.0 && std::isfinite(c.njx[f]) && std::isfinite(c.njy[f]);
                }
            }
            ok = ok && vmin > 0.0 && faces;
            std::printf("    level %dx%d: min area %.3e, total area %.6e (65x65: %.6e, rel. diff %.2e)\n", c.ni, c.nj,
                        vmin, area, fine_area, (area - fine_area) / fine_area);
        }
        check(ok, "M1 coarse_hierarchy(65x65, 4) = 9x9, 17x17, 33x33 (coarsest first), all areas and faces > 0");
    }

    // Number of levels: 64 -> 32 -> 16 -> 8 cells, i.e. 4 levels, coarsest 9x9 nodes.
    {
        const euler::MultigridSettings s; // default settings (levels = 0: as many as possible)
        const int n65 = euler::multigrid_levels(m65, s);
        const int n33 = euler::multigrid_levels(m33, s);
        euler::MultigridSettings s2; // levels limited to 2
        s2.levels = 2;
        const int n2 = euler::multigrid_levels(m65, s2);
        const std::vector<euler::Mesh> h = euler::coarse_hierarchy(m65, n65);
        check(n65 == 4 && n33 == 3 && n2 == 2 && h.front().ni == 9 && h.front().nj == 9,
              "M1 multigrid_levels: 65x65 -> %d levels (coarsest %dx%d nodes), 33x33 -> %d, levels = 2 -> %d", n65,
              h.front().ni, h.front().nj, n33, n2);
    }

    // 8 x 8 cells: cannot be coarsened (nic / 2 = 4 < 8), coarsen_mesh throws.
    {
        const euler::Mesh m8 = euler::coarse_hierarchy(m33, 3).front(); // 9x9 nodes
        bool threw = false;
        try {
            euler::coarsen_mesh(m8);
        }
        catch (const std::exception&) {
            threw = true;
        }
        check(m8.nic == 8 && m8.njc == 8 && !euler::can_coarsen(m8) && euler::can_coarsen(m33) && threw,
              "M1 can_coarsen(8x8 cells) = %d (expect 0), can_coarsen(33x33) = %d, coarsen_mesh(8x8) throws: %d",
              euler::can_coarsen(m8) ? 1 : 0, euler::can_coarsen(m33) ? 1 : 0, threw ? 1 : 0);
    }
}

// ===================================================================== M2-M5 transfer operators

// Smooth non-uniform test value of cell (i, j).
static double field(int i, int j)
{
    return 1.0 + 0.5 * std::sin(0.37 * i + 1.1 * j) + 0.01 * j;
}

static void test_transfer(const euler::Mesh& mf)
{
    const euler::Mesh mc = euler::coarsen_mesh(mf);
    std::vector<double> fine = cell_array(mf);   // fine cell array
    std::vector<double> coarse = cell_array(mc); // coarse cell array

    // M2a: a uniform field restricts to the same value (sum Omega w / sum Omega = w).
    {
        const double w = 2.75; // uniform value
        for (int j = 0; j < mf.njc; ++j) {
            for (int i = 0; i < mf.nic; ++i) {
                fine[mf.cell(i, j)] = w;
            }
        }
        euler::restrict_volume_average(mf, mc, fine, coarse);
        double err = 0.0; // largest |coarse - w|
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                err = std::max(err, std::abs(coarse[mc.cell(I, J)] - w));
            }
        }
        check(err < 1e-14, "M2 restrict_volume_average of a uniform field: max error %.1e", err);
    }

    /*
       M2b: conservation of the volume-weighted total (slide 49): with Omega_c = sum of the 4 fine
       areas, sum_c Omega_c W_c = sum_f Omega_f W_f; and each W_c is a convex combination of its 4
       fine values, so it lies between their min and max.
    */
    {
        for (int j = 0; j < mf.njc; ++j) {
            for (int i = 0; i < mf.nic; ++i) {
                fine[mf.cell(i, j)] = field(i, j);
            }
        }
        euler::restrict_volume_average(mf, mc, fine, coarse);
        double sf = 0.0;    // sum over fine cells of Omega w
        double sc = 0.0;    // sum over coarse cells of (sum of the 4 fine areas) w_c
        bool bounded = true; // coarse value within its 4 fine values
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                double vsum = 0.0; // sum of the 4 fine areas
                double lo = 1e300; // smallest fine value
                double hi = -1e300; // largest fine value
                for (int b = 0; b < 2; ++b) {
                    for (int a = 0; a < 2; ++a) {
                        const int c = mf.cell(2 * I + a, 2 * J + b);
                        vsum += mf.vol[c];
                        sf += mf.vol[c] * fine[c];
                        lo = std::min(lo, fine[c]);
                        hi = std::max(hi, fine[c]);
                    }
                }
                const double wc = coarse[mc.cell(I, J)];
                sc += vsum * wc;
                bounded = bounded && wc >= lo - 1e-15 && wc <= hi + 1e-15;
            }
        }
        const double rel = std::abs(sf - sc) / std::abs(sf);
        check(rel < 1e-13 && bounded,
              "M2 restrict_volume_average conserves sum Omega W: %.15e vs %.15e (rel. %.1e), within fine min/max: %d", sc,
              sf, rel, bounded ? 1 : 0);
    }

    // M3: restrict_sum keeps the total of a signed field (sum of the residuals of the 4 fine cells).
    {
        double sf = 0.0;   // sum of the fine values
        double sabs = 0.0; // sum of |fine values| (scale of the round-off)
        for (int j = 0; j < mf.njc; ++j) {
            for (int i = 0; i < mf.nic; ++i) {
                const double v = std::sin(0.37 * i + 1.1 * j) - 0.2;
                fine[mf.cell(i, j)] = v;
                sf += v;
                sabs += std::abs(v);
            }
        }
        euler::restrict_sum(mf, mc, fine, coarse);
        double sc = 0.0;  // sum of the coarse values
        double e4 = 0.0;  // largest |coarse - sum of its 4 fine|
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                const double v4 = fine[mf.cell(2 * I, 2 * J)] + fine[mf.cell(2 * I + 1, 2 * J)] +
                                  fine[mf.cell(2 * I, 2 * J + 1)] + fine[mf.cell(2 * I + 1, 2 * J + 1)];
                sc += coarse[mc.cell(I, J)];
                e4 = std::max(e4, std::abs(coarse[mc.cell(I, J)] - v4));
            }
        }
        check(std::abs(sc - sf) < 1e-13 * sabs && e4 < 1e-15,
              "M3 restrict_sum: sum coarse %.15e = sum fine %.15e (diff %.1e), coarse = sum of its 4 fine (%.1e)", sc,
              sf, sc - sf, e4);
    }

    // M4a: a constant coarse field gives the same constant (weights sum to 16/16).
    {
        const double w = -0.4; // constant value
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                coarse[mc.cell(I, J)] = w;
            }
        }
        euler::prolong_bilinear(mf, mc, coarse, fine);
        double err = 0.0; // largest |fine - w|
        for (int j = 0; j < mf.njc; ++j) {
            for (int i = 0; i < mf.nic; ++i) {
                err = std::max(err, std::abs(fine[mf.cell(i, j)] - w));
            }
        }
        check(err < 1e-15, "M4 prolong_bilinear of a constant: max error %.1e", err);
    }

    /*
       M4b: linear field in the coarse indices, w_c(I, J) = a + b I + c J.
       In index space the coarse cell I covers the fine cells 2I and 2I+1, so its centre is at fine
       position 2I + 1/2, i.e. the fine cell i sits at coarse coordinate X = (i - 1/2) / 2.
       1D weights 3/4 (own cell) and 1/4 (nearest neighbour):
         i = 2I   : 3/4 I + 1/4 (I - 1) = I - 1/4 = (2I - 1/2) / 2
         i = 2I+1 : 3/4 I + 1/4 (I + 1) = I + 1/4 = (2I + 1 - 1/2) / 2
       so the 9/3/3/1 tensor weights reproduce the linear function exactly:
         w_f(i, j) = a + b (i - 1/2) / 2 + c (j - 1/2) / 2,
       for 1 <= i <= nic_f - 2 (away from the periodic cut, where a linear field jumps) and
       1 <= j <= njc_f - 2. At the wall (j = 0) and the farfield (j = njc_f - 1) the neighbour
       row is clamped (copy of the first row): the j part is then the constant J = 0 / J = njc_c - 1.
    */
    {
        const double a = 0.3;  // constant part
        const double b = 0.05; // slope in I
        const double c = -0.02; // slope in J
        for (int J = 0; J < mc.njc; ++J) {
            for (int I = 0; I < mc.nic; ++I) {
                coarse[mc.cell(I, J)] = a + b * I + c * J;
            }
        }
        euler::prolong_bilinear(mf, mc, coarse, fine);
        double ein = 0.0; // largest error inside
        double ebd = 0.0; // largest error on the wall / farfield rows
        for (int j = 0; j < mf.njc; ++j) {
            for (int i = 1; i < mf.nic - 1; ++i) {
                const double X = 0.5 * (i - 0.5); // coarse coordinate of the fine cell
                double Y = 0.5 * (j - 0.5);
                if (j == 0) {
                    Y = 0.0;
                } else if (j == mf.njc - 1) {
                    Y = mc.njc - 1.0;
                }
                const double e = std::abs(fine[mf.cell(i, j)] - (a + b * X + c * Y));
                if (j == 0 || j == mf.njc - 1) {
                    ebd = std::max(ebd, e);
                } else {
                    ein = std::max(ein, e);
                }
            }
        }
        check(ein < 1e-14 && ebd < 1e-14,
              "M4 linear a + b I + c J -> a + b (i-1/2)/2 + c (j-1/2)/2 inside: max error %.1e; "
              "wall / farfield rows (clamped J): %.1e",
              ein, ebd);
    }

    /*
       M4c: periodicity across the TE cut. Coarse delta of 16 in cell (nic_c - 1, J0): the fine cells
       i = 0 (across the cut) see it as their i neighbour,
         fine (0, 2J0) = fine (0, 2J0+1) = 3 (3/16, side neighbour),
         fine (0, 2J0-1) = fine (0, 2J0+2) = 1 (1/16, corner neighbour),
       and its own fine cell (nic_f - 1, 2J0) gets 9. Same with the delta in (0, J0) seen from
       i = nic_f - 1. Each coarse cell spreads 2 x 2 = 4 times its value (sum of the weights it gives).
    */
    {
        const int J0 = mc.njc / 2; // delta row (interior)
        const int nf = mf.nic;     // fine cells in i
        double e1 = 0.0;           // errors, delta at I = nic_c - 1
        double e2 = 0.0;           // errors, delta at I = 0
        double tot = 0.0;          // sum of the fine values (expect 64)
        for (int pass = 0; pass < 2; ++pass) {
            std::fill(coarse.begin(), coarse.end(), 0.0);
            const int Id = pass == 0 ? mc.nic - 1 : 0; // delta column
            coarse[mc.cell(Id, J0)] = 16.0;
            euler::prolong_bilinear(mf, mc, coarse, fine);
            const int ia = pass == 0 ? 0 : nf - 1;         // fine column across the cut
            const int io = pass == 0 ? nf - 1 : 0;         // fine column of the delta's own side
            double e = 0.0;
            e = std::max(e, std::abs(fine[mf.cell(ia, 2 * J0)] - 3.0));
            e = std::max(e, std::abs(fine[mf.cell(ia, 2 * J0 + 1)] - 3.0));
            e = std::max(e, std::abs(fine[mf.cell(ia, 2 * J0 - 1)] - 1.0));
            e = std::max(e, std::abs(fine[mf.cell(ia, 2 * J0 + 2)] - 1.0));
            e = std::max(e, std::abs(fine[mf.cell(io, 2 * J0)] - 9.0));
            e = std::max(e, std::abs(fine[mf.cell(io, 2 * J0 + 1)] - 9.0));
            if (pass == 0) {
                e1 = e;
                for (int j = 0; j < mf.njc; ++j) {
                    for (int i = 0; i < mf.nic; ++i) {
                        tot += fine[mf.cell(i, j)];
                    }
                }
            } else {
                e2 = e;
            }
        }
        check(e1 < 1e-15 && e2 < 1e-15 && std::abs(tot - 64.0) < 1e-13,
              "M4 periodic across the cut: delta at I = nic-1 seen by i = 0 (err %.1e), delta at I = 0 seen by "
              "i = nic-1 (err %.1e), total spread %.3f (expect 64)",
              e1, e2, tot);
    }

    // M5: injection (slide 42) fine cell = its coarse cell, for every variable; restriction undoes it.
    {
        euler::Flow fc; // coarse flow
        euler::Flow ff; // fine flow
        euler::allocate_flow(mc, fc);
        euler::allocate_flow(mf, ff);
        for (int k = 0; k < NEQ; ++k) {
            for (int J = 0; J < mc.njc; ++J) {
                for (int I = 0; I < mc.nic; ++I) {
                    fc.W[k][mc.cell(I, J)] = field(I + 3 * k, J);
                }
            }
        }
        euler::inject_piecewise_constant(mf, mc, fc, ff);
        double e = 0.0;  // largest |fine - its coarse cell|
        double er = 0.0; // largest |restrict(inject(w)) - w|
        for (int k = 0; k < NEQ; ++k) {
            for (int j = 0; j < mf.njc; ++j) {
                for (int i = 0; i < mf.nic; ++i) {
                    e = std::max(e, std::abs(ff.W[k][mf.cell(i, j)] - fc.W[k][mc.cell(i / 2, j / 2)]));
                }
            }
            euler::restrict_volume_average(mf, mc, ff.W[k], coarse);
            for (int J = 0; J < mc.njc; ++J) {
                for (int I = 0; I < mc.nic; ++I) {
                    er = std::max(er, std::abs(coarse[mc.cell(I, J)] - fc.W[k][mc.cell(I, J)]));
                }
            }
        }
        check(e == 0.0 && er < 1e-14,
              "M5 inject_piecewise_constant: fine = its coarse cell (max diff %.1e), restriction of it = coarse (%.1e)",
              e, er);
    }
}

// ===================================================================== M6 forcing consistency

static void test_forcing(const euler::Mesh& m)
{
    const euler::Freestream fs = euler::make_freestream(0.5, 1.25);
    euler::Flow f; // solution
    euler::allocate_flow(m, f);
    euler::init_freestream(f, fs);

    // Converged fine solution: implicit BDF1 + JFNK-GMRES (slides 25, 39-40) to a drop of 1e-12.
    euler::ImplicitSettings is;
    is.tol_drop = 1e-12;
    is.print_every = 0;
    const euler::RunResult ri = euler::run_implicit(m, f, fs, is);
    std::printf("    implicit BDF1 + GMRES: %d it, drop %.2e, %.2f s, CL %.10f CD %.10f\n", ri.iterations,
                ri.res / ri.res0, ri.seconds, ri.forces.cl, ri.forces.cd);
    check(ri.converged, "M6 fine 33x33 solution converged by run_implicit (drop %.1e)", ri.res / ri.res0);
    const euler::Flow fconv = f; // converged solution

    /*
       What the forcing removes: on the coarse grid the restricted converged solution is NOT a solution
       of the coarse (first-order, coarser) equations, R_2h(I W_h) = O(truncation error) != 0. The
       forcing Q = I R_h - R_2h(I W_h) makes the coarse residual R_2h + Q = I R_h = 0 at W_2h^(0).
    */
    {
        const euler::Mesh mc = euler::coarsen_mesh(m);
        euler::Flow fc; // restricted solution on the coarse grid
        euler::allocate_flow(mc, fc);
        euler::init_freestream(fc, fs);
        for (int k = 0; k < NEQ; ++k) {
            euler::restrict_volume_average(m, mc, f.W[k], fc.W[k]);
        }
        euler::apply_bc(mc, fc, fs);
        euler::compute_residual(mc, fc, euler::Order::First);
        const double r2h = euler::density_residual_norm(mc, fc);
        euler::Flow ff = f; // fine residual of the converged solution
        euler::apply_bc(m, ff, fs);
        euler::compute_residual(m, ff, euler::Order::Second);
        const double rh = euler::density_residual_norm(m, ff);
        std::printf("    ||R_h(W_h)|| = %.2e, coarse ||R_2h(I W_h)|| without forcing = %.2e (truncation-error "
                    "difference the forcing cancels)\n",
                    rh, r2h);
    }

    // Reference scale: what 4 single-grid RK5 iterations change (2 cycles x (n_pre + n_post)).
    {
        euler::Flow fe = fconv; // copy for the explicit iterations
        euler::ExplicitSettings es;
        es.rk = euler::rk_scheme(euler::Order::Second, 5);
        es.max_iter = 4;
        es.tol_drop = 0.0;
        euler::run_explicit(m, fe, fs, es);
        std::printf("    4 single-grid RK5 iterations from it: max |dW| = %.2e\n", max_diff(m, fe, fconv));
    }

    // 2 W-cycles, RK smoother, then the SGS smoother: the solution must not move.
    for (int sm = 0; sm < 2; ++sm) {
        euler::MultigridSettings s = mg_settings();
        s.smoother = sm == 0 ? euler::MgSmoother::RK : euler::MgSmoother::SGS;
        s.gamma = 2;
        s.max_iter = 2;
        s.tol_drop = 0.0;
        f = fconv;
        const euler::RunResult r = euler::run_multigrid(m, f, fs, s);
        const double d = max_diff(m, f, fconv);
        check(d < 1e-10, "M6 forcing (slide 46): 2 W-cycles (%s smoother, %d levels) on the converged solution: "
                         "max |dW| = %.2e (expect ~0)",
              sm == 0 ? "RK" : "SGS", euler::multigrid_levels(m, s), d);
        (void)r;
    }
}

// ===================================================================== M7 / M8 convergence

// One line of the result table.
struct Row {
    std::string name;     // variant
    int cycles = 0;       // multigrid cycles (or single-grid iterations) done
    int fine_its = 0;     // fine-grid smoothing iterations (RK5 or implicit steps)
    double seconds = 0.0; // wall-clock time, setup and FMG included
    double drop = 0.0;    // final residual drop
    double res = 0.0;     // final ||R_rho||
    bool converged = false; // drop reached
    euler::Forces F;      // forces of the final solution
};

static void print_row(const Row& r)
{
    std::printf("    %-22s %7d %8d %9.2f %10.2e %10.2e %3s %15.10f %15.10f %15.10f\n", r.name.c_str(), r.cycles,
                r.fine_its, r.seconds, r.drop, r.res, r.converged ? "yes" : "NO", r.F.cl, r.F.cd, r.F.cm);
}

static void print_header()
{
    std::printf("    %-22s %7s %8s %9s %10s %10s %3s %15s %15s %15s\n", "variant", "cycles", "fine it", "seconds",
                "drop", "||R||", "cvg", "CL", "CD", "CM");
}

// Multigrid run from the free stream.
static Row run_mg(const euler::Mesh& m, const euler::Freestream& fs, const euler::MultigridSettings& s,
                  const std::string& name)
{
    euler::Flow f; // solution
    euler::allocate_flow(m, f);
    euler::init_freestream(f, fs);
    const auto t0 = std::chrono::steady_clock::now();
    const euler::RunResult r = euler::run_multigrid(m, f, fs, s);
    Row row;
    row.name = name;
    row.seconds = since(t0);
    row.cycles = r.converged ? r.iterations - 1 : r.iterations; // the last record stops before its cycle
    row.fine_its = row.cycles * (s.n_pre + s.n_post);
    row.drop = r.res / r.res0;
    row.res = r.res;
    row.converged = r.converged;
    row.F = r.forces;
    print_row(row);
    return row;
}

// Single-grid explicit RK5 run from the free stream (reference).
static Row run_sg(const euler::Mesh& m, const euler::Freestream& fs, double tol)
{
    euler::Flow f; // solution
    euler::allocate_flow(m, f);
    euler::init_freestream(f, fs);
    euler::ExplicitSettings es;
    es.rk = euler::rk_scheme(euler::Order::Second, 5);
    es.tol_drop = tol;
    es.max_iter = 400000;
    const auto t0 = std::chrono::steady_clock::now();
    const euler::RunResult r = euler::run_explicit(m, f, fs, es);
    Row row;
    row.name = "single grid RK5";
    row.seconds = since(t0);
    row.cycles = r.converged ? r.iterations - 1 : r.iterations;
    row.fine_its = row.cycles;
    row.drop = r.res / r.res0;
    row.res = r.res;
    row.converged = r.converged;
    row.F = r.forces;
    print_row(row);
    return row;
}

// Do the forces of a run agree with the reference to tol?
static bool same_forces(const Row& a, const Row& ref, double tol, double& dmax)
{
    dmax = std::max({std::abs(a.F.cl - ref.F.cl), std::abs(a.F.cd - ref.F.cd), std::abs(a.F.cm - ref.F.cm)});
    return dmax < tol;
}

static void test_convergence(const euler::Mesh& m, const std::string& gname)
{
    const euler::Freestream fs = euler::make_freestream(0.5, 1.25);
    const double tol = 1e-10;  // residual drop of every run
    const double ftol = 1e-8;  // agreement of CL, CD, CM with the single-grid solution
    std::printf("\n  M7 %s, M 0.5, alpha 1.25, drop %.0e, %d levels (2nd order fine, 1st order coarse, RK5 of slide 18)\n",
                gname.c_str(), tol, euler::multigrid_levels(m, mg_settings()));
    print_header();
    const Row sg = run_sg(m, fs, tol);

    std::vector<Row> rows; // multigrid variants
    {
        euler::MultigridSettings s = mg_settings();
        s.tol_drop = tol;
        s.gamma = 2;
        rows.push_back(run_mg(m, fs, s, "MG W RK"));
        s.gamma = 1;
        rows.push_back(run_mg(m, fs, s, "MG V RK"));
        s.gamma = 2;
        s.fmg_cycles = 10;
        rows.push_back(run_mg(m, fs, s, "MG W RK + FMG 10"));
        s.fmg_cycles = 0;
        s.smoother = euler::MgSmoother::SGS;
        s.sgs_cfl = 20.0;
        rows.push_back(run_mg(m, fs, s, "MG W SGS (cfl 20)"));
        s.gamma = 1;
        rows.push_back(run_mg(m, fs, s, "MG V SGS (cfl 20)"));
    }

    check(sg.converged, "M7 %s single-grid RK5 converged: %d it, %.2f s", gname.c_str(), sg.cycles, sg.seconds);
    for (const Row& r : rows) {
        if (r.name == "MG V RK") {
            /*
               Not a check: with the 2nd-order fine / 1st-order coarse operators the V-cycle with the
               RK smoother stalls at a spurious state (see the vdebug runs: it converges with the same
               order on every level, with n_pre or n_post = 0, or with half the fine CFL).
            */
            std::printf("  [INFO] M7 %s %-18s %s after %d cycles (drop %.2e, CL %.4f) - operator mismatch, "
                        "see V-cycle consistency check below\n",
                        gname.c_str(), r.name.c_str(), r.converged ? "converged" : "did NOT converge", r.cycles,
                        r.drop, r.F.cl);
            continue;
        }
        double d = 0.0; // largest force difference
        const bool same = same_forces(r, sg, ftol, d);
        check(r.converged && same, "M7 %s %-18s converged in %4d cycles (rate %.3f / cycle), %.2f s; max |dC| vs single grid %.1e",
              gname.c_str(), r.name.c_str(), r.cycles, std::pow(r.drop, 1.0 / std::max(r.cycles, 1)), r.seconds, d);
    }
    const Row& w = rows[0]; // W-cycle RK
    check(w.converged && 5 * w.fine_its < sg.fine_its && w.seconds < sg.seconds,
          "M7 %s W-cycle RK vs single grid: %d vs %d fine RK5 iterations (%.1fx fewer), %.2f vs %.2f s (%.1fx faster)",
          gname.c_str(), w.fine_its, sg.fine_its, static_cast<double>(sg.fine_its) / std::max(w.fine_its, 1),
          w.seconds, sg.seconds, sg.seconds / std::max(w.seconds, 1e-9));

    /*
       V-cycle machinery with consistent operators (first order on every level, RK5 first-order
       column): V and W cycles must converge to the same R_h = 0 solution (slide 46: the converged
       solution does not depend on the cycle).
    */
    {
        std::printf("    first order on every level:\n");
        euler::MultigridSettings s = mg_settings();
        s.tol_drop = tol;
        s.fine_order = euler::Order::First;
        s.fine_rk = euler::rk_scheme(euler::Order::First, 5);
        s.gamma = 1;
        const Row v1 = run_mg(m, fs, s, "MG V RK o1/o1");
        s.gamma = 2;
        const Row w1 = run_mg(m, fs, s, "MG W RK o1/o1");
        double d = 0.0; // force difference V vs W
        const bool same = same_forces(v1, w1, ftol, d);
        check(v1.converged && w1.converged && same,
              "M7 %s V-cycle RK, 1st order on every level: %d cycles (W: %d), max |dC| V vs W %.1e", gname.c_str(),
              v1.cycles, w1.cycles, d);
    }
}

static void test_transonic(const euler::Mesh& m)
{
    const euler::Freestream fs = euler::make_freestream(0.8, 1.25);
    std::printf("\n  M8 33x33, M 0.8, alpha 1.25, drop 1e-10\n");
    print_header();

    // Reference: implicit BDF1 + GMRES.
    Row ref;
    {
        euler::Flow f; // solution
        euler::allocate_flow(m, f);
        euler::init_freestream(f, fs);
        euler::ImplicitSettings is;
        is.tol_drop = 1e-12;
        is.print_every = 0;
        const auto t0 = std::chrono::steady_clock::now();
        const euler::RunResult r = euler::run_implicit(m, f, fs, is);
        ref.name = "implicit BDF1 GMRES";
        ref.seconds = since(t0);
        ref.cycles = r.converged ? r.iterations - 1 : r.iterations;
        ref.fine_its = ref.cycles;
        ref.drop = r.res / r.res0;
        ref.res = r.res;
        ref.converged = r.converged;
        ref.F = r.forces;
        print_row(ref);
    }

    euler::MultigridSettings s = mg_settings();
    s.tol_drop = 1e-10;
    s.max_iter = 3000;
    s.gamma = 2;
    const Row w = run_mg(m, fs, s, "MG W RK");
    s.smoother = euler::MgSmoother::SGS;
    const Row ws = run_mg(m, fs, s, "MG W SGS (cfl 20)");

    double d = 0.0; // force difference with the implicit solution
    const bool same = ref.converged && same_forces(w, ref, 1e-7, d);
    check(w.converged && same, "M8 transonic W-cycle RK converged in %d cycles, %.2f s, CL %.8f CD %.8f; vs implicit %.1e",
          w.cycles, w.seconds, w.F.cl, w.F.cd, d);
    double ds = 0.0; // same for the SGS smoother
    const bool same_s = ref.converged && same_forces(ws, ref, 1e-7, ds);
    check(ws.converged && same_s, "M8 transonic W-cycle SGS converged in %d cycles, %.2f s; vs implicit %.1e", ws.cycles,
          ws.seconds, ds);
}

/*
   Diagnostic runs of the V-cycle with the RK smoother (third argument "vdebug", first grid only,
   M 0.5, alpha 1.25, 600 cycles at most): which ingredient makes the default V-cycle stall
   (number of levels, smoothing counts, order of the fine / coarse operators, CFL, wall pressure,
   number of RK stages). No checks, a table only.
*/
static void debug_v(const euler::Mesh& m)
{
    const euler::Freestream fs = euler::make_freestream(0.5, 1.25);
    print_header();
    struct Case {
        const char* name; // variant
        int levels;       // number of levels
        int gamma;        // 1 V, 2 W
        int n_pre;        // pre-smoothing
        int n_post;       // post-smoothing
        int n_coarse;     // coarsest-level iterations
        int fo;           // fine order (1 or 2)
        int co;           // coarse order (1 or 2)
        double ccfl;      // coarse CFL factor (x the table value)
        double fcfl;      // fine CFL factor
        bool wext;        // extrapolated wall pressure
        int fst;          // fine RK stages
        int cst;          // coarse RK stages
    };
    const Case cases[] = {
        {"V lev3", 3, 1, 1, 1, 2, 2, 1, 1.0, 1.0, true, 5, 5},
        {"V lev3 o1/o1", 3, 1, 1, 1, 2, 1, 1, 1.0, 1.0, true, 5, 5},
        {"V lev3 o2/o2", 3, 1, 1, 1, 2, 2, 2, 1.0, 1.0, true, 5, 5},
        {"V lev3 ccfl0.5", 3, 1, 1, 1, 2, 2, 1, 0.5, 1.0, true, 5, 5},
        {"V lev3 fcfl0.5", 3, 1, 1, 1, 2, 2, 1, 1.0, 0.5, true, 5, 5},
        {"V lev3 mirror wall", 3, 1, 1, 1, 2, 2, 1, 1.0, 1.0, false, 5, 5},
        {"V lev3 RK3/RK3", 3, 1, 1, 1, 2, 2, 1, 1.0, 1.0, true, 3, 3},
        {"V lev2 nc20", 2, 1, 1, 1, 20, 2, 1, 1.0, 1.0, true, 5, 5},
        {"V lev2 nc20 o2/o2", 2, 1, 1, 1, 20, 2, 2, 1.0, 1.0, true, 5, 5},
        {"V lev2 nc20 o1/o1", 2, 1, 1, 1, 20, 1, 1, 1.0, 1.0, true, 5, 5},
        {"V lev2 nc1", 2, 1, 1, 1, 1, 2, 1, 1.0, 1.0, true, 5, 5},
        {"W lev3 o1/o1", 3, 2, 1, 1, 2, 1, 1, 1.0, 1.0, true, 5, 5},
    };
    for (const Case& c : cases) {
        euler::MultigridSettings s = mg_settings();
        s.tol_drop = 1e-10;
        s.max_iter = 600;
        s.levels = c.levels;
        s.gamma = c.gamma;
        s.n_pre = c.n_pre;
        s.n_post = c.n_post;
        s.n_coarse = c.n_coarse;
        s.fine_order = c.fo == 1 ? euler::Order::First : euler::Order::Second;
        s.coarse_order = c.co == 1 ? euler::Order::First : euler::Order::Second;
        s.fine_rk = euler::rk_scheme(s.fine_order, c.fst);
        s.coarse_rk = euler::rk_scheme(s.coarse_order, c.cst);
        s.fine_rk.cfl *= c.fcfl;
        s.coarse_rk.cfl *= c.ccfl;
        euler::WALL_PRESSURE_EXTRAPOLATE = c.wext;
        const Row r = run_mg(m, fs, s, c.name);
        (void)r;
        euler::WALL_PRESSURE_EXTRAPOLATE = true;
    }
}

int main(int argc, char** argv)
{
    if (argc >= 4 && std::string(argv[3]) == "vdebug") {
        try {
            debug_v(load_mesh(argv[1]));
        }
        catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 1;
        }
        return 0;
    }
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <33x33.x> <65x65.x>\n", argv[0]);
        return 1;
    }
    try {
        const euler::Mesh m33 = load_mesh(argv[1]);
        const euler::Mesh m65 = load_mesh(argv[2]);
        if (m33.ni != 33 || m65.ni != 65) {
            throw std::runtime_error("expected a 33x33 and a 65x65 grid");
        }
        std::printf("\n=== M1 coarse grids\n");
        test_coarsening(m33, m65);
        std::printf("\n=== M2-M5 transfer operators on 33x33\n");
        test_transfer(m33);
        std::printf("\n=== M2-M5 transfer operators on 65x65\n");
        test_transfer(m65);
        std::printf("\n=== M6 forcing consistency on 33x33\n");
        test_forcing(m33);
        std::printf("\n=== M7 convergence\n");
        test_convergence(m33, "33x33");
        test_convergence(m65, "65x65");
        std::printf("\n=== M8 transonic\n");
        test_transonic(m33);
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
