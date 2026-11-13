/*
   Step 7 check: implicit solvers (Roe_scheme_accel slides 25-40) against the explicit RK solution.

   K1 BDF1, BDF2 and ROW with JFNK-GMRES converge to machine zero in few iterations and give the
      same CL, CD, CM as explicit RK5 (slide 27: the steady solution does not depend on the
      implicit operator, only on R(W) = 0)
   K2 alpha = 0: the implicit solution stays mirror-symmetric, CL = CM = 0
   K3 transonic M = 0.8 and supersonic M = 1.5 converge
   K4 first-order residual (the Jacobian is then the exact one of the residual): Newton converges in
      very few iterations at a large CFL
   K5 SGS-only linear solver (slide 33): converges at a moderate CFL (information: iterations)

   Usage: test_implicit <mesh.x>   (33x33 recommended: the RK reference takes ~5 s there)
*/

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "implicit.h"
#include "mesh.h"
#include "solver.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>

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

// Largest mirror asymmetry of W: cell (i, j) vs (nic-1-i, j), v flips sign.
static double asymmetry(const euler::Mesh& m, const euler::Flow& f)
{
    const double sgn[euler::NEQ] = {1.0, 1.0, -1.0, 1.0};
    double d = 0.0;
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            for (int k = 0; k < euler::NEQ; ++k) {
                d = std::max(d, std::abs(f.W[k][m.cell(i, j)] - sgn[k] * f.W[k][m.cell(m.nic - 1 - i, j)]));
            }
        }
    }
    return d;
}

// One implicit run from the free stream.
static euler::RunResult implicit_run(const euler::Mesh& m, euler::Flow& f, const euler::Freestream& fs,
                                     euler::TimeScheme ts, euler::LinearSolver ls, euler::Order o, double cfl_max,
                                     int max_iter)
{
    euler::ImplicitSettings s;
    s.scheme = ts;
    s.linsolve = ls;
    s.order = o;
    s.cfl_max = cfl_max;
    s.max_iter = max_iter;
    s.tol_drop = 1e-12;
    s.print_every = 0;
    euler::init_freestream(f, fs);
    return euler::run_implicit(m, f, fs, s);
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mesh.x>\n", argv[0]);
        return 1;
    }
    try {
        euler::Mesh m = euler::read_plot3d(argv[1]);
        euler::compute_metrics(m);
        euler::Flow f;
        euler::allocate_flow(m, f);
        std::printf("=== %s (%d x %d nodes)\n", argv[1], m.ni, m.nj);

        // Reference: explicit RK5 to machine zero.
        const euler::Freestream fs = euler::make_freestream(0.5, 1.25);
        euler::ExplicitSettings es;
        es.rk = euler::rk_scheme(euler::Order::Second, 5);
        es.max_iter = 200000;
        es.tol_drop = 1e-12;
        euler::init_freestream(f, fs);
        const euler::RunResult rk = euler::run_explicit(m, f, fs, es);
        std::printf("  [INFO] RK5 reference: %d it, %.2f s, CL %.10f CD %.10f CM %.10f\n", rk.iterations, rk.seconds,
                    rk.forces.cl, rk.forces.cd, rk.forces.cm);

        // K1: the three implicit schemes with GMRES.
        const euler::TimeScheme schemes[3] = {euler::TimeScheme::BDF1, euler::TimeScheme::BDF2, euler::TimeScheme::ROW};
        const char* names[3] = {"BDF1", "BDF2", "ROW"};
        for (int k = 0; k < 3; ++k) {
            const euler::RunResult r = implicit_run(m, f, fs, schemes[k], euler::LinearSolver::GMRES, euler::Order::Second,
                                                    1e8, 100);
            const double d = std::max({std::abs(r.forces.cl - rk.forces.cl), std::abs(r.forces.cd - rk.forces.cd),
                                       std::abs(r.forces.cm - rk.forces.cm)});
            check(r.converged && r.iterations <= 40 && d < 1e-9,
                  "K1 %s + GMRES: drop %.1e in %d it (%.2f s, RK5 %d it), max |dCL,dCD,dCM| vs RK5 = %.1e", names[k],
                  r.res / r.res0, r.iterations, r.seconds, rk.iterations, d);
        }

        // K2: symmetry.
        {
            const euler::Freestream fs0 = euler::make_freestream(0.5, 0.0);
            const euler::RunResult r = implicit_run(m, f, fs0, euler::TimeScheme::BDF1, euler::LinearSolver::GMRES,
                                                    euler::Order::Second, 1e8, 100);
            const double asym = asymmetry(m, f);
            check(r.converged && asym < 1e-10 && std::abs(r.forces.cl) < 1e-10 && std::abs(r.forces.cm) < 1e-10,
                  "K2 alpha = 0: converged in %d it, asymmetry %.1e, CL %.1e, CM %.1e", r.iterations, asym, r.forces.cl,
                  r.forces.cm);
        }

        // K3: transonic and supersonic.
        for (double M : {0.8, 1.5}) {
            const euler::Freestream fsM = euler::make_freestream(M, M < 1.0 ? 1.25 : 0.0);
            const euler::RunResult r = implicit_run(m, f, fsM, euler::TimeScheme::BDF1, euler::LinearSolver::GMRES,
                                                    euler::Order::Second, 1e8, 100);
            check(r.converged, "K3 M = %.1f: drop %.1e in %d it (%.2f s), CL %.6f CD %.6f", M, r.res / r.res0,
                  r.iterations, r.seconds, r.forces.cl, r.forces.cd);
        }

        // K4: first-order residual -> the first-order Jacobian is its exact Jacobian (up to frozen |A|).
        {
            const euler::RunResult r = implicit_run(m, f, fs, euler::TimeScheme::BDF1, euler::LinearSolver::GMRES,
                                                    euler::Order::First, 1e8, 100);
            check(r.converged && r.iterations <= 20, "K4 first order, BDF1 + GMRES: drop %.1e in %d it", r.res / r.res0,
                  r.iterations);
        }

        // K5: SGS only (slide 33), moderate CFL.
        {
            const euler::RunResult r = implicit_run(m, f, fs, euler::TimeScheme::BDF1, euler::LinearSolver::SGS,
                                                    euler::Order::Second, 50.0, 3000);
            std::printf("  [INFO] K5 BDF1 + SGS only, CFL <= 50: drop %.1e in %d it (%.2f s), CL %.8f (RK5 %.8f)\n",
                        r.res / r.res0, r.iterations, r.seconds, r.forces.cl, rk.forces.cl);
        }
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
