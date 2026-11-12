/*
   Step 5 check: explicit Euler and multistage Runge-Kutta with local time stepping
   (Euler2D_basic_algo slides 17-18, 20; Roe_scheme_accel slide 20; post.pdf slide 4).

   G1 the stage coefficients and CFL numbers are the ones of the slide 18 table
   G2 local time step on a Cartesian cell, by hand: dt = sigma Omega / ((|u| + c) dy + (|v| + c) dx)
   G3 one explicit Euler iteration = W - dt/Omega R(W), computed by hand on a random field
   G4 slide 20 sanity check over many iterations: all BCs farfield + free stream -> nothing moves
   G5 alpha = 0 on the (exactly symmetric) grid: the solution stays mirror-symmetric and the lift is 0
   G6 convergence to machine zero (drop 1e-12) at M = 0.5, alpha = 1.25: 1st and 2nd order
   G7 converged solution: total enthalpy H = H_inf (steady adiabatic Euler), entropy at the wall,
      preliminary CL / CD (pressure integration on the wall faces; the real force module is step 6)

   Usage: test_explicit <mesh.x> [<mesh.x> ...]   (G5-G7 run on each grid given)
   Exit code = number of failed checks.
*/

#include "bc.h"
#include "flow.h"
#include "mesh.h"
#include "residual.h"
#include "solver.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <limits>
#include <random>
#include <vector>

using euler::GAMMA;
using euler::NEQ;
using euler::NG;

constexpr double PI = 3.14159265358979323846; // MSVC has no M_PI by default

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

/*
   Preliminary forces from the wall-face fluxes (only to judge the solution in this step):
   the momentum flux through a wall face is p_w n dS (n into the fluid), so the force on the
   airfoil is F = -sum (p_w - p_inf) n dS. CL, CD with q_inf = rho_inf V_inf^2 / 2 = M^2 / 2, chord 1.
*/
static void prelim_forces(const euler::Mesh& m, const euler::Flow& f, const euler::Freestream& fs, euler::Order o,
                          double& cl, double& cd)
{
    double Fx = 0.0; // force on the airfoil, x
    double Fy = 0.0; // force on the airfoil, y
    for (int i = 0; i < m.nic; ++i) {
        double F[NEQ];
        euler::j_face_flux(m, f, o, i, 0, F);
        const int fc = m.node(i, 0);
        Fx -= F[1] - fs.q[3] * m.njx[fc] * m.sj[fc];
        Fy -= F[2] - fs.q[3] * m.njy[fc] * m.sj[fc];
    }
    const double a = fs.alpha_deg * PI / 180.0;  // angle of attack
    const double q = 0.5 * fs.mach * fs.mach;    // dynamic pressure
    cl = (-Fx * std::sin(a) + Fy * std::cos(a)) / q;
    cd = (Fx * std::cos(a) + Fy * std::sin(a)) / q;
}

// ===================================================================== G1 - G3

static void test_basics()
{
    std::printf("\n=== G1-G3: scheme coefficients, time step, one explicit Euler iteration\n");

    // G1: slide 18 table.
    {
        const euler::RkScheme s25 = euler::rk_scheme(euler::Order::Second, 5);
        const euler::RkScheme s13 = euler::rk_scheme(euler::Order::First, 3);
        const euler::RkScheme s24 = euler::rk_scheme(euler::Order::Second, 4);
        const euler::RkScheme e1 = euler::rk_scheme(euler::Order::Second, 1);
        const bool ok = s25.stages == 5 && s25.cfl == 1.15 && s25.alpha[0] == 0.0695 && s25.alpha[1] == 0.1602 &&
                        s25.alpha[2] == 0.2898 && s25.alpha[3] == 0.5060 && s25.alpha[4] == 1.0 &&
                        s13.stages == 3 && s13.cfl == 1.5 && s13.alpha[0] == 0.1481 && s13.alpha[1] == 0.4 &&
                        s24.cfl == 0.92 && s24.alpha[2] == 0.5052 && e1.stages == 1 && e1.alpha[0] == 1.0 &&
                        e1.cfl < 1.0;
        bool threw = false;
        try {
            euler::rk_scheme(euler::Order::Second, 2);
        } catch (const std::exception&) {
            threw = true;
        }
        check(ok && threw, "G1 slide 18 coefficients (2nd-order 5-stage CFL 1.15: 0.0695 .. 1.0), bad stage count rejected");
    }

    // G2: one Cartesian cell, unit spacing (3x3 nodes built by hand, like test_metrics).
    {
        euler::Mesh m;
        m.ni = 3;
        m.nj = 3;
        m.nic = 2;
        m.njc = 2;
        m.ni_tot = m.ni + 2 * NG;
        m.nj_tot = m.nj + 2 * NG;
        m.nic_tot = m.nic + 2 * NG;
        m.njc_tot = m.njc + 2 * NG;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        m.x.assign(static_cast<size_t>(m.ni_tot) * m.nj_tot, nan);
        m.y.assign(m.x.size(), nan);
        for (int j = 0; j < m.nj; ++j) {
            for (int i = 0; i < m.ni; ++i) {
                m.x[m.node(i, j)] = 2.0 * i; // dx = 2
                m.y[m.node(i, j)] = 0.5 * j; // dy = 0.5
            }
        }
        euler::compute_metrics(m);
        euler::Flow f;
        euler::allocate_flow(m, f);
        const double q[NEQ] = {1.3, 0.4, -0.25, 0.9}; // any state
        double w[NEQ];
        euler::prim_to_cons(q, w);
        for (int k = 0; k < NEQ; ++k) {
            std::fill(f.W[k].begin(), f.W[k].end(), w[k]);
        }
        euler::compute_primitives(f);
        std::vector<double> dt;
        euler::local_time_step(m, f, 0.7, dt);
        const double c = std::sqrt(GAMMA * q[3] / q[0]); // speed of sound
        // i-faces have length dy = 0.5 and normal (1, 0); j-faces length dx = 2, normal (0, 1); area 1.
        const double expect = 0.7 * 1.0 / ((std::abs(q[1]) + c) * 0.5 + (std::abs(q[2]) + c) * 2.0);
        double worst = 0.0;
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                worst = std::max(worst, std::abs(dt[m.cell(i, j)] - expect));
            }
        }
        check(worst < 1e-15, "G2 local time step on a 2 x 0.5 Cartesian cell: dt = %.10f (by hand %.10f)", dt[m.cell(0, 0)],
              expect);
    }
}

// G3: one explicit Euler iteration of run_explicit = the update written by hand.
static void test_one_iteration(const euler::Mesh& m, euler::Flow& f)
{
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    const euler::Freestream fs = euler::make_freestream(0.6, 2.0);
    euler::init_freestream(f, fs);
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const double q[NEQ] = {1.0 + 0.1 * U(gen), 0.6 + 0.1 * U(gen), 0.1 * U(gen), 1.0 / GAMMA + 0.05 * U(gen)};
            double w[NEQ];
            euler::prim_to_cons(q, w);
            for (int k = 0; k < NEQ; ++k) {
                f.W[k][m.cell(i, j)] = w[k];
            }
        }
    }

    // By hand.
    euler::Flow g = f;
    euler::apply_bc(m, g, fs);
    euler::compute_residual(m, g, euler::Order::Second);
    std::vector<double> dt;
    euler::local_time_step(m, g, 0.9, dt);
    std::vector<double> expect[NEQ];
    for (int k = 0; k < NEQ; ++k) {
        expect[k] = g.W[k];
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                expect[k][c] -= dt[c] / m.vol[c] * g.R[k][c];
            }
        }
    }

    // run_explicit with explicit Euler and max_iter = 1: one residual, one time step, one update.
    euler::ExplicitSettings s;
    s.order = euler::Order::Second;
    s.rk = euler::rk_scheme(euler::Order::Second, 1);
    s.rk.cfl = 0.9;
    s.tol_drop = 0.0;
    s.max_iter = 1; // one iteration = residual, time step, one update
    euler::run_explicit(m, f, fs, s);
    double worst = 0.0;
    for (int k = 0; k < NEQ; ++k) {
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                worst = std::max(worst, std::abs(f.W[k][c] - expect[k][c]));
            }
        }
    }
    check(worst < 1e-14, "G3 one explicit Euler iteration = W - dt/Omega R by hand: max diff %.2e", worst);
}

// G4: slide 20 over many iterations.
static void test_nothing_moves(const euler::Mesh& m, euler::Flow& f)
{
    const euler::Freestream fs = euler::make_freestream(0.8, 1.25);
    euler::init_freestream(f, fs);
    std::vector<double> w0[NEQ];
    for (int k = 0; k < NEQ; ++k) {
        w0[k] = f.W[k];
    }
    euler::ExplicitSettings s;
    s.order = euler::Order::Second;
    s.rk = euler::rk_scheme(euler::Order::Second, 5);
    s.jmin = euler::JminBc::Farfield;
    s.tol_drop = -1.0; // never stop on the residual (it is 0 from the start)
    s.max_iter = 200;
    euler::run_explicit(m, f, fs, s);
    double worst = 0.0;
    for (int k = 0; k < NEQ; ++k) {
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                worst = std::max(worst, std::abs(f.W[k][c] - w0[k][c]));
            }
        }
    }
    check(worst < 1e-13, "G4 all BCs farfield + free stream, 200 RK5 iterations: max |W - W_inf| = %.2e (slide 20)", worst);
}

// ===================================================================== G5 - G7 on each grid

// Largest mirror asymmetry of W: cell (i, j) vs (nic-1-i, j), v flips sign. Relative to max |W|.
static double asymmetry(const euler::Mesh& m, const euler::Flow& f)
{
    const double sgn[NEQ] = {1.0, 1.0, -1.0, 1.0};
    double d = 0.0;
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            for (int k = 0; k < NEQ; ++k) {
                d = std::max(d, std::abs(f.W[k][m.cell(i, j)] - sgn[k] * f.W[k][m.cell(m.nic - 1 - i, j)]));
            }
        }
    }
    return d;
}

static void test_grid(const euler::Mesh& m, euler::Flow& f, int max_iter)
{
    // G5: alpha = 0.
    {
        const euler::Freestream fs = euler::make_freestream(0.5, 0.0);
        euler::init_freestream(f, fs);
        euler::ExplicitSettings s;
        s.order = euler::Order::Second;
        s.rk = euler::rk_scheme(euler::Order::Second, 5);
        s.max_iter = max_iter;
        s.tol_drop = 1e-10;
        const euler::ExplicitResult r = euler::run_explicit(m, f, fs, s);
        double cl = 0.0;
        double cd = 0.0;
        euler::compute_residual(m, f, s.order);
        prelim_forces(m, f, fs, s.order, cl, cd);
        const double asym = asymmetry(m, f);
        check(r.converged && asym < 1e-11 && std::abs(cl) < 1e-10,
              "G5 M 0.5 alpha 0, 2nd order: converged %s in %d it, asymmetry %.2e, CL %.2e (expect 0), CD %.5f",
              r.converged ? "yes" : "NO", r.iterations, asym, cl, cd);
    }

    // G6 / G7: M = 0.5, alpha = 1.25, first then second order.
    const euler::Freestream fs = euler::make_freestream(0.5, 1.25);
    const double H_inf = GAMMA / (GAMMA - 1.0) * fs.q[3] / fs.q[0] + 0.5 * fs.mach * fs.mach; // free-stream total enthalpy
    const double s_inf = fs.q[3] / std::pow(fs.q[0], GAMMA);                                  // free-stream entropy p/rho^gamma
    for (euler::Order o : {euler::Order::First, euler::Order::Second}) {
        const char* name = o == euler::Order::First ? "1st" : "2nd";
        euler::init_freestream(f, fs);
        euler::ExplicitSettings s;
        s.order = o;
        s.rk = euler::rk_scheme(o, 5);
        s.max_iter = max_iter;
        s.tol_drop = 1e-12;
        const euler::ExplicitResult r = euler::run_explicit(m, f, fs, s);
        check(r.converged && !r.diverged, "G6 %s order, M 0.5 alpha 1.25, RK5 CFL %.2f: drop %.2e in %d it (%.1f s)", name,
              s.rk.cfl, r.res / r.res0, r.iterations, r.seconds);

        euler::compute_residual(m, f, o);
        double Hmax = 0.0; // max |H / H_inf - 1| over the field
        double smax = 0.0; // max |s / s_inf - 1| on the wall row
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                const double rho = f.Q[0][c];
                const double u = f.Q[1][c];
                const double v = f.Q[2][c];
                const double p = f.Q[3][c];
                const double H = GAMMA / (GAMMA - 1.0) * p / rho + 0.5 * (u * u + v * v);
                Hmax = std::max(Hmax, std::abs(H / H_inf - 1.0));
                if (j == 0) {
                    smax = std::max(smax, std::abs(p / std::pow(rho, GAMMA) / s_inf - 1.0));
                }
            }
        }
        double cl = 0.0;
        double cd = 0.0;
        prelim_forces(m, f, fs, o, cl, cd);
        std::printf("  [INFO] G7 %s order: max |H/H_inf - 1| = %.2e, max wall |s/s_inf - 1| = %.2e, "
                    "preliminary CL %.5f CD %.5f\n",
                    name, Hmax, smax, cl, cd);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mesh.x> [<mesh.x> ...]\n", argv[0]);
        return 1;
    }
    try {
        test_basics();
        for (int a = 1; a < argc; ++a) {
            euler::Mesh m = euler::read_plot3d(argv[a]);
            euler::compute_metrics(m);
            euler::Flow f;
            euler::allocate_flow(m, f);
            std::printf("\n=== %s (%d x %d nodes)\n", argv[a], m.ni, m.nj);
            test_one_iteration(m, f);
            test_nothing_moves(m, f);
            test_grid(m, f, 200000);
        }
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
