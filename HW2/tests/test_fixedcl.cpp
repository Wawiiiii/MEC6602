/*
   Step 9a check: fixed-CL controller (post.pdf slide 9) on a model "flow", no Euler solver.

   Model of one nonlinear iteration: the residual drops by 0.5 per iteration; an alpha change makes it
   jump in proportion to the change (a small correction of the farfield direction is a small
   perturbation of the converged flow), norm = max(norm, |dalpha| / 1 deg). CL lags behind its
   converged value CL_ss(alpha) by an error proportional to the residual, CL = CL_ss(alpha) + 0.05 norm.
   The flow is "converged" at a residual of 1e-12.

   L1 linear CL_ss: reaches CL_t = 0.4 to cl_tol, alpha -> the exact alpha*, at most 12 alpha changes
   L2 nonlinear CL_ss (lift curve bending over, like transonic): reaches the target too
   L3 start above the target: the first perturbation decreases alpha (sign of the secant step, A29)
   L4 start exactly at the target alpha: the slide's perturbation (dalpha0) is still applied once,
      then the secant brings alpha back to alpha*
   L5 max_dalpha clamps a step
   L6 controller off: never changes alpha
   L7 (with a grid) the real solver, BDF1 + GMRES, M 0.8, CL_t 0.4 from alpha 1.25: converged with
      |CL - 0.4| <= 1e-6, alpha changed; a plain fixed-alpha run at the alpha found gives the same CL
      to 1e-8; with cl_target off alpha is left unchanged

   Usage: test_fixedcl [mesh.x]   (33x33 recommended for L7; without a grid only L1-L6 run)
*/

#include "bc.h"
#include "fixedcl.h"
#include "flow.h"
#include "implicit.h"
#include "mesh.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>

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
    (ok ? g_pass : g_fail)++;
}

struct ModelRun {
    double alpha = 0.0;      // final alpha, degrees
    double cl = 0.0;         // final CL
    int iterations = 0;      // model iterations
    int changes = 0;         // alpha changes
    double largest_step = 0; // largest |alpha change|
    double first_step = 0;   // first alpha change (signed)
    bool converged = false;  // residual at 1e-12 and CL on target
};

// Run the model flow with the controller until convergence or 2000 iterations.
static ModelRun run_model(const euler::FixedClSettings& s, double alpha0, const std::function<double(double)>& cl_ss)
{
    euler::FixedClController ctl(s, alpha0);
    ModelRun r;
    double norm = 1.0; // model residual
    for (int it = 1; it <= 2000; ++it) {
        const double cl = cl_ss(ctl.alpha()) + 0.05 * norm; // CL of the partly converged flow
        const bool flow_converged = norm <= 1e-12;
        r.iterations = it;
        r.cl = cl;
        if (flow_converged && (!s.on || ctl.on_target(cl))) {
            r.converged = true;
            break;
        }
        const double a_before = ctl.alpha();
        if (ctl.update(norm, cl, flow_converged)) {
            const double step = ctl.alpha() - a_before;
            if (ctl.changes() == 1) {
                r.first_step = step;
            }
            r.largest_step = std::fmax(r.largest_step, std::abs(step));
            norm = std::fmax(norm, std::abs(step)); // the residual jumps after an alpha change (slide 10)
        } else {
            norm *= 0.5;
        }
    }
    r.alpha = ctl.alpha();
    r.changes = ctl.changes();
    return r;
}

int main(int argc, char** argv)
{
    euler::FixedClSettings s;
    s.on = true;
    s.cl_target = 0.4;

    std::printf("L1 linear lift curve CL_ss = 0.11 (alpha + 0.3)\n");
    {
        auto cl_ss = [](double a) { return 0.11 * (a + 0.3); };
        const double a_star = 0.4 / 0.11 - 0.3; // exact alpha of the target
        const ModelRun r = run_model(s, 0.0, cl_ss);
        std::printf("    alpha %.10f (exact %.10f), CL %.10f, %d iterations, %d alpha changes\n", r.alpha, a_star,
                    r.cl, r.iterations, r.changes);
        check(r.converged, "converged (residual 1e-12 and |CL - CL_t| <= 1e-6)");
        check(std::abs(r.alpha - a_star) < 1e-4, "alpha = exact alpha* to 1e-4 (error %.2e)", std::abs(r.alpha - a_star));
        check(r.changes <= 12, "at most 12 alpha changes (%d)", r.changes);
    }

    std::printf("L2 nonlinear lift curve CL_ss = 0.12 alpha - 0.004 alpha^2 + 0.05\n");
    {
        auto cl_ss = [](double a) { return 0.12 * a - 0.004 * a * a + 0.05; };
        const ModelRun r = run_model(s, 1.0, cl_ss);
        std::printf("    alpha %.10f, CL %.10f, CL_ss(alpha) %.10f, %d iterations, %d alpha changes\n", r.alpha, r.cl,
                    cl_ss(r.alpha), r.iterations, r.changes);
        check(r.converged, "converged");
        check(std::abs(cl_ss(r.alpha) - 0.4) < 1e-5, "CL_ss(alpha) = CL_t to 1e-5 (error %.2e)", std::abs(cl_ss(r.alpha) - 0.4));
        check(r.changes <= 10, "at most 10 alpha changes (%d)", r.changes);
    }

    std::printf("L3 start above the target (alpha 6, CL ~ 0.69)\n");
    {
        auto cl_ss = [](double a) { return 0.11 * (a + 0.3); };
        const ModelRun r = run_model(s, 6.0, cl_ss);
        check(r.first_step < 0.0, "first perturbation decreases alpha (%+.3f deg)", r.first_step);
        check(r.converged && std::abs(r.alpha - (0.4 / 0.11 - 0.3)) < 1e-4, "converges to alpha* (%.6f)", r.alpha);
    }

    std::printf("L4 start at the target alpha\n");
    {
        auto cl_ss = [](double a) { return 0.11 * (a + 0.3); };
        const ModelRun r = run_model(s, 0.4 / 0.11 - 0.3, cl_ss);
        check(r.converged && std::abs(r.alpha - (0.4 / 0.11 - 0.3)) < 1e-4, "back to alpha* (%.6f, %d changes)", r.alpha,
              r.changes);
    }

    std::printf("L5 step clamp max_dalpha = 0.5 deg\n");
    {
        euler::FixedClSettings c = s;
        c.max_dalpha = 0.5;
        auto cl_ss = [](double a) { return 0.11 * (a + 0.3); };
        const ModelRun r = run_model(c, -3.0, cl_ss);
        check(r.largest_step <= 0.5 + 1e-15, "largest step %.3f deg <= 0.5", r.largest_step);
        check(r.converged && std::abs(r.alpha - (0.4 / 0.11 - 0.3)) < 1e-4, "still converges (%.6f, %d changes)", r.alpha,
              r.changes);
    }

    std::printf("L6 controller off\n");
    {
        euler::FixedClSettings c = s;
        c.on = false;
        auto cl_ss = [](double a) { return 0.11 * (a + 0.3); };
        const ModelRun r = run_model(c, 1.0, cl_ss);
        check(r.changes == 0 && r.alpha == 1.0, "alpha unchanged (%.3f, %d changes)", r.alpha, r.changes);
    }

    if (argc > 1) {
        std::printf("L7 real solver %s, M 0.8, CL_t 0.4, BDF1 + GMRES\n", argv[1]);
        euler::Mesh m = euler::read_plot3d(argv[1]);
        euler::compute_metrics(m);
        euler::ImplicitSettings is; // implicit settings, defaults of the solver
        is.gmres_restart = 50;
        is.gmres_max_iter = 150;
        is.max_iter = 300;
        is.print_every = 0;
        // Run from the free stream at alpha_deg, fixed CL on or off.
        auto run = [&](double alpha_deg, bool fixed) {
            euler::ImplicitSettings s2 = is;
            s2.fixed_cl.on = fixed;
            s2.fixed_cl.cl_target = 0.4;
            euler::Flow f;
            euler::allocate_flow(m, f);
            const euler::Freestream fs = euler::make_freestream(0.8, alpha_deg);
            euler::init_freestream(f, fs);
            return euler::run_implicit(m, f, fs, s2);
        };
        const euler::RunResult rc = run(1.25, true);
        std::printf("    fixed CL: alpha %.10f deg, CL %.10f, %d iterations, %d alpha changes, drop %.1e\n", rc.alpha_deg,
                    rc.forces.cl, rc.iterations, rc.alpha_changes, rc.res / rc.res0);
        check(rc.converged && std::abs(rc.forces.cl - 0.4) <= 1e-6, "converged with |CL - 0.4| = %.1e <= 1e-6",
              std::abs(rc.forces.cl - 0.4));
        check(rc.alpha_changes > 0 && rc.alpha_deg != 1.25, "alpha moved (%d changes)", rc.alpha_changes);
        const euler::RunResult ra = run(rc.alpha_deg, false);
        std::printf("    fixed alpha %.10f: CL %.10f, %d iterations\n", rc.alpha_deg, ra.forces.cl, ra.iterations);
        check(ra.converged && std::abs(ra.forces.cl - rc.forces.cl) < 1e-8, "fixed-alpha run at that alpha: same CL (diff %.1e)",
              std::abs(ra.forces.cl - rc.forces.cl));
        check(ra.alpha_changes == 0 && ra.alpha_deg == rc.alpha_deg, "cl_target off: alpha unchanged");
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
