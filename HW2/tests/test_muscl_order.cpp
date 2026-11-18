/*
   Order of the MUSCL face values (include/muscl.h), 1D, no Euler solver (2026-11-18).

   Question for the report: Blazek (3rd ed., p. 91) calls kappa^ = 1/3 "a second-order upwind-biased
   scheme with lower truncation error"; a Taylor expansion of the cell averages gives face values
   exact to O(h^3) for kappa^ = 1/3 only. This test measures that directly on the reconstruction
   alone.

   Data: u(x) = sin(2 pi x) + 0.5 cos(4 pi x) on [0, 1], periodic, N uniform cells of size h = 1/N,
   cell averages computed exactly. Each face x_{I+1/2} gets the left state from cells I-1, I, I+1
   (muscl_left) and the right state from cells I, I+1, I+2 (muscl_right); both are compared with the
   exact point value u(x_{I+1/2}). The function has smooth extrema, where a limiter may clip.

   Errors: maximum and mean (L1) over all faces of |U_L - u| and |U_R - u|. Observed order between
   N and 2N: log2(e_N / e_2N). Grids N = 32 ... 2048.

   Schemes:
     S1 unlimited kappa^ = 1/3                  expected 3 (Taylor expansion)
     S2 Hemker-Koren, eps = (5 h)^3 (solver)    L1 -> 3, max ~2.5 (see below)
     S3 Hemker-Koren, eps = h^4 (slide: Omega^2) expected below 3 in the max norm (clips the extrema)
     S4 unlimited kappa^ = 0                    expected 2
     S5 Van Albada (kappa^ = 0), eps = (5 h)^3   expected 2
   eps is given the 2D meaning of the solver: limiter_eps(area) with area = h^2.

   Checks (PASS / FAIL): S1 order >= 2.8 on the two finest pairs (max and L1); S2 L1 order >= 2.7 (max
   reported only); S4 and S5 order in [1.8, 2.2] on the two finest pairs.

   S2, first run (2026-11-18): L1 2.80 / 2.82 / 2.84 (rising to 3), max 2.49. At an extremum
   a, b = O(h^2), so a^2 << eps = (5 h)^3 and the slope is the unlimited kappa^ = 1/3 one; far from it
   the limiter formula reproduces kappa^ = 1/3 to O(h^3). In between, a band of width ~ sqrt(h) where
   a^2 ~ eps has an error ~ h^2.5: it sets the max norm, its weight in the L1 norm vanishes. With the
   slide's eps = h^4 (S3) the extrema themselves are clipped: max order 2.00.
   Observed drag orders at M 0.5 alpha 0 (RESULTS.md): 3.1 limiter off, 2.6-2.8 Hemker-Koren,
   2.3 -> 2 for kappa^ = 0: the same ordering as the face values here.

   Usage: test_muscl_order
*/

#include "muscl.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <functional>
#include <vector>

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

static const double PI = 3.14159265358979323846;

// Exact point value of the test function.
static double u_exact(double x)
{
    return std::sin(2.0 * PI * x) + 0.5 * std::cos(4.0 * PI * x);
}

// Exact average of the test function over [xl, xr].
static double u_average(double xl, double xr)
{
    const double s = (std::cos(2.0 * PI * xl) - std::cos(2.0 * PI * xr)) / (2.0 * PI); // integral of sin(2 pi x)
    const double c = 0.5 * (std::sin(4.0 * PI * xr) - std::sin(4.0 * PI * xl)) / (4.0 * PI); // integral of 0.5 cos(4 pi x)
    return (s + c) / (xr - xl);
}

struct Errors {
    double max = 0.0; // largest face-value error
    double l1 = 0.0;  // mean face-value error
};

// Face-value errors of one scheme on N cells. left/right: the reconstruction (cells, eps) -> state.
using Recon = std::function<double(double, double, double, double)>;

static Errors face_errors(int N, const Recon& left, const Recon& right, const std::function<double(double)>& eps_of_h)
{
    const double h = 1.0 / N;          // cell size
    std::vector<double> ub(N);         // exact cell averages
    for (int i = 0; i < N; ++i) {
        ub[i] = u_average(i * h, (i + 1) * h);
    }
    auto at = [&](int i) { return ub[((i % N) + N) % N]; }; // periodic access
    const double eps = eps_of_h(h);    // limiter eps (same for every cell, uniform grid)
    Errors e;
    for (int i = 0; i < N; ++i) {
        const double uf = u_exact((i + 1) * h);                    // exact value at face i+1/2
        const double eL = std::abs(left(at(i - 1), at(i), at(i + 1), eps) - uf);
        const double eR = std::abs(right(at(i), at(i + 1), at(i + 2), eps) - uf);
        e.max = std::fmax(e.max, std::fmax(eL, eR));
        e.l1 += 0.5 * (eL + eR) / N;
    }
    return e;
}

struct Scheme {
    const char* name;                       // label in the output
    Recon left;                             // left state of a face
    Recon right;                            // right state of a face
    std::function<double(double)> eps;      // eps as a function of h
    double lo;                              // smallest accepted order (two finest pairs)
    double hi;                              // largest accepted order
    bool max_checked = true;                // false: only the L1 order is checked
};

int main()
{
    using namespace euler;
    const auto eps_solver = [](double h) { return limiter_eps(h * h); }; // (5 h)^3, as in the solver
    const auto eps_slide = [](double h) { return h * h * h * h; };      // Omega^2 = h^4 (slide 16)
    const auto eps_none = [](double) { return 0.0; };                   // unused by the unlimited schemes

    // Unlimited kappa scheme with a given kappa^ (kappa_slope reads the global MUSCL_KAPPA).
    auto kappa_left = [](double k) {
        return Recon([k](double um1, double u0, double up1, double) {
            MUSCL_KAPPA = k;
            return u0 + 0.5 * kappa_slope(up1 - u0, u0 - um1);
        });
    };
    auto kappa_right = [](double k) {
        return Recon([k](double u0, double up1, double up2, double) {
            MUSCL_KAPPA = k;
            return up1 - 0.5 * kappa_slope(up1 - u0, up2 - up1);
        });
    };
    const Recon hk_left = [](double um1, double u0, double up1, double eps) { return u0 + 0.5 * hk_slope(up1 - u0, u0 - um1, eps); };
    const Recon hk_right = [](double u0, double up1, double up2, double eps) { return up1 - 0.5 * hk_slope(up1 - u0, up2 - up1, eps); };
    const Recon va_left = [](double um1, double u0, double up1, double eps) { return u0 + 0.5 * va_slope(up1 - u0, u0 - um1, eps); };
    const Recon va_right = [](double u0, double up1, double up2, double eps) { return up1 - 0.5 * va_slope(up1 - u0, up2 - up1, eps); };

    const std::vector<Scheme> schemes = {
        {"S1 unlimited kappa^ = 1/3", kappa_left(1.0 / 3.0), kappa_right(1.0 / 3.0), eps_none, 2.8, 99.0},
        {"S2 Hemker-Koren, eps = (5h)^3", hk_left, hk_right, eps_solver, 2.7, 99.0, false}, // L1 only (max ~2.5, header)
        {"S3 Hemker-Koren, eps = h^4 (slide)", hk_left, hk_right, eps_slide, -99.0, 99.0}, // reported only
        {"S4 unlimited kappa^ = 0", kappa_left(0.0), kappa_right(0.0), eps_none, 1.8, 2.2},
        {"S5 Van Albada, eps = (5h)^3", va_left, va_right, eps_solver, 1.8, 2.2},
    };
    const std::vector<int> grids = {32, 64, 128, 256, 512, 1024, 2048}; // number of cells N

    for (const Scheme& s : schemes) {
        std::printf("%s\n", s.name);
        std::printf("  %6s %12s %6s %12s %6s\n", "N", "max error", "order", "L1 error", "order");
        std::vector<Errors> e;
        for (int N : grids) {
            e.push_back(face_errors(N, s.left, s.right, s.eps));
        }
        std::vector<double> pmax(grids.size(), 0.0); // order of the max error, N/2 -> N
        std::vector<double> pl1(grids.size(), 0.0);  // order of the L1 error
        for (size_t k = 0; k < grids.size(); ++k) {
            if (k > 0) {
                pmax[k] = std::log2(e[k - 1].max / e[k].max);
                pl1[k] = std::log2(e[k - 1].l1 / e[k].l1);
                std::printf("  %6d %12.4e %6.2f %12.4e %6.2f\n", grids[k], e[k].max, pmax[k], e[k].l1, pl1[k]);
            } else {
                std::printf("  %6d %12.4e %6s %12.4e %6s\n", grids[k], e[k].max, "", e[k].l1, "");
            }
        }
        const size_t n = grids.size();
        if (s.lo > -50.0) {
            auto in = [&](double p) { return p >= s.lo && p <= s.hi; }; // order inside the accepted range
            const bool ok_max = !s.max_checked || (in(pmax[n - 1]) && in(pmax[n - 2]));
            const bool ok = ok_max && in(pl1[n - 1]) && in(pl1[n - 2]);
            check(ok, "order on the two finest pairs in [%.1f, %.1f]: max %.2f %.2f%s, L1 %.2f %.2f", s.lo,
                  s.hi > 50.0 ? 9.9 : s.hi, pmax[n - 2], pmax[n - 1], s.max_checked ? "" : " (not checked)",
                  pl1[n - 2], pl1[n - 1]);
        }
    }
    MUSCL_KAPPA = 1.0 / 3.0; // restore the default

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
