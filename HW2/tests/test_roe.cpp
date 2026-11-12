/*
   Step 3 check: Roe flux + MUSCL / Hemker-Koren + residual (Roe_scheme_accel slides 5-16).

   A. Roe flux alone (no mesh)
      A1 consistency          F_Roe(W, W) = Fc(W)                       (slide 5, property 2)
      A2 conservation         F(L, R, n) = -F(R, L, -n)                 (what leaves L enters R)
      A3 upwinding            all waves going one way -> F = Fc(L) or Fc(R)  (slide 6; tests property 1)
      A4 rotation             the flux turns with the frame
      A5 stationary shock     Rankine-Hugoniot states: V~ - c~ = 0 and F = Fc(L) without the fix (slide 10)
      A6 stationary slip line exact flux (0, p nx, p ny, 0) without the fix
      A7 Harten correction    values and continuity (slide 11)
   B. MUSCL + Hemker-Koren alone (1D data)
      B1 linear data reproduced exactly at the face
      B2 left/right mirror symmetry
      B3 smooth form of slide 16 = psi(r) form of slide 15, both states
      B4 a step stays a step (no new values next to a discontinuity)
      B5 monotone data: no face value outside its two neighbour cells (slide 14); extremum overshoot printed
      B6 order of accuracy of the face values on smooth data (expect 2)
   C. Residual on the airfoil grids
      C1 uniform flow -> residual 0, first and second order          (slide 12 "to verify")
      C2 i-cut: face i = 0 and face i = ni-1 give the same flux       (random field)
      C3 conservation: sum of the cell residuals = net boundary flux  (random field)
      C4 order of the residual on a smooth field vs. the exact divergence (needs >= 3 nested grids)
         -> L1 order checked (1 and 2); the second-order L2 order is printed, it is lower (see mms_errors)
      timing of one second-order residual on the finest grid given

   Usage: test_roe <mesh.x> [<mesh.x> ...]   nested grids, coarse to fine.
   Exit code = number of failed checks.
*/

#include "flow.h"
#include "mesh.h"
#include "muscl.h"
#include "residual.h"
#include "roe.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <random>
#include <string>
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

// Largest |a_k - b_k| over the 4 components.
static double max_abs_diff(const double a[NEQ], const double b[NEQ])
{
    double d = 0.0;
    for (int k = 0; k < NEQ; ++k) {
        d = std::max(d, std::abs(a[k] - b[k]));
    }
    return d;
}

// Size of a flux vector, used to turn absolute differences into relative ones (never below 1).
static double flux_scale(const double F[NEQ])
{
    double s = 1.0;
    for (int k = 0; k < NEQ; ++k) {
        s = std::max(s, std::abs(F[k]));
    }
    return s;
}

// Random primitive state: rho in [0.3, 2], p in [0.2, 2], velocity of size up to mach_max * c.
static void random_state(std::mt19937& gen, double q[NEQ], double mach_max)
{
    std::uniform_real_distribution<double> U(0.0, 1.0);
    q[0] = 0.3 + 1.7 * U(gen);
    q[3] = 0.2 + 1.8 * U(gen);
    const double c = std::sqrt(GAMMA * q[3] / q[0]);   // speed of sound of this state
    const double speed = mach_max * c * U(gen);       // velocity magnitude
    const double angle = 2.0 * PI * U(gen);           // velocity direction
    q[1] = speed * std::cos(angle);
    q[2] = speed * std::sin(angle);
}

// Random unit normal.
static void random_normal(std::mt19937& gen, double& nx, double& ny)
{
    std::uniform_real_distribution<double> U(0.0, 2.0 * PI);
    const double a = U(gen); // normal direction
    nx = std::cos(a);
    ny = std::sin(a);
}

/*
   Roe average written again here (slide 9), independently of roe.cpp, so the tests can
   check where the waves go: returns V~ - c~, V~ + c~ and c~.
*/
static void roe_waves(const double qL[NEQ], const double qR[NEQ], double nx, double ny,
                      double& lam1, double& lam5, double& c_t)
{
    const double HL = GAMMA / (GAMMA - 1.0) * qL[3] / qL[0] + 0.5 * (qL[1] * qL[1] + qL[2] * qL[2]); // left H
    const double HR = GAMMA / (GAMMA - 1.0) * qR[3] / qR[0] + 0.5 * (qR[1] * qR[1] + qR[2] * qR[2]); // right H
    const double sL = std::sqrt(qL[0]); // sqrt(rho_L)
    const double sR = std::sqrt(qR[0]); // sqrt(rho_R)
    const double u = (sL * qL[1] + sR * qR[1]) / (sL + sR); // u~
    const double v = (sL * qL[2] + sR * qR[2]) / (sL + sR); // v~
    const double H = (sL * HL + sR * HR) / (sL + sR);       // H~
    c_t = std::sqrt((GAMMA - 1.0) * (H - 0.5 * (u * u + v * v)));
    const double V = u * nx + v * ny; // V~
    lam1 = V - c_t;
    lam5 = V + c_t;
}

// ===================================================================== A. Roe flux

static void test_roe_flux()
{
    std::printf("\n=== A. Roe flux (slides 5-12)\n");
    std::mt19937 gen(12345);

    // A1: no jump -> no dissipation -> the physical flux.
    {
        double worst = 0.0; // largest relative difference
        for (int n = 0; n < 10000; ++n) {
            double q[NEQ];
            double nx, ny;
            random_state(gen, q, 3.0);
            random_normal(gen, nx, ny);
            double Fr[NEQ];
            double Fp[NEQ];
            euler::roe_flux(q, q, nx, ny, Fr);
            euler::physical_flux(q, nx, ny, Fp);
            worst = std::max(worst, max_abs_diff(Fr, Fp) / flux_scale(Fp));
        }
        check(worst < 1e-14, "A1 consistency F(W,W) = Fc(W), 10000 random states: max rel. diff %.2e", worst);
    }

    // A2: swapping the sides and flipping the normal flips the flux.
    {
        double worst = 0.0; // largest relative difference
        for (int n = 0; n < 10000; ++n) {
            double qL[NEQ];
            double qR[NEQ];
            double nx, ny;
            random_state(gen, qL, 3.0);
            random_state(gen, qR, 3.0);
            random_normal(gen, nx, ny);
            double F1[NEQ];
            double F2[NEQ];
            euler::roe_flux(qL, qR, nx, ny, F1);
            euler::roe_flux(qR, qL, -nx, -ny, F2);
            for (int k = 0; k < NEQ; ++k) {
                F2[k] = -F2[k];
            }
            worst = std::max(worst, max_abs_diff(F1, F2) / flux_scale(F1));
        }
        check(worst < 1e-13, "A2 conservation F(L,R,n) = -F(R,L,-n), 10000 random pairs: max rel. diff %.2e", worst);
    }

    /*
       A3: when all Roe eigenvalues are positive (and above delta, so the fix is inactive),
       |A_Roe| = A_Roe and F = 1/2 [FL + FR - (FR - FL)] = FL, but only if the Roe matrix
       satisfies property 1: FR - FL = A_Roe (WR - WL). So this checks the whole |A| dW
       assembly against property 1, on pairs with large jumps. Same thing backwards -> FR.
    */
    {
        int tested = 0;     // number of pairs that met the condition
        double worst = 0.0; // largest relative difference
        for (int n = 0; n < 20000; ++n) {
            double qL[NEQ];
            double qR[NEQ];
            double nx, ny;
            random_state(gen, qL, 4.0);
            random_state(gen, qR, 4.0);
            random_normal(gen, nx, ny);
            double lam1, lam5, c_t;
            roe_waves(qL, qR, nx, ny, lam1, lam5, c_t);
            const double delta = euler::HARTEN_FRACTION * c_t; // entropy fix width
            double F[NEQ];
            double Fexp[NEQ];
            if (lam1 > delta) { // all waves go from L to R
                euler::roe_flux(qL, qR, nx, ny, F);
                euler::physical_flux(qL, nx, ny, Fexp);
            } else if (lam5 < -delta) { // all waves go from R to L
                euler::roe_flux(qL, qR, nx, ny, F);
                euler::physical_flux(qR, nx, ny, Fexp);
            } else {
                continue;
            }
            ++tested;
            worst = std::max(worst, max_abs_diff(F, Fexp) / flux_scale(Fexp));
        }
        check(tested > 100 && worst < 1e-12,
              "A3 supersonic upwinding F = Fc(upwind side), %d random pairs: max rel. diff %.2e", tested, worst);
    }

    // A4: rotating the two states and the normal by theta must rotate the momentum flux by theta.
    {
        double worst = 0.0; // largest relative difference
        for (int n = 0; n < 10000; ++n) {
            double qL[NEQ];
            double qR[NEQ];
            double nx, ny;
            random_state(gen, qL, 3.0);
            random_state(gen, qR, 3.0);
            random_normal(gen, nx, ny);
            const double th = 2.0 * PI * std::uniform_real_distribution<double>(0.0, 1.0)(gen); // rotation angle
            const double ct = std::cos(th);
            const double st = std::sin(th);
            double qLr[NEQ] = {qL[0], ct * qL[1] - st * qL[2], st * qL[1] + ct * qL[2], qL[3]}; // rotated left
            double qRr[NEQ] = {qR[0], ct * qR[1] - st * qR[2], st * qR[1] + ct * qR[2], qR[3]}; // rotated right
            double F[NEQ];
            double Fr[NEQ];
            euler::roe_flux(qL, qR, nx, ny, F);
            euler::roe_flux(qLr, qRr, ct * nx - st * ny, st * nx + ct * ny, Fr);
            const double Fexp[NEQ] = {F[0], ct * F[1] - st * F[2], st * F[1] + ct * F[2], F[3]}; // F turned by theta
            worst = std::max(worst, max_abs_diff(Fr, Fexp) / flux_scale(F));
        }
        check(worst < 1e-13, "A4 rotation invariance, 10000 random cases: max rel. diff %.2e", worst);
    }

    /*
       A5: stationary normal shock, M1 = 2 (Rankine-Hugoniot). Slide 10: the Roe average makes
       the shock a single wave with zero speed (V~ - c~ = 0), so |A| dW = 0 and F = Fc(L) = Fc(R).
       The Harten fix gives that wave |lambda| = delta/2 > 0, so with the fix the flux differs
       slightly: printed for information.
    */
    {
        const double M1 = 2.0;                   // upstream Mach number
        const double rho1 = 1.0;                 // upstream density
        const double p1 = 1.0 / GAMMA;           // upstream pressure (c1 = 1)
        const double u1 = M1;                    // upstream velocity
        const double rho2 = rho1 * (GAMMA + 1.0) * M1 * M1 / ((GAMMA - 1.0) * M1 * M1 + 2.0); // RH density
        const double p2 = p1 * (1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (M1 * M1 - 1.0));        // RH pressure
        const double u2 = u1 * rho1 / rho2;      // mass conservation

        // Oblique frame: shock normal along (nx, ny), plus the same tangential velocity on both sides.
        const double nx = 0.6;  // shock normal x
        const double ny = 0.8;  // shock normal y
        const double vt = 0.37; // tangential velocity (continuous across a shock)
        const double qL[NEQ] = {rho1, u1 * nx - vt * ny, u1 * ny + vt * nx, p1}; // upstream state
        const double qR[NEQ] = {rho2, u2 * nx - vt * ny, u2 * ny + vt * nx, p2}; // downstream state

        double FL[NEQ];
        double FR[NEQ];
        euler::physical_flux(qL, nx, ny, FL);
        euler::physical_flux(qR, nx, ny, FR);
        check(max_abs_diff(FL, FR) < 1e-14, "A5 RH states give Fc(L) = Fc(R) (sanity of the test): %.2e",
              max_abs_diff(FL, FR));

        double lam1, lam5, c_t;
        roe_waves(qL, qR, nx, ny, lam1, lam5, c_t);
        check(std::abs(lam1) < 1e-14, "A5 Roe-averaged shock speed V~ - c~ = %.2e (expect 0)", lam1);

        double F0[NEQ];
        double F1[NEQ];
        euler::roe_flux(qL, qR, nx, ny, F0, 0.0);
        euler::roe_flux(qL, qR, nx, ny, F1);
        check(max_abs_diff(F0, FL) < 1e-13, "A5 without the fix the shock is exact: |F - Fc(L)| = %.2e",
              max_abs_diff(F0, FL));
        std::printf("         info: with the fix (delta = c~/10): |F - Fc(L)| = %.3e\n", max_abs_diff(F1, FL));
    }

    /*
       A6: stationary slip line (contact + shear): V = 0 on both sides, same p, different rho
       and tangential velocity. Exact flux = (0, p nx, p ny, 0). Only the V~ waves carry the jump,
       |V~| = 0, so Roe is exact without the fix. Slide 11 also applies the fix to |V~|, so with
       the fix the slip line gets some dissipation: printed for information.
    */
    {
        const double nx = -0.28; // face normal x
        const double ny = 0.96;  // face normal y
        const double p = 0.8;    // common pressure
        const double qL[NEQ] = {1.0, 0.5 * -ny, 0.5 * nx, p};  // tangential speed +0.5
        const double qR[NEQ] = {0.4, -0.3 * -ny, -0.3 * nx, p}; // tangential speed -0.3
        const double Fexp[NEQ] = {0.0, p * nx, p * ny, 0.0};   // exact flux
        double F0[NEQ];
        double F1[NEQ];
        euler::roe_flux(qL, qR, nx, ny, F0, 0.0);
        euler::roe_flux(qL, qR, nx, ny, F1);
        check(max_abs_diff(F0, Fexp) < 1e-15, "A6 slip line exact without the fix: |F - exact| = %.2e",
              max_abs_diff(F0, Fexp));
        std::printf("         info: with the fix: |F - exact| = %.3e\n", max_abs_diff(F1, Fexp));
    }

    /*
       A7: Harten correction (slide 11). Just below delta, h(lambda) - lambda = (delta - lambda)^2 / (2 delta),
       so the fix joins |lambda| with the same value AND the same slope at |lambda| = delta.
    */
    {
        const double d = 0.1;                // delta
        const double lam = d * (1.0 - 1e-6); // just below delta
        const bool ok = std::abs(euler::harten(0.0, d) - 0.5 * d) < 1e-16                           // |lambda| = 0 -> delta/2
                     && std::abs(euler::harten(d, d) - d) < 1e-16                                   // continuous at delta
                     && std::abs(euler::harten(lam, d) - lam - (d - lam) * (d - lam) / (2.0 * d)) < 1e-17 // slope 1 at delta
                     && euler::harten(0.3, d) == 0.3                                                // untouched above delta
                     && euler::harten(0.0, 0.0) == 0.0;                                             // fix off, no 0/0
        check(ok, "A7 Harten: h(0)=%.3f (d/2=%.3f), h(d)=%.3f, h(0.3)=%.3f, h(0) with d=0 -> %.1f",
              euler::harten(0.0, d), 0.5 * d, euler::harten(d, d), euler::harten(0.3, d), euler::harten(0.0, 0.0));
    }
}

// ===================================================================== B. MUSCL

// Slide 15 form: Psi(r) = 1/2 [(1 + k) r + (1 - k)] psi(r), k = 1/3, psi(r) = 3r / (2r^2 - r + 2).
static double Psi_hk(double r)
{
    const double k = 1.0 / 3.0; // kappa^
    return 0.5 * ((1.0 + k) * r + (1.0 - k)) * 3.0 * r / (2.0 * r * r - r + 2.0);
}

static void test_muscl()
{
    std::printf("\n=== B. MUSCL + Hemker-Koren (slides 13-16)\n");
    std::mt19937 gen(777);
    std::uniform_real_distribution<double> U(-1.0, 1.0);

    // B1: on linear data U_I = 3 + 2 I, both states at face I+1/2 must be U = 3 + 2 (I + 1/2) = 4 (I = 0).
    {
        const double UL = euler::muscl_left(1.0, 3.0, 5.0, 0.0);  // cells -1, 0, 1
        const double UR = euler::muscl_right(3.0, 5.0, 7.0, 0.0); // cells 0, 1, 2
        check(std::abs(UL - 4.0) < 1e-15 && std::abs(UR - 4.0) < 1e-15,
              "B1 linear data: UL = %.15f, UR = %.15f (expect 4)", UL, UR);
    }

    /*
       B2: mirror the data (I -> -I). The right state of the original face must equal the
       left state of the mirrored face: muscl_right(A, B, C) = muscl_left(C, B, A).
       This is the check that tells the a/b order of the right state apart.
    */
    {
        double worst = 0.0; // largest difference
        for (int n = 0; n < 10000; ++n) {
            const double A = U(gen);
            const double B = U(gen);
            const double C = U(gen);
            const double eps = 1e-6 * (U(gen) + 1.0);
            worst = std::max(worst, std::abs(euler::muscl_right(A, B, C, eps) - euler::muscl_left(C, B, A, eps)));
        }
        check(worst < 1e-15, "B2 left/right mirror symmetry, 10000 random triples: max diff %.2e", worst);
    }

    // B3: slide 16 smooth form (eps = 0) against the slide 15 psi(r) form, for both states.
    {
        double worstL = 0.0; // left state, largest difference
        double worstR = 0.0; // right state, largest difference
        for (int n = 0; n < 10000; ++n) {
            const double Um1 = U(gen); // U_{I-1}
            const double U0 = U(gen);  // U_I
            const double Up1 = U(gen); // U_{I+1}
            const double Up2 = U(gen); // U_{I+2}
            const double rL = (Up1 - U0) / (U0 - Um1);  // r_L (slide 15)
            const double rR = (Up1 - U0) / (Up2 - Up1); // r_R (slide 15)
            const double UL15 = U0 + 0.5 * Psi_hk(rL) * (U0 - Um1);   // UL, slide 15
            const double UR15 = Up1 - 0.5 * Psi_hk(rR) * (Up2 - Up1); // UR, slide 15
            worstL = std::max(worstL, std::abs(euler::muscl_left(Um1, U0, Up1, 0.0) - UL15));
            worstR = std::max(worstR, std::abs(euler::muscl_right(U0, Up1, Up2, 0.0) - UR15));
        }
        check(worstL < 1e-12 && worstR < 1e-12,
              "B3 slide 16 form = slide 15 psi(r) form: left %.2e, right %.2e", worstL, worstR);

        // Same check with the slide 16 labels taken literally for the right state (a_R = D+U_{I+1}, b_R = D-U_{I+1}).
        double worstLit = 0.0; // largest difference with the literal labels
        for (int n = 0; n < 1000; ++n) {
            const double U0 = U(gen);
            const double Up1 = U(gen);
            const double Up2 = U(gen);
            const double rR = (Up1 - U0) / (Up2 - Up1);
            const double UR15 = Up1 - 0.5 * Psi_hk(rR) * (Up2 - Up1);
            const double URlit = Up1 - 0.5 * euler::hk_slope(Up2 - Up1, Up1 - U0, 0.0);
            worstLit = std::max(worstLit, std::abs(URlit - UR15));
        }
        std::printf("         info: slide 16 a_R/b_R labels taken literally differ from slide 15 by up to %.2e\n",
                    worstLit);
    }

    // B4: a step 0 0 | 1 1 stays a step: the flat side has a zero difference -> zero slope.
    {
        const double UL = euler::muscl_left(0.0, 0.0, 1.0, 1e-12);
        const double UR = euler::muscl_right(0.0, 1.0, 1.0, 1e-12);
        check(std::abs(UL) < 1e-11 && std::abs(UR - 1.0) < 1e-11, "B4 step 0 0 | 1 1: UL = %.2e, UR = %.12f", UL, UR);
    }

    /*
       B5: monotonicity (slide 14). On monotone data (a and b of the same sign) the face value
       must lie between the two cells next to the face: no new extremum. At a local extremum
       (a b < 0) the smooth slide 16 form is not switched off (see muscl.h), so a small overshoot
       is possible there: its largest size is printed for information.
    */
    {
        double mono = 0.0;    // largest excursion on monotone stencils (must be <= 0)
        double extrem = 0.0;  // largest excursion at a local extremum, relative to the local jump
        int n_mono = 0;       // monotone stencils seen
        for (int n = 0; n < 100000; ++n) {
            const double Um1 = U(gen);
            const double U0 = U(gen);
            const double Up1 = U(gen);
            const double Up2 = U(gen);
            const double lo = std::min(U0, Up1); // the two cells next to face I+1/2
            const double hi = std::max(U0, Up1);
            const double UL = euler::muscl_left(Um1, U0, Up1, 1e-12);
            const double UR = euler::muscl_right(U0, Up1, Up2, 1e-12);
            const double exL = std::max(lo - UL, UL - hi); // excursion of the left state
            const double exR = std::max(lo - UR, UR - hi); // excursion of the right state
            const bool monoL = (Up1 - U0) * (U0 - Um1) > 0.0;  // left stencil monotone
            const bool monoR = (Up1 - U0) * (Up2 - Up1) > 0.0; // right stencil monotone
            if (monoL) {
                mono = std::max(mono, exL);
                ++n_mono;
            } else {
                extrem = std::max(extrem, exL / std::max(std::abs(Up1 - U0), std::abs(U0 - Um1)));
            }
            if (monoR) {
                mono = std::max(mono, exR);
                ++n_mono;
            } else {
                extrem = std::max(extrem, exR / std::max(std::abs(Up1 - U0), std::abs(Up2 - Up1)));
            }
        }
        check(mono <= 1e-15, "B5 monotone data, %d stencils: face value never outside its two cells (max excursion %.2e)",
              n_mono, mono);
        std::printf("  [INFO] B5 at a local extremum: overshoot up to %.3f x the largest local jump; (0.5, 1, 0) -> UL = %.4f\n",
                    extrem, euler::muscl_left(0.5, 1.0, 0.0, 1e-12));
    }

    /*
       B6: order of the face values. f(x) = sin(2 pi x) + 0.5 cos(6 pi x) at x_I = (I + 1/2) h,
       face values compared with f((I + 1) h). eps = h^4 (a 2D cell area squared ~ h^4).
       Periodic grid, N = 32 .. 1024, mean (L1) and max errors over all faces.
       The max error sits at an extremum of f, where the limiter acts and where its size depends
       on where the points fall around the extremum: its pair-to-pair order jumps around 2.
       So the check is on the L1 order between each pair and on the max-norm order over the
       whole range (first to last grid).
    */
    {
        auto f = [](double x) { return std::sin(2.0 * PI * x) + 0.5 * std::cos(6.0 * PI * x); };
        double prev1 = 0.0;   // L1 error on the previous grid
        double first_inf = 0.0; // max error on the first grid
        double last_inf = 0.0;  // max error on the last grid
        double pmin1 = 1e300; // smallest L1 order between two grids
        std::printf("         N     L1 error  order     max error  order   (left and right states together)\n");
        for (int N = 32; N <= 1024; N *= 2) {
            const double h = 1.0 / N;
            const double eps = h * h * h * h;
            double e1 = 0.0;   // L1 error (mean over the 2N face states)
            double einf = 0.0; // max error
            for (int I = 0; I < N; ++I) {
                const double x = (I + 0.5) * h; // centre of cell I
                const double exact = f(x + 0.5 * h);
                const double eL = std::abs(euler::muscl_left(f(x - h), f(x), f(x + h), eps) - exact);
                const double eR = std::abs(euler::muscl_right(f(x), f(x + h), f(x + 2.0 * h), eps) - exact);
                e1 += (eL + eR) / (2.0 * N);
                einf = std::max(einf, std::max(eL, eR));
            }
            if (N == 32) {
                std::printf("      %4d   %.3e            %.3e\n", N, e1, einf);
                first_inf = einf;
            } else {
                const double p1 = std::log2(prev1 / e1);
                pmin1 = std::min(pmin1, p1);
                std::printf("      %4d   %.3e  %5.2f     %.3e  %5.2f\n", N, e1, p1, einf, std::log2(last_inf / einf));
            }
            prev1 = e1;
            last_inf = einf;
        }
        const double pinf = std::log2(first_inf / last_inf) / 5.0; // max-norm order over 5 refinements
        check(pmin1 > 1.9 && pinf > 1.9, "B6 face values second order: lowest L1 order %.2f, max-norm order 32->1024 %.2f",
              pmin1, pinf);
    }
}

// ===================================================================== C. residual on the grids

// Fill every stored cell (ghosts included) with one primitive state.
static void fill_uniform(euler::Flow& f, const double q[NEQ])
{
    double w[NEQ];
    euler::prim_to_cons(q, w);
    for (int k = 0; k < NEQ; ++k) {
        std::fill(f.W[k].begin(), f.W[k].end(), w[k]);
    }
}

// Copy the physical cells across the i-cut into the i ghost cells (what step 4 will do).
static void fill_icut_ghosts(const euler::Mesh& m, euler::Flow& f)
{
    for (int k = 0; k < NEQ; ++k) {
        for (int j = 0; j < m.njc; ++j) {
            for (int g = 1; g <= NG; ++g) {
                f.W[k][m.cell(-g, j)] = f.W[k][m.cell(m.nic - g, j)];
                f.W[k][m.cell(m.nic - 1 + g, j)] = f.W[k][m.cell(g - 1, j)];
            }
        }
    }
}

// Largest |R| over the physical cells and the 4 equations, divided by the cell perimeter.
static double max_residual_per_perimeter(const euler::Mesh& m, const euler::Flow& f, bool& has_nan)
{
    double worst = 0.0;
    has_nan = false;
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const double per = m.si[m.node(i, j)] + m.si[m.node(i + 1, j)] + m.sj[m.node(i, j)] + m.sj[m.node(i, j + 1)];
            for (int k = 0; k < NEQ; ++k) {
                const double r = f.R[k][m.cell(i, j)];
                if (std::isnan(r)) {
                    has_nan = true;
                }
                worst = std::max(worst, std::abs(r) / per);
            }
        }
    }
    return worst;
}

// C1: free stream everywhere -> zero residual (slide 12), for a few Mach numbers and angles.
static void test_uniform(const euler::Mesh& m, euler::Flow& f)
{
    const double cases[3][2] = {{0.5, 1.25}, {0.8, -3.0}, {1.5, 10.0}}; // (Mach, alpha in degrees)
    double worst = 0.0;
    bool any_nan = false;
    for (const auto& c : cases) {
        const double a = c[1] * PI / 180.0; // angle of attack
        const double q[NEQ] = {1.0, c[0] * std::cos(a), c[0] * std::sin(a), 1.0 / GAMMA}; // free stream, c = 1
        fill_uniform(f, q);
        for (euler::Order order : {euler::Order::First, euler::Order::Second}) {
            euler::compute_residual(m, f, order);
            bool has_nan = false;
            worst = std::max(worst, max_residual_per_perimeter(m, f, has_nan));
            any_nan = any_nan || has_nan;
        }
    }
    check(!any_nan && worst < 1e-13,
          "C1 uniform flow (M 0.5/0.8/1.5, 1st+2nd order): max |R|/perimeter = %.2e%s", worst,
          any_nan ? "  NaN found!" : "");
}

// Random but positive field in the physical cells; i ghosts across the cut; j ghosts = copies of the edge rows.
static void fill_random(const euler::Mesh& m, euler::Flow& f, std::mt19937& gen)
{
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const double q[NEQ] = {1.0 + 0.3 * U(gen), 0.4 * U(gen), 0.4 * U(gen), 0.7 + 0.2 * U(gen)};
            double w[NEQ];
            euler::prim_to_cons(q, w);
            for (int k = 0; k < NEQ; ++k) {
                f.W[k][m.cell(i, j)] = w[k];
            }
        }
    }
    fill_icut_ghosts(m, f);
    for (int k = 0; k < NEQ; ++k) {
        for (int i = 0; i < m.nic; ++i) {
            for (int g = 1; g <= NG; ++g) {
                f.W[k][m.cell(i, -g)] = f.W[k][m.cell(i, 0)];
                f.W[k][m.cell(i, m.njc - 1 + g)] = f.W[k][m.cell(i, m.njc - 1)];
            }
        }
    }
}

/*
   C2: face i = 0 and face i = ni-1 are the same physical face (the cut). With the ghosts
   copied across the cut, both must give the same flux.
   C3: interior fluxes cancel in pairs, so the sum of the physical residuals equals the net
   flux through the boundary faces (wall j = 0, farfield j = nj-1, the two sides of the cut).
*/
static void test_random_field(const euler::Mesh& m, euler::Flow& f)
{
    std::mt19937 gen(4242);
    fill_random(m, f, gen);

    for (euler::Order order : {euler::Order::First, euler::Order::Second}) {
        const char* name = order == euler::Order::First ? "1st" : "2nd";
        euler::compute_residual(m, f, order); // also fills f.Q

        double cut = 0.0; // largest flux difference across the cut
        for (int j = 0; j < m.njc; ++j) {
            double Fa[NEQ];
            double Fb[NEQ];
            euler::i_face_flux(m, f, order, 0, j, Fa);
            euler::i_face_flux(m, f, order, m.ni - 1, j, Fb);
            cut = std::max(cut, max_abs_diff(Fa, Fb) / flux_scale(Fa));
        }
        check(cut < 1e-12, "C2 %s order: flux at i = 0 equals flux at i = ni-1 (the cut): max rel. diff %.2e", name, cut);

        double worst = 0.0; // largest relative conservation error over the 4 equations
        for (int k = 0; k < NEQ; ++k) {
            double sumR = 0.0;   // sum of the physical residuals
            double sumB = 0.0;   // net boundary flux
            double scale = 0.0;  // sum of |boundary contributions|, for a relative error
            double F[NEQ];
            for (int j = 0; j < m.njc; ++j) {
                for (int i = 0; i < m.nic; ++i) {
                    sumR += f.R[k][m.cell(i, j)];
                }
            }
            for (int i = 0; i < m.nic; ++i) {
                euler::j_face_flux(m, f, order, i, 0, F);      // wall: physical cell on the right -> minus
                sumB -= F[k];
                scale += std::abs(F[k]);
                euler::j_face_flux(m, f, order, i, m.nj - 1, F); // farfield: physical cell on the left -> plus
                sumB += F[k];
                scale += std::abs(F[k]);
            }
            for (int j = 0; j < m.njc; ++j) {
                euler::i_face_flux(m, f, order, 0, j, F);      // cut, physical cell on the right
                sumB -= F[k];
                scale += std::abs(F[k]);
                euler::i_face_flux(m, f, order, m.ni - 1, j, F); // cut, physical cell on the left
                sumB += F[k];
                scale += std::abs(F[k]);
            }
            worst = std::max(worst, std::abs(sumR - sumB) / scale);
        }
        check(worst < 1e-12, "C3 %s order: sum of residuals = net boundary flux: max rel. diff %.2e", name, worst);
    }
}

// Smooth manufactured field (not a solution of Euler, just smooth): primitive variables at (x, y).
static void mms_prim(double x, double y, double q[NEQ])
{
    q[0] = 1.0 + 0.2 * std::sin(1.3 * x + 0.4) * std::cos(0.9 * y - 0.2);
    q[1] = 0.5 + 0.1 * std::cos(0.8 * x) * std::sin(1.1 * y + 0.3);
    q[2] = 0.2 + 0.1 * std::sin(0.7 * x - 0.5) * std::cos(1.2 * y);
    q[3] = 1.0 / GAMMA + 0.1 * std::cos(1.1 * x + 0.2) * std::cos(0.6 * y + 0.1);
}

// Exact divergence dF/dx + dG/dy of the manufactured field (central differences of the exact flux, error ~1e-10).
static void mms_div(double x, double y, double d[NEQ])
{
    const double h = 1e-5; // differencing step
    double qa[NEQ];
    double qb[NEQ];
    double Fa[NEQ];
    double Fb[NEQ];
    double Ga[NEQ];
    double Gb[NEQ];
    mms_prim(x + h, y, qa);
    mms_prim(x - h, y, qb);
    euler::physical_flux(qa, 1.0, 0.0, Fa); // F = flux through a face of normal (1, 0)
    euler::physical_flux(qb, 1.0, 0.0, Fb);
    mms_prim(x, y + h, qa);
    mms_prim(x, y - h, qb);
    euler::physical_flux(qa, 0.0, 1.0, Ga); // G = flux through a face of normal (0, 1)
    euler::physical_flux(qb, 0.0, 1.0, Gb);
    for (int k = 0; k < NEQ; ++k) {
        d[k] = (Fa[k] - Fb[k]) / (2.0 * h) + (Ga[k] - Gb[k]) / (2.0 * h);
    }
}

/*
   C4: R / Omega -> div(Fc) when the mesh is refined. Field set at the cell centres
   (average of the 4 corners); j ghost centres extrapolated linearly from the interior
   (like the professor's example code); i ghosts copied across the cut.
   Error measured in cells whose centre is within 3 chords of mid-chord (the farfield
   cells are far too large for this smooth field), area-weighted L1 and L2 norms.

   Expected orders: first order -> 1. Second order -> 2 in L1, but less in L2: where one
   variable of the field has an extremum along a grid line, a and b of the limiter have
   opposite signs and the limiter moves the face value by O(h^2) differently on the two
   faces of a cell, so those cells keep an O(h) error. They form lines, O(1/h) cells, which
   is O(h^2) in L1 but only O(h^1.5) in L2. (Checked while writing this test: with the
   limiter switched off, eps -> infinity, which gives the unlimited kappa = 1/3 slope
   (b + 2a)/3, the L2 order goes back to 2.)
   err[o][0] = L1, err[o][1] = L2, o = 0 first order, o = 1 second order.
*/
static void mms_errors(const euler::Mesh& m, euler::Flow& f, double err[2][2])
{
    std::vector<double> xc(m.vol.size(), 0.0); // cell centre x (ghosts in j included)
    std::vector<double> yc(m.vol.size(), 0.0); // cell centre y
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const int c = m.cell(i, j);
            xc[c] = 0.25 * (m.x[m.node(i, j)] + m.x[m.node(i + 1, j)] + m.x[m.node(i + 1, j + 1)] + m.x[m.node(i, j + 1)]);
            yc[c] = 0.25 * (m.y[m.node(i, j)] + m.y[m.node(i + 1, j)] + m.y[m.node(i + 1, j + 1)] + m.y[m.node(i, j + 1)]);
        }
    }
    for (int i = 0; i < m.nic; ++i) {
        for (int g = 1; g <= NG; ++g) {
            const int in = m.cell(i, -g + 1);      // next cell towards the interior
            const int in2 = m.cell(i, -g + 2);     // the one after
            const int gc = m.cell(i, -g);          // wall ghost
            xc[gc] = 2.0 * xc[in] - xc[in2];
            yc[gc] = 2.0 * yc[in] - yc[in2];
            const int out = m.cell(i, m.njc - 2 + g);  // next cell towards the interior (farfield side)
            const int out2 = m.cell(i, m.njc - 3 + g); // the one after
            const int gf = m.cell(i, m.njc - 1 + g);   // farfield ghost
            xc[gf] = 2.0 * xc[out] - xc[out2];
            yc[gf] = 2.0 * yc[out] - yc[out2];
        }
    }
    for (int j = -NG; j < m.njc + NG; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const int c = m.cell(i, j);
            double q[NEQ];
            double w[NEQ];
            mms_prim(xc[c], yc[c], q);
            euler::prim_to_cons(q, w);
            for (int k = 0; k < NEQ; ++k) {
                f.W[k][c] = w[k];
            }
        }
    }
    fill_icut_ghosts(m, f);

    for (int pass = 0; pass < 2; ++pass) {
        const euler::Order order = pass == 0 ? euler::Order::First : euler::Order::Second;
        euler::compute_residual(m, f, order);
        double sum1 = 0.0; // sum of |err| * area
        double sum2 = 0.0; // sum of err^2 * area
        double area = 0.0; // sum of area
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int c = m.cell(i, j);
                if (std::hypot(xc[c] - 0.5, yc[c]) > 3.0) {
                    continue;
                }
                double d[NEQ];
                mms_div(xc[c], yc[c], d);
                for (int k = 0; k < NEQ; ++k) {
                    const double e = f.R[k][c] / m.vol[c] - d[k];
                    sum1 += std::abs(e) * m.vol[c];
                    sum2 += e * e * m.vol[c];
                }
                area += m.vol[c];
            }
        }
        err[pass][0] = sum1 / area;
        err[pass][1] = std::sqrt(sum2 / area);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mesh.x> [<mesh.x> ...]  (nested grids, coarse to fine)\n", argv[0]);
        return 1;
    }

    try {
        test_roe_flux();
        test_muscl();

        std::printf("\n=== C. Residual on the grids\n");
        std::vector<int> sizes;  // ni of each grid
        std::vector<std::array<double, 4>> errs; // C4 errors per grid: L1 1st, L2 1st, L1 2nd, L2 2nd
        euler::Mesh last;        // finest grid, kept for the timing
        euler::Flow flast;
        for (int a = 1; a < argc; ++a) {
            euler::Mesh m = euler::read_plot3d(argv[a]);
            euler::compute_metrics(m);
            euler::Flow f;
            euler::allocate_flow(m, f);
            std::printf("--- %s (%d x %d nodes)\n", argv[a], m.ni, m.nj);
            test_uniform(m, f);
            test_random_field(m, f);
            double err[2][2];
            mms_errors(m, f, err);
            sizes.push_back(m.ni);
            errs.push_back({err[0][0], err[0][1], err[1][0], err[1][1]});
            last = std::move(m);
            flast = std::move(f);
        }

        std::printf("\n--- C4 residual error |R/Omega - div Fc| (cells within 3 chords), order between grids in ()\n");
        std::printf("    grid    1st order L1     1st order L2     2nd order L1     2nd order L2\n");
        for (size_t g = 0; g < sizes.size(); ++g) {
            std::printf("    %4d ", sizes[g]);
            for (int c = 0; c < 4; ++c) {
                if (g == 0) {
                    std::printf("   %.3e       ", errs[g][c]);
                } else {
                    std::printf("   %.3e (%4.2f)", errs[g][c], std::log2(errs[g - 1][c] / errs[g][c]));
                }
            }
            std::printf("\n");
        }
        if (sizes.size() >= 3) {
            const size_t n = sizes.size() - 1;
            const double p1 = std::log2(errs[n - 1][0] / errs[n][0]);      // first order, L1, finest pair
            const double p2 = std::log2(errs[n - 1][2] / errs[n][2]);      // second order, L1, finest pair
            const double p2b = std::log2(errs[n - 2][2] / errs[n - 1][2]); // second order, L1, the pair before
            check(p1 > 0.9, "C4 first-order residual: L1 order %.2f on the finest pair (expect 1)", p1);
            check(p2 > 1.85 && p2b > 1.85,
                  "C4 second-order residual: L1 order %.2f and %.2f on the two finest pairs (expect 2)", p2b, p2);
            std::printf("  [INFO] C4 second-order L2 order on the finest pair: %.2f (limiter at smooth extrema, see mms_errors)\n",
                        std::log2(errs[n - 1][3] / errs[n][3]));
        } else {
            std::printf("  [SKIP] C4 needs at least 3 nested grids\n");
        }

        // Timing of one second-order residual on the finest grid (uniform flow, no NaN anywhere).
        {
            const double q[NEQ] = {1.0, 0.5, 0.0, 1.0 / GAMMA};
            fill_uniform(flast, q);
            const int nrep = 20; // residual evaluations timed
            const auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < nrep; ++r) {
                euler::compute_residual(last, flast, euler::Order::Second);
            }
            const auto t1 = std::chrono::steady_clock::now();
            const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / nrep;
            std::printf("\n  [INFO] one 2nd-order residual on %d x %d: %.2f ms (%.1f ns per cell)\n", last.ni, last.nj, ms,
                        1e6 * ms / (static_cast<double>(last.nic) * last.njc));
        }
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }

    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
