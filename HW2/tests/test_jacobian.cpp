/*
   Implicit operator check: flux Jacobians and first-order block matrix (Roe_scheme_accel slides 25-30).

   A. Physical flux Jacobian A = dFc/dW
      A1 A vs central finite differences of physical_flux in conservative variables (Mach up to 3)
      A2 homogeneity A(W) W = Fc(W) (the Euler flux is homogeneous of degree 1)
   B. Roe matrix |A~| = T~ |Lambda~| T~^-1
      B1 1/2 (FL + FR) - 1/2 |A~| (WR - WL) = roe_flux (Harten on and off)
      B2 all eigenvalues > delta (supersonic along n) -> |A~| (WR - WL) = FR - FL (Roe property)
      B3 qL = qR supersonic -> |A~| x = A x;  qL = qR subsonic, no fix -> |A~|^2 = A^2
   C. Frozen-|A~| face Jacobians at qL = qR vs finite differences of roe_flux (slide 27)
   D. Assembled first-order Jacobian on the airfoil grids
      D1 uniform free stream (jmin Farfield and Wall): M v vs finite-difference directional
         derivative of the first-order residual (ghosts follow the interior through apply_bc),
         per region: wall row, farfield row, TE cut columns, interior
      D2 (info) same on a randomly perturbed free stream: error of the frozen-|A~| approximation
      D3 diag_add lands exactly on the diagonal of D and nowhere else

   Usage: test_jacobian <mesh.x> [<mesh.x> ...]
   Exit code = number of failed checks.
*/

#include "bc.h"
#include "blockmat.h"
#include "flow.h"
#include "jacobian.h"
#include "mesh.h"
#include "residual.h"
#include "roe.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <random>
#include <string>
#include <vector>

using euler::GAMMA;
using euler::NEQ;

constexpr double PI = 3.14159265358979323846; // MSVC has no M_PI by default
constexpr int NB = NEQ * NEQ;                 // entries of one 4x4 block

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

static std::mt19937 g_rng(12345); // fixed seed: reproducible tests

// Uniform random number in [a, b].
static double rnd(double a, double b)
{
    std::uniform_real_distribution<double> d(a, b);
    return d(g_rng);
}

// Random primitive state: rho in [0.3, 3], p in [0.3, 3]/gamma, speed = Mach * c, Mach in [0, mach_max].
static void random_state(double mach_max, double q[NEQ])
{
    q[0] = rnd(0.3, 3.0);
    q[3] = rnd(0.3, 3.0) / GAMMA;
    const double c = std::sqrt(GAMMA * q[3] / q[0]); // speed of sound
    const double mach = rnd(0.0, mach_max);          // Mach number
    const double th = rnd(0.0, 2.0 * PI);            // flow direction
    q[1] = mach * c * std::cos(th);
    q[2] = mach * c * std::sin(th);
}

// Random unit normal.
static void random_normal(double& nx, double& ny)
{
    const double th = rnd(0.0, 2.0 * PI); // normal direction
    nx = std::cos(th);
    ny = std::sin(th);
}

// Physical flux as a function of the conservative state.
static void flux_of_w(const double w[NEQ], double nx, double ny, double F[NEQ])
{
    double q[NEQ]; // primitive state
    euler::cons_to_prim(w, q);
    euler::physical_flux(q, nx, ny, F);
}

// y = A x for one 4x4 row-major block.
static void mat4_vec(const double A[NB], const double x[NEQ], double y[NEQ])
{
    for (int r = 0; r < NEQ; ++r) {
        y[r] = 0.0;
        for (int c = 0; c < NEQ; ++c) {
            y[r] += A[r * NEQ + c] * x[c];
        }
    }
}

// C = A B for 4x4 row-major blocks.
static void mat4_mul(const double A[NB], const double B[NB], double C[NB])
{
    for (int r = 0; r < NEQ; ++r) {
        for (int c = 0; c < NEQ; ++c) {
            C[r * NEQ + c] = 0.0;
            for (int k = 0; k < NEQ; ++k) {
                C[r * NEQ + c] += A[r * NEQ + k] * B[k * NEQ + c];
            }
        }
    }
}

// Largest |a_k| over n entries.
static double max_abs(const double* a, int n)
{
    double v = 0.0;
    for (int k = 0; k < n; ++k) {
        v = std::max(v, std::abs(a[k]));
    }
    return v;
}

// Largest |a_k - b_k| over n entries.
static double max_abs_diff(const double* a, const double* b, int n)
{
    double v = 0.0;
    for (int k = 0; k < n; ++k) {
        v = std::max(v, std::abs(a[k] - b[k]));
    }
    return v;
}

/*
   Central finite-difference Jacobian of a flux G(w) (G = any function of the conservative state):
   J[r * NEQ + c] = (G_r(w + h e_c) - G_r(w - h e_c)) / (2 h), h = h_rel (|w_c| + rho).
*/
template <class Func>
static void fd_jacobian(Func G, const double w[NEQ], double h_rel, double J[NB])
{
    for (int c = 0; c < NEQ; ++c) {
        double wp[NEQ]; // w + h e_c
        double wm[NEQ]; // w - h e_c
        for (int k = 0; k < NEQ; ++k) {
            wp[k] = w[k];
            wm[k] = w[k];
        }
        const double h = h_rel * (std::abs(w[c]) + std::abs(w[0])); // step
        wp[c] += h;
        wm[c] -= h;
        double Fp[NEQ]; // G(w + h e_c)
        double Fm[NEQ]; // G(w - h e_c)
        G(wp, Fp);
        G(wm, Fm);
        for (int r = 0; r < NEQ; ++r) {
            J[r * NEQ + c] = (Fp[r] - Fm[r]) / (wp[c] - wm[c]);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// A. Physical flux Jacobian
// ---------------------------------------------------------------------------------------------
static void test_flux_jacobian()
{
    std::printf("A. Physical flux Jacobian A = dFc/dW\n");
    double err_fd = 0.0;   // worst relative error vs finite differences
    double err_homo = 0.0; // worst relative error of A W - Fc
    for (int t = 0; t < 500; ++t) {
        double q[NEQ]; // random state
        double w[NEQ]; // same, conservative
        double nx;     // normal x
        double ny;     // normal y
        random_state(3.0, q);
        random_normal(nx, ny);
        euler::prim_to_cons(q, w);

        double A[NB];   // analytic Jacobian
        double Afd[NB]; // finite-difference Jacobian
        euler::flux_jacobian(q, nx, ny, A);
        fd_jacobian([&](const double* ww, double* F) { flux_of_w(ww, nx, ny, F); }, w, 1e-6, Afd);
        err_fd = std::max(err_fd, max_abs_diff(A, Afd, NB) / max_abs(A, NB));

        double AW[NEQ]; // A(W) W
        double F[NEQ];  // Fc(W)
        mat4_vec(A, w, AW);
        euler::physical_flux(q, nx, ny, F);
        err_homo = std::max(err_homo, max_abs_diff(AW, F, NEQ) / max_abs(F, NEQ));
    }
    check(err_fd < 1e-8, "A1 A vs finite differences, 500 states Mach <= 3: max rel err %.2e", err_fd);
    check(err_homo < 1e-13, "A2 homogeneity A(W) W = Fc(W): max rel err %.2e", err_homo);
}

// ---------------------------------------------------------------------------------------------
// B. Roe matrix
// ---------------------------------------------------------------------------------------------

// Roe average normal velocity V~ and speed of sound c~ (same formulas as roe.cpp).
static void roe_vc(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double& V, double& c)
{
    const double HL = GAMMA / (GAMMA - 1.0) * qL[3] / qL[0] + 0.5 * (qL[1] * qL[1] + qL[2] * qL[2]); // left H
    const double HR = GAMMA / (GAMMA - 1.0) * qR[3] / qR[0] + 0.5 * (qR[1] * qR[1] + qR[2] * qR[2]); // right H
    const double sL = std::sqrt(qL[0]);             // sqrt rho_L
    const double sR = std::sqrt(qR[0]);             // sqrt rho_R
    const double wL = sL / (sL + sR);               // left weight
    const double wR = sR / (sL + sR);               // right weight
    const double u = wL * qL[1] + wR * qR[1];       // u~
    const double v = wL * qL[2] + wR * qR[2];       // v~
    const double H = wL * HL + wR * HR;             // H~
    c = std::sqrt((GAMMA - 1.0) * (H - 0.5 * (u * u + v * v)));
    V = u * nx + v * ny;
}

// Random state moving along +n at a normal Mach number in [mn_min, mn_max], tangential Mach in [-1, 1].
static void random_normal_state(double nx, double ny, double mn_min, double mn_max, double q[NEQ])
{
    q[0] = rnd(0.3, 3.0);
    q[3] = rnd(0.3, 3.0) / GAMMA;
    const double c = std::sqrt(GAMMA * q[3] / q[0]); // speed of sound
    const double mn = rnd(mn_min, mn_max);           // normal Mach number
    const double mt = rnd(-1.0, 1.0);                // tangential Mach number
    q[1] = c * (mn * nx - mt * ny);
    q[2] = c * (mn * ny + mt * nx);
}

static void test_roe_matrix()
{
    std::printf("B. Roe matrix |A~|\n");

    // B1: |A~| dW is the roe_flux dissipation, with and without the Harten fix.
    for (double hf : {euler::HARTEN_FRACTION, 0.0}) {
        double err = 0.0; // worst relative error
        for (int t = 0; t < 500; ++t) {
            double qL[NEQ]; // left state
            double qR[NEQ]; // right state
            double nx;      // normal x
            double ny;      // normal y
            random_state(3.0, qL);
            random_state(3.0, qR);
            random_normal(nx, ny);

            double wL[NEQ];   // left, conservative
            double wR[NEQ];   // right, conservative
            double dW[NEQ];   // WR - WL
            double absA[NB];  // |A~|
            double D[NEQ];    // |A~| dW
            double FL[NEQ];   // Fc(WL)
            double FR[NEQ];   // Fc(WR)
            double Fm[NEQ];   // flux rebuilt from |A~|
            double Froe[NEQ]; // roe_flux
            euler::prim_to_cons(qL, wL);
            euler::prim_to_cons(qR, wR);
            for (int k = 0; k < NEQ; ++k) {
                dW[k] = wR[k] - wL[k];
            }
            euler::roe_abs_matrix(qL, qR, nx, ny, absA, hf);
            mat4_vec(absA, dW, D);
            euler::physical_flux(qL, nx, ny, FL);
            euler::physical_flux(qR, nx, ny, FR);
            for (int k = 0; k < NEQ; ++k) {
                Fm[k] = 0.5 * (FL[k] + FR[k] - D[k]);
            }
            euler::roe_flux(qL, qR, nx, ny, Froe, hf);
            const double scale = std::max({max_abs(FL, NEQ), max_abs(FR, NEQ), max_abs(D, NEQ)}); // magnitude
            err = std::max(err, max_abs_diff(Fm, Froe, NEQ) / scale);
        }
        check(err < 1e-13, "B1 1/2(FL+FR) - 1/2|A~|dW = roe_flux, harten_fraction %.1f, 500 pairs: max rel err %.2e",
              hf, err);
    }

    // B2: supersonic Roe state along n (all eigenvalues > delta) -> |A~| = A~, |A~| dW = FR - FL.
    {
        double err = 0.0; // worst relative error
        int used = 0;     // pairs that satisfied the eigenvalue condition
        for (int t = 0; t < 500; ++t) {
            double nx; // normal x
            double ny; // normal y
            random_normal(nx, ny);
            double qL[NEQ]; // left state
            double qR[NEQ]; // right state
            random_normal_state(nx, ny, 1.5, 3.0, qL);
            random_normal_state(nx, ny, 1.5, 3.0, qR);
            double V; // V~
            double c; // c~
            roe_vc(qL, qR, nx, ny, V, c);
            if (V - c <= euler::HARTEN_FRACTION * c) {
                continue;
            }
            ++used;

            double wL[NEQ];  // left, conservative
            double wR[NEQ];  // right, conservative
            double dW[NEQ];  // WR - WL
            double absA[NB]; // |A~|
            double D[NEQ];   // |A~| dW
            double FL[NEQ];  // Fc(WL)
            double FR[NEQ];  // Fc(WR)
            double dF[NEQ];  // FR - FL
            euler::prim_to_cons(qL, wL);
            euler::prim_to_cons(qR, wR);
            for (int k = 0; k < NEQ; ++k) {
                dW[k] = wR[k] - wL[k];
            }
            euler::roe_abs_matrix(qL, qR, nx, ny, absA);
            mat4_vec(absA, dW, D);
            euler::physical_flux(qL, nx, ny, FL);
            euler::physical_flux(qR, nx, ny, FR);
            for (int k = 0; k < NEQ; ++k) {
                dF[k] = FR[k] - FL[k];
            }
            const double scale = std::max({max_abs(FL, NEQ), max_abs(FR, NEQ)}); // magnitude
            err = std::max(err, max_abs_diff(D, dF, NEQ) / scale);
        }
        check(used > 100 && err < 1e-13, "B2 all eigenvalues > delta: |A~| dW = FR - FL, %d pairs: max rel err %.2e",
              used, err);
    }

    // B3a: equal supersonic states -> |A~| = A(q).
    {
        double err = 0.0; // worst relative error
        for (int t = 0; t < 500; ++t) {
            double nx; // normal x
            double ny; // normal y
            random_normal(nx, ny);
            double q[NEQ]; // state, normal Mach in [1.2, 3] (V - c >= 0.17 c > delta)
            random_normal_state(nx, ny, 1.2, 3.0, q);
            double absA[NB]; // |A~|(q, q)
            double A[NB];    // A(q)
            euler::roe_abs_matrix(q, q, nx, ny, absA);
            euler::flux_jacobian(q, nx, ny, A);
            double x[NEQ];  // random vector
            double y1[NEQ]; // |A~| x
            double y2[NEQ]; // A x
            for (int k = 0; k < NEQ; ++k) {
                x[k] = rnd(-1.0, 1.0);
            }
            mat4_vec(absA, x, y1);
            mat4_vec(A, x, y2);
            err = std::max(err, max_abs_diff(y1, y2, NEQ) / max_abs(A, NB));
        }
        check(err < 1e-13, "B3a qL = qR supersonic: |A~| x = A(q) x, 500 states: max rel err %.2e", err);
    }

    // B3b: equal subsonic states, no Harten fix -> |A~|^2 = T Lambda^2 T^-1 = A^2.
    {
        double err = 0.0; // worst relative error
        for (int t = 0; t < 500; ++t) {
            double nx; // normal x
            double ny; // normal y
            random_normal(nx, ny);
            double q[NEQ]; // subsonic state
            random_state(0.95, q);
            double absA[NB];  // |A~|(q, q)
            double A[NB];     // A(q)
            double absA2[NB]; // |A~|^2
            double A2[NB];    // A^2
            euler::roe_abs_matrix(q, q, nx, ny, absA, 0.0);
            euler::flux_jacobian(q, nx, ny, A);
            mat4_mul(absA, absA, absA2);
            mat4_mul(A, A, A2);
            err = std::max(err, max_abs_diff(absA2, A2, NB) / max_abs(A2, NB));
        }
        check(err < 1e-12, "B3b qL = qR subsonic, no fix: |A~|^2 = A^2, 500 states: max rel err %.2e", err);
    }
}

// ---------------------------------------------------------------------------------------------
// C. Face Jacobians
// ---------------------------------------------------------------------------------------------
static void test_face_jacobians()
{
    std::printf("C. Frozen-|A~| face Jacobians at qL = qR\n");
    double errL = 0.0; // worst relative error of dF/dW_L
    double errR = 0.0; // worst relative error of dF/dW_R
    for (int t = 0; t < 500; ++t) {
        double q[NEQ]; // common state
        double w[NEQ]; // same, conservative
        double nx;     // normal x
        double ny;     // normal y
        random_state(3.0, q);
        random_normal(nx, ny);
        euler::prim_to_cons(q, w);

        double dFdL[NB]; // analytic, frozen |A~|
        double dFdR[NB];
        euler::roe_face_jacobians(q, q, nx, ny, dFdL, dFdR);

        double JL[NB]; // finite differences of roe_flux w.r.t. W_L
        double JR[NB]; // and w.r.t. W_R
        fd_jacobian(
            [&](const double* ww, double* F) {
                double qq[NEQ]; // perturbed left state
                euler::cons_to_prim(ww, qq);
                euler::roe_flux(qq, q, nx, ny, F);
            },
            w, 1e-6, JL);
        fd_jacobian(
            [&](const double* ww, double* F) {
                double qq[NEQ]; // perturbed right state
                euler::cons_to_prim(ww, qq);
                euler::roe_flux(q, qq, nx, ny, F);
            },
            w, 1e-6, JR);
        const double scale = std::max(max_abs(dFdL, NB), max_abs(dFdR, NB)); // magnitude
        errL = std::max(errL, max_abs_diff(dFdL, JL, NB) / scale);
        errR = std::max(errR, max_abs_diff(dFdR, JR, NB) / scale);
    }
    check(errL < 1e-6, "C  dF/dW_L vs finite differences of roe_flux, 500 states: max rel err %.2e", errL);
    check(errR < 1e-6, "C  dF/dW_R vs finite differences of roe_flux, 500 states: max rel err %.2e", errR);
}

// ---------------------------------------------------------------------------------------------
// D. Assembled Jacobian on a grid
// ---------------------------------------------------------------------------------------------

// y = M x for the block matrix (blockmat.h layout), local helper of the test.
static void test_matvec(const euler::BlockMatrix& M, const std::vector<double>& x, std::vector<double>& y)
{
    y.assign(x.size(), 0.0);
    for (int j = 0; j < M.njc; ++j) {
        for (int i = 0; i < M.nic; ++i) {
            const int p = M.index(i, j); // block row
            double t[NEQ];               // one block product
            mat4_vec(&M.D[p * NB], &x[p * NEQ], t);
            for (int k = 0; k < NEQ; ++k) {
                y[p * NEQ + k] += t[k];
            }
            mat4_vec(&M.Wb[p * NB], &x[M.index(M.west(i), j) * NEQ], t);
            for (int k = 0; k < NEQ; ++k) {
                y[p * NEQ + k] += t[k];
            }
            mat4_vec(&M.Eb[p * NB], &x[M.index(M.east(i), j) * NEQ], t);
            for (int k = 0; k < NEQ; ++k) {
                y[p * NEQ + k] += t[k];
            }
            if (j > 0) {
                mat4_vec(&M.Sb[p * NB], &x[M.index(i, j - 1) * NEQ], t);
                for (int k = 0; k < NEQ; ++k) {
                    y[p * NEQ + k] += t[k];
                }
            }
            if (j < M.njc - 1) {
                mat4_vec(&M.Nb[p * NB], &x[M.index(i, j + 1) * NEQ], t);
                for (int k = 0; k < NEQ; ++k) {
                    y[p * NEQ + k] += t[k];
                }
            }
        }
    }
}

/*
   First-order residual of the field base + s v (physical cells only, ghosts refilled by
   apply_bc so that they follow the interior), compact numbering.
*/
static void residual_of(const euler::Mesh& m, const euler::Flow& base, const euler::Freestream& fs,
                        euler::JminBc jmin, const std::vector<double>& v, double s, std::vector<double>& R)
{
    euler::Flow g = base; // working copy
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const int p = j * m.nic + i; // compact index
            for (int k = 0; k < NEQ; ++k) {
                g.W[k][m.cell(i, j)] += s * v[p * NEQ + k];
            }
        }
    }
    euler::apply_bc(m, g, fs, jmin);
    euler::compute_residual(m, g, euler::Order::First);
    R.assign(static_cast<size_t>(m.nic) * m.njc * NEQ, 0.0);
    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            const int p = j * m.nic + i; // compact index
            for (int k = 0; k < NEQ; ++k) {
                R[p * NEQ + k] = g.R[k][m.cell(i, j)];
            }
        }
    }
}

// Max error and max reference value of one region of the grid.
struct RegionErr {
    double err = 0.0; // max |M v - FD|
    double ref = 0.0; // max |FD|
    double rel() const { return ref > 0.0 ? err / ref : err; }
};

/*
   Compare M v with the central finite-difference directional derivative
   (R1(W + eps v) - R1(W - eps v)) / (2 eps) on the field f (ghosts filled).
   Regions: 0 wall row j = 0, 1 farfield row j = njc-1, 2 TE cut columns i = 0, nic-1, 3 interior.
*/
static void compare_matvec(const euler::Mesh& m, const euler::Flow& f, const euler::Freestream& fs,
                           euler::JminBc jmin, RegionErr reg[4])
{
    euler::Flow g = f; // field with up-to-date primitives
    euler::compute_primitives(g);
    euler::BlockMatrix M; // assembled Jacobian
    euler::assemble_first_order_jacobian(m, g, fs, jmin, {}, M);

    const int n = m.nic * m.njc * NEQ; // number of unknowns
    std::vector<double> v(n);          // random direction
    const double scale[NEQ] = {1.0, 1.0, 1.0, 2.0}; // typical size of each conservative variable
    for (int r = 0; r < n; ++r) {
        v[r] = rnd(-1.0, 1.0) * scale[r % NEQ];
    }
    std::vector<double> y; // M v
    test_matvec(M, v, y);

    const double eps = 1e-6;  // finite-difference step along v
    std::vector<double> Rp;   // R1(W + eps v)
    std::vector<double> Rm;   // R1(W - eps v)
    residual_of(m, g, fs, jmin, v, eps, Rp);
    residual_of(m, g, fs, jmin, v, -eps, Rm);

    for (int j = 0; j < m.njc; ++j) {
        for (int i = 0; i < m.nic; ++i) {
            int region = 3; // interior by default
            if (j == 0) {
                region = 0;
            } else if (j == m.njc - 1) {
                region = 1;
            } else if (i == 0 || i == m.nic - 1) {
                region = 2;
            }
            const int p = j * m.nic + i; // compact index
            for (int k = 0; k < NEQ; ++k) {
                const double fd = (Rp[p * NEQ + k] - Rm[p * NEQ + k]) / (2.0 * eps); // FD derivative
                reg[region].err = std::max(reg[region].err, std::abs(y[p * NEQ + k] - fd));
                reg[region].ref = std::max(reg[region].ref, std::abs(fd));
            }
        }
    }
}

static const char* REGION_NAME[4] = {"wall/jmin row", "farfield row ", "TE cut cols  ", "interior     "};

static void test_grid(const std::string& path)
{
    std::printf("D. Assembled first-order Jacobian on %s\n", path.c_str());
    euler::Mesh m = euler::read_plot3d(path);
    euler::compute_metrics(m);
    euler::Flow f;
    euler::allocate_flow(m, f);
    std::printf("   %d x %d cells\n", m.nic, m.njc);

    const euler::Freestream fs = euler::make_freestream(0.8, 1.25); // M = 0.8, alpha = 1.25 deg

    for (euler::JminBc jmin : {euler::JminBc::Farfield, euler::JminBc::Wall}) {
        const char* bcname = jmin == euler::JminBc::Wall ? "Wall" : "Farfield"; // jmin name

        // D1: uniform free stream, interior faces have equal states -> exact.
        euler::init_freestream(f, fs);
        euler::apply_bc(m, f, fs, jmin);
        RegionErr reg[4];
        compare_matvec(m, f, fs, jmin, reg);
        for (int r = 0; r < 4; ++r) {
            check(reg[r].rel() < 1e-6, "D1 free stream, jmin %-8s, %s: max |Mv - FD| / max |FD| = %.2e (max |FD| %.2e)",
                  bcname, REGION_NAME[r], reg[r].rel(), reg[r].ref);
        }

        // D2 (info): free stream + 1% random noise on the interior cells.
        euler::init_freestream(f, fs);
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                for (int k = 0; k < NEQ; ++k) {
                    double& w = f.W[k][m.cell(i, j)]; // conservative variable k of the cell
                    w += 0.01 * rnd(-1.0, 1.0) * (k == 3 ? 2.0 : 1.0);
                }
            }
        }
        euler::apply_bc(m, f, fs, jmin);
        RegionErr regn[4];
        compare_matvec(m, f, fs, jmin, regn);
        for (int r = 0; r < 4; ++r) {
            std::printf("  [INFO] D2 1%% noise,   jmin %-8s, %s: max |Mv - FD| / max |FD| = %.2e (frozen |A~|)\n",
                        bcname, REGION_NAME[r], regn[r].rel());
        }
    }

    // D3: diag_add on the diagonal of D only.
    {
        euler::init_freestream(f, fs);
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                for (int k = 0; k < NEQ; ++k) {
                    f.W[k][m.cell(i, j)] *= 1.0 + 0.01 * rnd(-1.0, 1.0);
                }
            }
        }
        euler::apply_bc(m, f, fs, euler::JminBc::Wall);
        euler::compute_primitives(f);
        const int nc = m.nic * m.njc;  // number of physical cells
        std::vector<double> d(nc);     // diagonal additions
        for (double& x : d) {
            x = rnd(0.1, 10.0);
        }
        euler::BlockMatrix M0; // without diag_add
        euler::BlockMatrix M1; // with diag_add
        euler::assemble_first_order_jacobian(m, f, fs, euler::JminBc::Wall, {}, M0);
        euler::assemble_first_order_jacobian(m, f, fs, euler::JminBc::Wall, d, M1);
        double err_diag = 0.0; // worst |(D1 - D0) - d| / (|D0| + d) on the diagonal
        bool others_equal = M0.Wb == M1.Wb && M0.Eb == M1.Eb && M0.Sb == M1.Sb && M0.Nb == M1.Nb; // bitwise
        for (int p = 0; p < nc; ++p) {
            for (int r = 0; r < NEQ; ++r) {
                for (int c = 0; c < NEQ; ++c) {
                    const int e = p * NB + r * NEQ + c; // entry index
                    if (r == c) {
                        const double diff = (M1.D[e] - M0.D[e]) - d[p]; // should be round-off
                        err_diag = std::max(err_diag, std::abs(diff) / (std::abs(M0.D[e]) + d[p]));
                    } else if (M1.D[e] != M0.D[e]) {
                        others_equal = false;
                    }
                }
            }
        }
        check(err_diag < 1e-15 && others_equal,
              "D3 diag_add on the diagonal of D only: rel err %.2e, off-diagonal / neighbour blocks unchanged: %s",
              err_diag, others_equal ? "yes" : "NO");
    }

    // Size check of diag_add.
    {
        bool threw = false; // did a wrong-size diag_add throw?
        euler::BlockMatrix M;
        try {
            euler::assemble_first_order_jacobian(m, f, fs, euler::JminBc::Wall, std::vector<double>(3, 1.0), M);
        } catch (const std::exception&) {
            threw = true;
        }
        check(threw, "D3 diag_add of the wrong size is rejected");
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::printf("usage: test_jacobian <mesh.x> [<mesh.x> ...]\n");
        return 1;
    }
    try {
        test_flux_jacobian();
        test_roe_matrix();
        test_face_jacobians();
        for (int a = 1; a < argc; ++a) {
            test_grid(argv[a]);
        }
    } catch (const std::exception& e) {
        std::printf("  [FAIL] exception: %s\n", e.what());
        ++g_fail;
    }
    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
