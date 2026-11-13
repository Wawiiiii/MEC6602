/*
   Linear solver check (Roe_scheme_accel slides 32-40), no mesh needed: the BlockMatrix objects are filled directly.

   A. block_matvec = dense product of the matrix assembled from the blocks (periodic i, no S at j = 0, no N at j = njc-1)
   B. invert4 / factor_diagonal / diagonal_solve: D x = b to round-off, singular block detected
   C. block SGS (slides 33-35) on a block-diagonally-dominant matrix:
      C1 residual decreases monotonically with the number of sweeps (omega = 1 and under-relaxed omega = 0.7)
      C2 many sweeps -> dense LU solution
      C3 3 sweeps + 2 more sweeps = 5 sweeps (x really is the initial guess)
      C4 block-diagonal matrix: one sweep is exact
   D. GMRES with right preconditioning (slide 39) on a nonsymmetric, weakly dominant matrix:
      no preconditioner vs 2 SGS sweeps, restart 30 and restart 5, compared with dense LU
   E. GMRES special cases: b = 0, exact initial guess, happy breakdown, max_iter reached
   F. matrix-free operators: dense-matrix lambda, finite-difference Jacobian-vector product (JFNK, slide 40)
   G. timing of one matvec and one SGS sweep on a 1024 x 1024 and 2048 x 2048 cell matrix

   Exit code = number of failed checks.
*/

#include "blockmat.h"
#include "linsolve.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <new>
#include <random>
#include <stdexcept>
#include <vector>

using euler::BlockMatrix;
using euler::GmresResult;
using euler::NEQ;

constexpr int NB = NEQ * NEQ; // entries in one 4x4 block

static int g_pass = 0; // number of passed checks
static int g_fail = 0; // number of failed checks

// Print one check line and count it.
static void check(bool ok, const char* fmt, ...)
{
    std::printf(ok ? "  [PASS] " : "  [FAIL] ");
    va_list args; // printf-style arguments
    va_start(args, fmt);
    std::vprintf(fmt, args);
    va_end(args);
    std::printf("\n");
    if (ok) {
        ++g_pass;
    }
    else {
        ++g_fail;
    }
}

// Elapsed milliseconds since t0.
static double ms_since(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// ---- matrix generators ----

/*
   Random matrix shaped like an implicit CFD operator, block-diagonally dominant:
     D  = ddiag I + U(-1, 1) doff     (diagonal block: V/dt I + dR/dW)
     Wb, Eb, Sb, Nb = U(-1, 1) dnb    (neighbour blocks)
   South blocks at j = 0 and north blocks at j = njc-1 stay zero (they are not used anyway).
*/
static void make_dominant(BlockMatrix& M, int nic, int njc, unsigned seed, double ddiag, double doff, double dnb)
{
    std::mt19937_64 gen(seed);                             // random generator
    std::uniform_real_distribution<double> U(-1.0, 1.0);   // uniform in [-1, 1]
    M.resize(nic, njc);
    const size_t nc = static_cast<size_t>(M.ncell());      // number of cells
    for (size_t p = 0; p < nc; ++p) {
        const int j = static_cast<int>(p) / nic;           // row index of the cell
        for (int k = 0; k < NB; ++k) {
            M.D[p * NB + k] = doff * U(gen) + ((k % (NEQ + 1) == 0) ? ddiag : 0.0);
            M.Wb[p * NB + k] = dnb * U(gen);
            M.Eb[p * NB + k] = dnb * U(gen);
            M.Sb[p * NB + k] = (j > 0) ? dnb * U(gen) : 0.0;
            M.Nb[p * NB + k] = (j < njc - 1) ? dnb * U(gen) : 0.0;
        }
    }
}

/*
   Nonsymmetric, only weakly dominant matrix (upwind convection-diffusion with a small "V/dt" term sigma):
     Wb = -cw I + pert U,  Eb = -ce I + pert U,  Sb = -cs I + pert U,  Nb = -cn I + pert U
     D  = (cw + ce + cs + cn + sigma) I + pert U
   cw != ce and cs != cn make it nonsymmetric; small sigma makes it hard for unpreconditioned GMRES.
*/
static void make_convdiff(BlockMatrix& M, int nic, int njc, unsigned seed, double sigma, double pert)
{
    const double cw = 1.5; // west coupling (upwind side)
    const double ce = 0.5; // east coupling
    const double cs = 1.0; // south coupling
    const double cn = 0.4; // north coupling
    std::mt19937_64 gen(seed);                           // random generator
    std::uniform_real_distribution<double> U(-1.0, 1.0); // uniform in [-1, 1]
    M.resize(nic, njc);
    const size_t nc = static_cast<size_t>(M.ncell());    // number of cells
    for (size_t p = 0; p < nc; ++p) {
        const int j = static_cast<int>(p) / nic;         // row index of the cell
        for (int k = 0; k < NB; ++k) {
            const double id = (k % (NEQ + 1) == 0) ? 1.0 : 0.0; // identity entry
            M.D[p * NB + k] = (cw + ce + cs + cn + sigma) * id + pert * U(gen);
            M.Wb[p * NB + k] = -cw * id + pert * U(gen);
            M.Eb[p * NB + k] = -ce * id + pert * U(gen);
            M.Sb[p * NB + k] = (j > 0) ? -cs * id + pert * U(gen) : 0.0;
            M.Nb[p * NB + k] = (j < njc - 1) ? -cn * id + pert * U(gen) : 0.0;
        }
    }
}

static std::vector<double> random_vector(size_t n, unsigned seed)
{
    std::mt19937_64 gen(seed);                           // random generator
    std::uniform_real_distribution<double> U(-1.0, 1.0); // uniform in [-1, 1]
    std::vector<double> v(n);                            // result
    for (double& e : v) {
        e = U(gen);
    }
    return v;
}

// ---- dense reference ----

// Dense n x n row-major matrix assembled from the blocks with the rules of blockmat.h (blocks added, so nic = 2 works).
static std::vector<double> dense_from_blocks(const BlockMatrix& M)
{
    const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // matrix size
    std::vector<double> A(n * n, 0.0);                      // dense matrix
    auto add_block = [&](const std::vector<double>& B, size_t p, size_t q) {
        for (int r = 0; r < NEQ; ++r) {
            for (int c = 0; c < NEQ; ++c) {
                A[(p * NEQ + r) * n + q * NEQ + c] += B[p * NB + r * NEQ + c];
            }
        }
    };
    for (int j = 0; j < M.njc; ++j) {
        for (int i = 0; i < M.nic; ++i) {
            const size_t p = static_cast<size_t>(M.index(i, j)); // block row
            add_block(M.D, p, p);
            add_block(M.Wb, p, static_cast<size_t>(M.index(M.west(i), j)));
            add_block(M.Eb, p, static_cast<size_t>(M.index(M.east(i), j)));
            if (j > 0) {
                add_block(M.Sb, p, static_cast<size_t>(M.index(i, j - 1)));
            }
            if (j < M.njc - 1) {
                add_block(M.Nb, p, static_cast<size_t>(M.index(i, j + 1)));
            }
        }
    }
    return A;
}

static void dense_matvec(const std::vector<double>& A, const std::vector<double>& x, std::vector<double>& y)
{
    const size_t n = x.size(); // matrix size
    y.assign(n, 0.0);
    for (size_t r = 0; r < n; ++r) {
        double s = 0.0; // row sum
        const double* a = &A[r * n]; // row r
        for (size_t c = 0; c < n; ++c) {
            s += a[c] * x[c];
        }
        y[r] = s;
    }
}

// Solve A x = b by dense LU with partial pivoting (A is copied).
static std::vector<double> dense_solve(std::vector<double> A, std::vector<double> b)
{
    const size_t n = b.size(); // matrix size
    for (size_t c = 0; c < n; ++c) {
        size_t piv = c; // pivot row
        for (size_t r = c + 1; r < n; ++r) {
            if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) {
                piv = r;
            }
        }
        if (piv != c) {
            for (size_t k = 0; k < n; ++k) {
                std::swap(A[c * n + k], A[piv * n + k]);
            }
            std::swap(b[c], b[piv]);
        }
        const double inv = 1.0 / A[c * n + c]; // 1 / pivot
        for (size_t r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] * inv; // elimination factor
            if (f == 0.0) {
                continue;
            }
            double* ar = &A[r * n];       // row r
            const double* ac = &A[c * n]; // pivot row
            for (size_t k = c; k < n; ++k) {
                ar[k] -= f * ac[k];
            }
            b[r] -= f * b[c];
        }
    }
    std::vector<double> x(n); // solution
    for (size_t r = n; r-- > 0;) {
        double s = b[r]; // right-hand side minus known terms
        for (size_t k = r + 1; k < n; ++k) {
            s -= A[r * n + k] * x[k];
        }
        x[r] = s / A[r * n + r];
    }
    return x;
}

// ||a - b|| / ||b||
static double rel_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double num = 0.0; // ||a - b||^2
    double den = 0.0; // ||b||^2
    for (size_t k = 0; k < a.size(); ++k) {
        num += (a[k] - b[k]) * (a[k] - b[k]);
        den += b[k] * b[k];
    }
    return std::sqrt(num / den);
}

// ||b - M x|| / ||b||
static double rel_res(const BlockMatrix& M, const std::vector<double>& b, const std::vector<double>& x)
{
    std::vector<double> r; // linear residual
    euler::block_residual(M, b, x, r);
    return euler::norm2(r) / euler::norm2(b);
}

// ---- tests ----

static void test_matvec(int nic, int njc)
{
    BlockMatrix M; // test matrix
    make_dominant(M, nic, njc, 11u + static_cast<unsigned>(nic), 3.0, 1.0, 1.0);
    // Fill the unused S (j = 0) and N (j = njc-1) blocks with garbage: they must be ignored.
    for (int i = 0; i < nic; ++i) {
        for (int k = 0; k < NB; ++k) {
            M.Sb[static_cast<size_t>(M.index(i, 0)) * NB + k] = 1e6;
            M.Nb[static_cast<size_t>(M.index(i, njc - 1)) * NB + k] = 1e6;
        }
    }
    const std::vector<double> A = dense_from_blocks(M);                                // dense reference
    const std::vector<double> x = random_vector(static_cast<size_t>(M.ncell()) * NEQ, 5u); // input
    std::vector<double> y;                                                             // block product
    std::vector<double> yd;                                                            // dense product
    euler::block_matvec(M, x, y);
    dense_matvec(A, x, yd);
    const double e = rel_diff(y, yd); // relative difference
    check(e < 1e-14, "A  matvec %3dx%-3d  ||y_block - y_dense|| / ||y_dense|| = %.2e", nic, njc, e);
}

static void test_diagonal()
{
    const int nic = 32; // cells in i
    const int njc = 24; // cells in j
    BlockMatrix M;      // only D is used
    make_dominant(M, nic, njc, 21u, 5.0, 1.0, 0.0);
    std::vector<double> Dinv; // inverse diagonal blocks
    euler::factor_diagonal(M, Dinv);
    const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // vector size
    const std::vector<double> b = random_vector(n, 22u);   // right-hand side
    std::vector<double> x;                                  // D^-1 b
    euler::diagonal_solve(Dinv, b, x);
    double worst = 0.0;  // worst relative block residual ||D_p x_p - b_p|| / ||b_p||
    double worstI = 0.0; // worst |D_p Dinv_p - I| entry
    for (int p = 0; p < M.ncell(); ++p) {
        double num = 0.0; // ||D_p x_p - b_p||^2
        double den = 0.0; // ||b_p||^2
        for (int r = 0; r < NEQ; ++r) {
            double s = 0.0; // (D_p x_p)_r
            for (int c = 0; c < NEQ; ++c) {
                s += M.D[static_cast<size_t>(p) * NB + r * NEQ + c] * x[static_cast<size_t>(p) * NEQ + c];
                double prod = 0.0; // (D_p Dinv_p)_rc
                for (int k = 0; k < NEQ; ++k) {
                    prod += M.D[static_cast<size_t>(p) * NB + r * NEQ + k] * Dinv[static_cast<size_t>(p) * NB + k * NEQ + c];
                }
                worstI = std::max(worstI, std::fabs(prod - (r == c ? 1.0 : 0.0)));
            }
            const double br = b[static_cast<size_t>(p) * NEQ + r]; // b_p component r
            num += (s - br) * (s - br);
            den += br * br;
        }
        worst = std::max(worst, std::sqrt(num / den));
    }
    check(worst < 1e-14, "B1 D x = b, worst block ||D x - b|| / ||b|| = %.2e  (%d blocks)", worst, M.ncell());
    check(worstI < 1e-14, "B2 worst |D Dinv - I| entry = %.2e", worstI);

    // A block that needs pivoting (zero (0,0) entry) and a permutation-like block.
    const double P[NB] = {0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 3, 0, 0, 4, 1}; // needs row swaps
    double Pinv[NB];                                                       // its inverse
    bool ok = euler::invert4(P, Pinv);                                     // inversion succeeded
    double err = 0.0;                                                      // worst |P Pinv - I| entry
    for (int r = 0; r < NEQ; ++r) {
        for (int c = 0; c < NEQ; ++c) {
            double s = 0.0; // (P Pinv)_rc
            for (int k = 0; k < NEQ; ++k) {
                s += P[r * NEQ + k] * Pinv[k * NEQ + c];
            }
            err = std::max(err, std::fabs(s - (r == c ? 1.0 : 0.0)));
        }
    }
    check(ok && err < 1e-15, "B3 block with zero leading entry (pivoting): |P Pinv - I| = %.2e", err);

    // Singular block: row 2 = row 0 + row 1.
    const double S[NB] = {1, 2, 3, 4, 0, 1, 5, 2, 1, 3, 8, 6, 7, 0, 1, 1}; // row 2 = row 0 + row 1
    double Sinv[NB];                                                       // not used
    const bool sing = !euler::invert4(S, Sinv);                            // detected as singular
    BlockMatrix Z;                                                         // all-zero matrix
    Z.resize(3, 2);
    bool threw = false; // factor_diagonal threw on the zero blocks
    try {
        std::vector<double> Zinv; // not used
        euler::factor_diagonal(Z, Zinv);
    }
    catch (const std::runtime_error&) {
        threw = true;
    }
    check(sing && threw, "B4 singular block -> invert4 false (%d), zero blocks -> factor_diagonal throws (%d)",
          sing ? 1 : 0, threw ? 1 : 0);
}

static void test_sgs(int nic, int njc, bool dense_check)
{
    BlockMatrix M; // block-diagonally-dominant test matrix
    make_dominant(M, nic, njc, 31u + static_cast<unsigned>(njc), 8.0, 1.0, 0.4);
    std::vector<double> Dinv; // inverse diagonal blocks
    euler::factor_diagonal(M, Dinv);
    const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // vector size
    const std::vector<double> b = random_vector(n, 32u);   // right-hand side

    for (double omega : {1.0, 0.7}) {
        std::printf("     SGS %dx%d omega = %.1f, ||b - M x|| / ||b|| after k sweeps from 0:", nic, njc, omega);
        double prev = 1.0;    // residual with 0 sweeps (x = 0)
        bool mono = true;     // decreasing at every sweep
        double last = 1.0;    // residual after the last sweep shown
        for (int k = 1; k <= 12; ++k) {
            std::vector<double> x(n, 0.0); // fresh initial guess
            euler::block_sgs(M, Dinv, b, x, k, omega);
            last = rel_res(M, b, x);
            if ((k - 1) % 3 == 0) {
                std::printf(" %d:%.1e", k, last);
            }
            mono = mono && (last < prev || prev < 1e-14); // once at round-off the residual just wobbles
            prev = last;
        }
        std::printf("\n");
        check(mono && last < 1e-6,
              "C1 SGS %dx%d omega = %.1f: monotone decrease over 12 sweeps (down to round-off 1e-14), final %.2e",
              nic, njc, omega, last);
    }

    std::vector<double> x(n, 0.0); // many sweeps
    euler::block_sgs(M, Dinv, b, x, 60);
    const double r60 = rel_res(M, b, x); // residual after 60 sweeps
    if (dense_check) {
        const std::vector<double> xd = dense_solve(dense_from_blocks(M), b); // direct solution
        const double e = rel_diff(x, xd);                                    // distance to it
        check(e < 1e-12, "C2 SGS %dx%d 60 sweeps: residual %.2e, ||x - x_LU|| / ||x_LU|| = %.2e", nic, njc, r60, e);
    }
    else {
        check(r60 < 1e-13, "C2 SGS %dx%d 60 sweeps: residual %.2e", nic, njc, r60);
    }

    // 3 sweeps then 2 more from the result = 5 sweeps from zero (bitwise).
    std::vector<double> x5(n, 0.0); // 5 sweeps at once
    std::vector<double> x32(n, 0.0); // 3 + 2 sweeps
    euler::block_sgs(M, Dinv, b, x5, 5, 0.8);
    euler::block_sgs(M, Dinv, b, x32, 3, 0.8);
    euler::block_sgs(M, Dinv, b, x32, 2, 0.8);
    const double d = rel_diff(x32, x5); // must be exactly 0
    check(d == 0.0, "C3 SGS %dx%d 3 + 2 sweeps = 5 sweeps: difference %.2e", nic, njc, d);
}

static void test_sgs_block_diagonal()
{
    BlockMatrix M; // D only
    make_dominant(M, 7, 5, 41u, 4.0, 1.0, 0.0);
    std::vector<double> Dinv; // inverse diagonal blocks
    euler::factor_diagonal(M, Dinv);
    const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // vector size
    const std::vector<double> b = random_vector(n, 42u);   // right-hand side
    std::vector<double> x(n, 0.0);                         // one sweep
    euler::block_sgs(M, Dinv, b, x, 1);
    const double r = rel_res(M, b, x); // must be round-off
    check(r < 1e-15, "C4 block-diagonal matrix, one SGS sweep exact: residual %.2e", r);
}

static void test_gmres(int nic, int njc, bool dense_check)
{
    BlockMatrix M; // nonsymmetric weakly dominant matrix
    make_convdiff(M, nic, njc, 51u + static_cast<unsigned>(nic), 0.05, 0.3);
    std::vector<double> Dinv; // inverse diagonal blocks
    euler::factor_diagonal(M, Dinv);
    const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // vector size
    const std::vector<double> b = random_vector(n, 52u);   // right-hand side
    std::vector<double> xd;                                 // dense LU solution
    if (dense_check) {
        xd = dense_solve(dense_from_blocks(M), b);
    }

    // Operator and preconditioners as functions.
    const euler::LinearOperator A = [&](const std::vector<double>& v, std::vector<double>& Av) {
        euler::block_matvec(M, v, Av);
    };
    const euler::LinearOperator Id = [](const std::vector<double>& r, std::vector<double>& z) { z = r; };
    const euler::LinearOperator Sgs2 = [&](const std::vector<double>& r, std::vector<double>& z) {
        z.assign(r.size(), 0.0); // SGS from zero: a fixed linear operator
        euler::block_sgs(M, Dinv, r, z, 2);
    };
    // SGS on its own (for comparison): residual after 10 sweeps.
    {
        std::vector<double> xs(n, 0.0); // plain SGS iterate
        euler::block_sgs(M, Dinv, b, xs, 10);
        std::printf("     %dx%d: plain SGS, 10 sweeps -> residual %.2e (weak dominance: slow)\n", nic, njc, rel_res(M, b, xs));
    }

    struct Case {
        const char* name;               // label
        const euler::LinearOperator* P; // preconditioner
        int restart;                    // Krylov size
    };
    const Case cases[] = {{"no precond, restart 30", &Id, 30},
                          {"SGS(2),    restart 30", &Sgs2, 30},
                          {"no precond, restart 10", &Id, 10},
                          {"SGS(2),    restart  5", &Sgs2, 5}};
    int it_none = 0; // iterations without preconditioner (restart 30)
    int it_sgs = 0;  // iterations with SGS (restart 30)
    for (const Case& c : cases) {
        std::vector<double> x(n, 0.0); // solution
        const auto t0 = std::chrono::steady_clock::now();
        const GmresResult g = euler::gmres(A, *c.P, b, x, c.restart, 20000, 1e-10);
        const double t = ms_since(t0);                 // solve time
        const double rtrue = rel_res(M, b, x);         // recomputed residual
        const double e = dense_check ? rel_diff(x, xd) : 0.0; // error vs LU
        const bool ok = g.converged && rtrue <= 1e-10 && std::fabs(rtrue - g.rel_residual) <= 1e-3 * rtrue &&
                        (!dense_check || e < 1e-7);
        check(ok, "D  GMRES %3dx%-3d %s: %5d its, %7.1f ms, residual %.2e (reported %.2e)%s%.2e", nic, njc, c.name,
              g.iterations, t, rtrue, g.rel_residual, dense_check ? ", ||x - x_LU||/||x_LU|| = " : ", -", e);
        if (c.P == &Id && c.restart == 30) {
            it_none = g.iterations;
        }
        if (c.P == &Sgs2 && c.restart == 30) {
            it_sgs = g.iterations;
        }
    }
    check(it_sgs * 3 <= it_none, "D  GMRES %dx%d SGS(2) preconditioning: %d vs %d iterations (ratio %.1f)", nic, njc,
          it_sgs, it_none, static_cast<double>(it_none) / std::max(it_sgs, 1));

    // F. Matrix-free: the same operator through the dense matrix (different code path).
    if (dense_check) {
        const std::vector<double> Ad = dense_from_blocks(M); // dense matrix
        const euler::LinearOperator Adense = [&](const std::vector<double>& v, std::vector<double>& Av) {
            dense_matvec(Ad, v, Av);
        };
        std::vector<double> x1(n, 0.0); // block operator
        std::vector<double> x2(n, 0.0); // dense operator
        const GmresResult g1 = euler::gmres(A, Sgs2, b, x1, 30, 2000, 1e-10);
        const GmresResult g2 = euler::gmres(Adense, Sgs2, b, x2, 30, 2000, 1e-10);
        const double d = rel_diff(x2, x1); // difference of the two solutions
        check(g2.converged && std::abs(g1.iterations - g2.iterations) <= 1 && d < 1e-9,
              "F1 GMRES %dx%d dense-matrix lambda vs block lambda: %d vs %d its, ||x2 - x1||/||x1|| = %.2e", nic, njc,
              g2.iterations, g1.iterations, d);
    }

    /*
       F2. JFNK-style operator (slide 40): J v ~ (R(W + eps v) - R(W)) / eps with R(W) = M W - b,
       eps = sqrt(machine eps) (1 + ||W||) / ||v||. Only residual evaluations, no matrix in GMRES.
    */
    {
        const std::vector<double> W0 = random_vector(n, 53u); // linearisation state
        std::vector<double> R0;                                // R(W0)
        euler::block_residual(M, b, W0, R0);                   // b - M W0 = -R(W0)
        const double nW = euler::norm2(W0);                    // ||W0||
        int nres = 0;                                          // residual evaluations
        const euler::LinearOperator Jfd = [&](const std::vector<double>& v, std::vector<double>& Jv) {
            const double nv = euler::norm2(v); // ||v||
            if (nv == 0.0) {
                Jv.assign(v.size(), 0.0);
                return;
            }
            const double eps = std::sqrt(2.2e-16) * (1.0 + nW) / nv; // finite-difference step
            std::vector<double> Wp(W0);                               // W0 + eps v
            euler::axpy(eps, v, Wp);
            std::vector<double> Rp; // b - M (W0 + eps v)
            euler::block_residual(M, b, Wp, Rp);
            ++nres;
            Jv.resize(v.size());
            for (size_t k = 0; k < v.size(); ++k) {
                Jv[k] = -(Rp[k] - R0[k]) / eps; // R = -(b - M W)
            }
        };
        std::vector<double> xf(n, 0.0); // matrix-free solution
        std::vector<double> xe(n, 0.0); // exact-operator solution
        const GmresResult gf = euler::gmres(Jfd, Sgs2, b, xf, 30, 2000, 1e-6);
        const GmresResult ge = euler::gmres(A, Sgs2, b, xe, 30, 2000, 1e-6);
        const double d = rel_diff(xf, xe);  // difference of the two solutions
        const double r = rel_res(M, b, xf); // true residual of the matrix-free solution
        check(gf.converged && r < 1e-5 && std::abs(gf.iterations - ge.iterations) <= 2,
              "F2 GMRES %dx%d finite-difference J v (JFNK): %d its (exact op %d), %d residual evals, true residual %.2e, "
              "||x_fd - x||/||x|| = %.2e",
              nic, njc, gf.iterations, ge.iterations, nres, r, d);
    }
}

static void test_gmres_special()
{
    const int nic = 7; // cells in i
    const int njc = 5; // cells in j
    BlockMatrix M;     // general test matrix
    make_convdiff(M, nic, njc, 61u, 0.1, 0.3);
    const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // vector size
    const euler::LinearOperator A = [&](const std::vector<double>& v, std::vector<double>& Av) {
        euler::block_matvec(M, v, Av);
    };
    const euler::LinearOperator Id = [](const std::vector<double>& r, std::vector<double>& z) { z = r; };

    // E1: b = 0 -> x = 0 at once (x starts as garbage).
    {
        const std::vector<double> b(n, 0.0); // zero right-hand side
        std::vector<double> x = random_vector(n, 62u); // garbage initial guess
        const GmresResult g = euler::gmres(A, Id, b, x, 20, 100, 1e-10);
        const double nx = euler::norm2(x); // must be 0
        check(g.converged && g.iterations == 0 && nx == 0.0, "E1 b = 0: %d its, ||x|| = %.1e", g.iterations, nx);
    }
    // E2: exact initial guess -> 0 iterations.
    {
        const std::vector<double> xe = random_vector(n, 63u); // chosen solution
        std::vector<double> b;                                 // b = M xe
        euler::block_matvec(M, xe, b);
        std::vector<double> x(xe); // start at the solution
        const GmresResult g = euler::gmres(A, Id, b, x, 20, 100, 1e-10);
        check(g.converged && g.iterations == 0, "E2 exact initial guess: %d its, residual %.2e", g.iterations,
              g.rel_residual);
    }
    // E3: happy breakdown: same 4x4 block in every cell, same b_p everywhere -> Krylov space of dimension <= 4.
    {
        BlockMatrix B; // block-diagonal, identical blocks
        make_dominant(B, nic, njc, 64u, 3.0, 1.0, 0.0);
        for (size_t p = 1; p < static_cast<size_t>(B.ncell()); ++p) {
            for (int k = 0; k < NB; ++k) {
                B.D[p * NB + k] = B.D[k];
            }
        }
        const std::vector<double> b0 = random_vector(NEQ, 65u); // one cell's right-hand side
        std::vector<double> b(n);                                // repeated in every cell
        for (size_t k = 0; k < n; ++k) {
            b[k] = b0[k % NEQ];
        }
        const euler::LinearOperator AB = [&](const std::vector<double>& v, std::vector<double>& Av) {
            euler::block_matvec(B, v, Av);
        };
        std::vector<double> x(n, 0.0); // solution
        const GmresResult g = euler::gmres(AB, Id, b, x, 30, 100, 1e-14);
        const double r = rel_res(B, b, x); // true residual
        const bool finite = std::isfinite(euler::norm2(x)); // no 0/0 from the breakdown
        check(g.converged && g.iterations <= NEQ && r < 1e-13 && finite,
              "E3 happy breakdown (Krylov dim <= 4): %d its, residual %.2e", g.iterations, r);

        // Exact preconditioner (D^-1 on a block-diagonal matrix): one iteration.
        std::vector<double> Dinv; // inverse blocks of B
        euler::factor_diagonal(B, Dinv);
        const euler::LinearOperator PD = [&](const std::vector<double>& rr, std::vector<double>& z) {
            euler::diagonal_solve(Dinv, rr, z);
        };
        std::vector<double> x1(n, 0.0); // solution
        const GmresResult g1 = euler::gmres(AB, PD, b, x1, 30, 100, 1e-14);
        const double r1 = rel_res(B, b, x1); // true residual
        check(g1.converged && g1.iterations == 1 && r1 < 1e-14,
              "E4 exact preconditioner -> 1 iteration: %d its, residual %.2e", g1.iterations, r1);
    }
    // E5: max_iter reached -> converged = false, iterations = max_iter, residual still decreased.
    {
        const std::vector<double> b = random_vector(n, 66u); // right-hand side
        std::vector<double> x(n, 0.0);                        // solution
        const GmresResult g = euler::gmres(A, Id, b, x, 4, 10, 1e-14);
        const double r = rel_res(M, b, x); // true residual
        check(!g.converged && g.iterations == 10 && r < 1.0 && std::fabs(r - g.rel_residual) < 1e-12,
              "E5 max_iter = 10, restart 4: converged = %d, %d its, residual %.2e (reported %.2e)", g.converged ? 1 : 0,
              g.iterations, r, g.rel_residual);
    }
}

// Fast fill for the big timing matrices (LCG, values in [-1, 1)).
static void fill_fast(std::vector<double>& v, std::uint64_t seed, double scale, double diag)
{
    std::uint64_t s = seed; // generator state
    for (size_t k = 0; k < v.size(); ++k) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        const double u = static_cast<double>(s >> 11) * (2.0 / 9007199254740992.0) - 1.0; // in [-1, 1)
        v[k] = scale * u + ((k % NB) % (NEQ + 1) == 0 ? diag : 0.0);
    }
}

static void test_timing(int nic, int njc)
{
    try {
        BlockMatrix M; // big matrix
        auto t0 = std::chrono::steady_clock::now();
        M.resize(nic, njc);
        fill_fast(M.D, 1, 1.0, 10.0);
        fill_fast(M.Wb, 2, 0.4, 0.0);
        fill_fast(M.Eb, 3, 0.4, 0.0);
        fill_fast(M.Sb, 4, 0.4, 0.0);
        fill_fast(M.Nb, 5, 0.4, 0.0);
        const double tfill = ms_since(t0); // allocation + fill time
        const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // vector size
        std::vector<double> b(n);  // right-hand side
        fill_fast(b, 6, 1.0, 0.0);
        std::vector<double> x(n, 0.0); // SGS iterate
        std::vector<double> y(n, 0.0); // matvec output
        std::vector<double> Dinv;      // inverse diagonal blocks

        t0 = std::chrono::steady_clock::now();
        euler::factor_diagonal(M, Dinv);
        const double tfac = ms_since(t0); // factor time

        const int nrep = 3; // repetitions for the averages
        t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < nrep; ++r) {
            euler::block_matvec(M, b, y);
        }
        const double tmv = ms_since(t0) / nrep; // one matvec

        t0 = std::chrono::steady_clock::now();
        euler::block_sgs(M, Dinv, b, x, nrep);
        const double tsgs = ms_since(t0) / nrep; // one sweep (forward + backward)

        const double rr = rel_res(M, b, x); // sanity: the sweeps did reduce the residual
        const double gb = (5.0 * NB * M.ncell() * 8.0) / 1e9; // matrix memory
        check(rr < 0.1, "G  timing %dx%d (%.2f GB matrix, fill %.0f ms): factor_diagonal %.1f ms, matvec %.1f ms, "
                        "SGS sweep (fwd+bwd) %.1f ms, residual after %d sweeps %.1e",
              nic, njc, gb, tfill, tfac, tmv, tsgs, nrep, rr);
    }
    catch (const std::bad_alloc&) {
        std::printf("  [SKIP] timing %dx%d: not enough memory\n", nic, njc);
    }
}

int main()
{
    try {
        std::printf("A. block_matvec vs dense\n");
        test_matvec(7, 5);
        test_matvec(32, 24);
        test_matvec(2, 3); // west = east neighbour: both blocks add on the same cell
        test_matvec(1, 4); // the cell is its own west and east neighbour

        std::printf("B. diagonal blocks\n");
        test_diagonal();

        std::printf("C. block symmetric Gauss-Seidel\n");
        test_sgs(7, 5, true);
        test_sgs(32, 24, true);
        test_sgs_block_diagonal();

        std::printf("D/F. GMRES (right preconditioning) and matrix-free operators\n");
        test_gmres(7, 5, true);
        test_gmres(32, 24, true);

        std::printf("E. GMRES special cases\n");
        test_gmres_special();

        std::printf("G. timing\n");
        test_timing(1024, 1024);
        test_timing(2048, 2048);
    }
    catch (const std::exception& e) {
        std::printf("  [FAIL] exception: %s\n", e.what());
        ++g_fail;
    }
    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail;
}
