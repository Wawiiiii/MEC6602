#include "linsolve.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace euler {

    constexpr int NB = NEQ * NEQ; // number of entries in one 4x4 block

    // ---- small vector helpers ----

    double dot(const std::vector<double>& a, const std::vector<double>& b)
    {
        const size_t n = a.size(); // vector length
        double s = 0.0;            // running sum
        for (size_t k = 0; k < n; ++k) {
            s += a[k] * b[k];
        }
        return s;
    }

    double norm2(const std::vector<double>& a)
    {
        return std::sqrt(dot(a, a));
    }

    void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y)
    {
        const size_t n = x.size(); // vector length
        for (size_t k = 0; k < n; ++k) {
            y[k] += alpha * x[k];
        }
    }

    // ---- 4x4 block kernels (row-major blocks, NEQ = 4) ----

    // y = y + A x for one block.
    static inline void block_add(const double* A, const double* x, double* y)
    {
        y[0] += A[0] * x[0] + A[1] * x[1] + A[2] * x[2] + A[3] * x[3];
        y[1] += A[4] * x[0] + A[5] * x[1] + A[6] * x[2] + A[7] * x[3];
        y[2] += A[8] * x[0] + A[9] * x[1] + A[10] * x[2] + A[11] * x[3];
        y[3] += A[12] * x[0] + A[13] * x[1] + A[14] * x[2] + A[15] * x[3];
    }

    // y = y - A x for one block.
    static inline void block_sub(const double* A, const double* x, double* y)
    {
        y[0] -= A[0] * x[0] + A[1] * x[1] + A[2] * x[2] + A[3] * x[3];
        y[1] -= A[4] * x[0] + A[5] * x[1] + A[6] * x[2] + A[7] * x[3];
        y[2] -= A[8] * x[0] + A[9] * x[1] + A[10] * x[2] + A[11] * x[3];
        y[3] -= A[12] * x[0] + A[13] * x[1] + A[14] * x[2] + A[15] * x[3];
    }

    // y = A x for one block (y must not alias x).
    static inline void block_mul(const double* A, const double* x, double* y)
    {
        y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2] + A[3] * x[3];
        y[1] = A[4] * x[0] + A[5] * x[1] + A[6] * x[2] + A[7] * x[3];
        y[2] = A[8] * x[0] + A[9] * x[1] + A[10] * x[2] + A[11] * x[3];
        y[3] = A[12] * x[0] + A[13] * x[1] + A[14] * x[2] + A[15] * x[3];
    }

    // Throw if a vector does not have the size NEQ * nic * njc of the matrix.
    static void check_size(const BlockMatrix& M, const std::vector<double>& v, const char* what)
    {
        const size_t n = static_cast<size_t>(M.ncell()) * NEQ; // expected length
        if (v.size() != n) {
            throw std::runtime_error(std::string("linsolve: wrong size of ") + what);
        }
    }

    // ---- block operator ----

    void block_matvec(const BlockMatrix& M, const std::vector<double>& x, std::vector<double>& y)
    {
        check_size(M, x, "x in block_matvec");
        if (&x == &y) {
            throw std::runtime_error("linsolve: block_matvec needs y different from x");
        }
        y.resize(x.size());
        const int nic = M.nic; // cells in i
        const int njc = M.njc; // cells in j
        // Each output row is independent: OpenMP parallel loop over the j rows.
#pragma omp parallel for schedule(static)
        for (int j = 0; j < njc; ++j) {
            for (int i = 0; i < nic; ++i) {
                const size_t p = static_cast<size_t>(M.index(i, j));        // this cell
                const size_t pw = static_cast<size_t>(M.index(M.west(i), j)); // west neighbour (periodic)
                const size_t pe = static_cast<size_t>(M.index(M.east(i), j)); // east neighbour (periodic)
                double* yp = &y[p * NEQ];                                    // output of row p
                block_mul(&M.D[p * NB], &x[p * NEQ], yp);
                block_add(&M.Wb[p * NB], &x[pw * NEQ], yp);
                block_add(&M.Eb[p * NB], &x[pe * NEQ], yp);
                if (j > 0) {
                    block_add(&M.Sb[p * NB], &x[(p - nic) * NEQ], yp); // south (i, j-1); none at the wall
                }
                if (j < njc - 1) {
                    block_add(&M.Nb[p * NB], &x[(p + nic) * NEQ], yp); // north (i, j+1); none at the farfield
                }
            }
        }
    }

    void block_residual(const BlockMatrix& M, const std::vector<double>& b, const std::vector<double>& x,
                        std::vector<double>& r)
    {
        check_size(M, b, "b in block_residual");
        block_matvec(M, x, r);
        const size_t n = b.size(); // vector length
        for (size_t k = 0; k < n; ++k) {
            r[k] = b[k] - r[k];
        }
    }

    // ---- diagonal blocks ----

    bool invert4(const double A[NEQ * NEQ], double Ainv[NEQ * NEQ])
    {
        /*
           Gauss-Jordan on the augmented matrix [a | Ainv], a = copy of A, Ainv starting as the identity.
           Column c: pick the row with the largest |a(r, c)| for r >= c (partial pivoting), swap it up,
           scale it to a unit pivot, then eliminate column c from all other rows. At the end a = I, Ainv = A^-1.
        */
        double a[NB];       // working copy of A
        double amax = 0.0;  // largest |entry| of A (scale for the singularity test)
        for (int k = 0; k < NB; ++k) {
            a[k] = A[k];
            Ainv[k] = 0.0;
            amax = std::max(amax, std::fabs(A[k]));
        }
        const double tiny = 1e-14 * amax; // a pivot this small means numerically singular (condition ~ 1e14)
        for (int k = 0; k < NEQ; ++k) {
            Ainv[k * NEQ + k] = 1.0;
        }
        for (int c = 0; c < NEQ; ++c) {
            int piv = c;                          // pivot row
            double big = std::fabs(a[c * NEQ + c]); // largest |entry| found in column c
            for (int r = c + 1; r < NEQ; ++r) {
                const double v = std::fabs(a[r * NEQ + c]); // candidate magnitude
                if (v > big) {
                    big = v;
                    piv = r;
                }
            }
            if (!(big > tiny) || !std::isfinite(big)) {
                return false; // singular, numerically singular, or NaN / Inf block
            }
            if (piv != c) {
                for (int k = 0; k < NEQ; ++k) {
                    std::swap(a[c * NEQ + k], a[piv * NEQ + k]);
                    std::swap(Ainv[c * NEQ + k], Ainv[piv * NEQ + k]);
                }
            }
            const double s = 1.0 / a[c * NEQ + c]; // scale of the pivot row
            for (int k = 0; k < NEQ; ++k) {
                a[c * NEQ + k] *= s;
                Ainv[c * NEQ + k] *= s;
            }
            for (int r = 0; r < NEQ; ++r) {
                if (r == c) {
                    continue;
                }
                const double f = a[r * NEQ + c]; // multiple of the pivot row removed from row r
                if (f != 0.0) {
                    for (int k = 0; k < NEQ; ++k) {
                        a[r * NEQ + k] -= f * a[c * NEQ + k];
                        Ainv[r * NEQ + k] -= f * Ainv[c * NEQ + k];
                    }
                }
            }
        }
        return true;
    }

    void factor_diagonal(const BlockMatrix& M, std::vector<double>& Dinv)
    {
        const int nc = M.ncell(); // number of cells
        Dinv.resize(static_cast<size_t>(nc) * NB);
        /*
        Blocks are independent: OpenMP parallel loop. An exception cannot leave an OpenMP region,
        so each block records success and the first failure is reported after the loop.
        */
        std::vector<char> ok(static_cast<size_t>(nc), 1); // 1 = block inverted
#pragma omp parallel for schedule(static)
        for (int p = 0; p < nc; ++p) {
            const size_t off = static_cast<size_t>(p) * NB; // start of block p
            ok[p] = invert4(&M.D[off], &Dinv[off]) ? 1 : 0;
        }
        for (int p = 0; p < nc; ++p) {
            if (!ok[p]) {
                throw std::runtime_error("linsolve: singular diagonal block at cell i = " +
                                         std::to_string(p % M.nic) + ", j = " + std::to_string(p / M.nic));
            }
        }
    }

    void diagonal_solve(const std::vector<double>& Dinv, const std::vector<double>& b, std::vector<double>& x)
    {
        if (Dinv.size() * NEQ != b.size() * NB) {
            throw std::runtime_error("linsolve: diagonal_solve sizes of Dinv and b do not match");
        }
        if (&x == &b) {
            throw std::runtime_error("linsolve: diagonal_solve needs x different from b");
        }
        x.resize(b.size());
        const size_t nc = b.size() / NEQ; // number of cells
        for (size_t p = 0; p < nc; ++p) {
            block_mul(&Dinv[p * NB], &b[p * NEQ], &x[p * NEQ]);
        }
    }

    // ---- block symmetric Gauss-Seidel ----

    /*
       One Gauss-Seidel visit of cell (i, j) (slide 34 with 4x4 blocks):
         r = b_p - Wb x_west - Eb x_east - Sb x_south - Nb x_north   (latest available x)
         x_p <- x_p + omega (D_p^-1 r - x_p)
    */
    static inline void sgs_visit(const BlockMatrix& M, const std::vector<double>& Dinv,
                                 const std::vector<double>& b, std::vector<double>& x,
                                 int i, int j, double omega)
    {
        const int nic = M.nic;                                         // cells in i
        const size_t p = static_cast<size_t>(M.index(i, j));          // this cell
        const size_t pw = static_cast<size_t>(M.index(M.west(i), j)); // west neighbour (periodic)
        const size_t pe = static_cast<size_t>(M.index(M.east(i), j)); // east neighbour (periodic)
        double r[NEQ];                                                // right-hand side of the cell
        double xs[NEQ];                                               // Gauss-Seidel value D^-1 r
        for (int k = 0; k < NEQ; ++k) {
            r[k] = b[p * NEQ + k];
        }
        block_sub(&M.Wb[p * NB], &x[pw * NEQ], r);
        block_sub(&M.Eb[p * NB], &x[pe * NEQ], r);
        if (j > 0) {
            block_sub(&M.Sb[p * NB], &x[(p - nic) * NEQ], r);
        }
        if (j < M.njc - 1) {
            block_sub(&M.Nb[p * NB], &x[(p + nic) * NEQ], r);
        }
        block_mul(&Dinv[p * NB], r, xs);
        double* xp = &x[p * NEQ]; // unknowns of this cell
        if (omega == 1.0) {
            for (int k = 0; k < NEQ; ++k) {
                xp[k] = xs[k];
            }
        }
        else {
            for (int k = 0; k < NEQ; ++k) {
                xp[k] += omega * (xs[k] - xp[k]);
            }
        }
    }

    void block_sgs(const BlockMatrix& M, const std::vector<double>& Dinv, const std::vector<double>& b,
                   std::vector<double>& x, int sweeps, double omega)
    {
        check_size(M, b, "b in block_sgs");
        check_size(M, x, "x in block_sgs");
        if (Dinv.size() != static_cast<size_t>(M.ncell()) * NB) {
            throw std::runtime_error("linsolve: block_sgs needs Dinv from factor_diagonal");
        }
        const int nic = M.nic; // cells in i
        const int njc = M.njc; // cells in j
        for (int s = 0; s < sweeps; ++s) {
            // Forward pass (slide 35: lower blocks at k+1, upper blocks at k): i fastest, then j.
            for (int j = 0; j < njc; ++j) {
                for (int i = 0; i < nic; ++i) {
                    sgs_visit(M, Dinv, b, x, i, j, omega);
                }
            }
            // Backward pass: reverse order, the roles of lower and upper blocks swap.
            for (int j = njc - 1; j >= 0; --j) {
                for (int i = nic - 1; i >= 0; --i) {
                    sgs_visit(M, Dinv, b, x, i, j, omega);
                }
            }
        }
    }

    // ---- restarted GMRES with right preconditioning ----

    GmresResult gmres(const LinearOperator& A, const LinearOperator& Minv, const std::vector<double>& b,
                      std::vector<double>& x, int restart, int max_iter, double rtol)
    {
        GmresResult res;               // returned statistics
        const size_t n = b.size();     // system size
        const int m = std::max(restart, 1); // Krylov basis size per cycle
        if (x.size() != n) {
            x.assign(n, 0.0);
        }

        // b = 0: the solution is x = 0.
        const double bnorm = norm2(b); // ||b||
        if (bnorm == 0.0) {
            std::fill(x.begin(), x.end(), 0.0);
            res.converged = true;
            return res;
        }

        /*
           Storage allocated once for the whole call:
             V      m+1 Krylov vectors v_1 .. v_{m+1} (slide 39)
             H      (m+1) x m Hessenberg matrix, column-major: H(i, j) = H[i + j (m+1)]
             cs, sn Givens rotations that make H upper triangular
             g      rotated right-hand side beta e_1; |g(k)| = residual norm after k iterations
        */
        std::vector<std::vector<double>> V(static_cast<size_t>(m) + 1, std::vector<double>(n, 0.0));
        std::vector<double> H(static_cast<size_t>(m + 1) * m, 0.0); // Hessenberg matrix
        std::vector<double> cs(m, 0.0);                             // Givens cosines
        std::vector<double> sn(m, 0.0);                             // Givens sines
        std::vector<double> g(static_cast<size_t>(m) + 1, 0.0);     // rotated beta e_1
        std::vector<double> y(m, 0.0);                              // least-squares solution
        std::vector<double> w(n, 0.0);                              // work: A Minv v_j, then residual
        std::vector<double> z(n, 0.0);                              // work: Minv v_j, then the correction
        auto h = [&](int i, int j) -> double& { return H[static_cast<size_t>(i) + static_cast<size_t>(j) * (m + 1)]; };

        int total = 0; // inner iterations done
        for (;;) {
            // True residual r0 = b - A x (stored in w) and convergence test.
            A(x, w);
            for (size_t k = 0; k < n; ++k) {
                w[k] = b[k] - w[k];
            }
            const double beta = norm2(w); // ||r0||
            res.iterations = total;
            res.rel_residual = beta / bnorm;
            res.converged = res.rel_residual <= rtol;
            if (res.converged || total >= max_iter || !std::isfinite(beta)) {
                return res;
            }

            // v_1 = r0 / beta, g = beta e_1 (slide 39).
            for (size_t k = 0; k < n; ++k) {
                V[0][k] = w[k] / beta;
            }
            std::fill(g.begin(), g.end(), 0.0);
            g[0] = beta;

            int kdim = 0; // Krylov dimension reached in this cycle
            for (int j = 0; j < m && total < max_iter; ++j) {
                // w = A Minv v_j (right preconditioning).
                Minv(V[j], z);
                A(z, w);
                const double wnorm0 = norm2(w); // ||w|| before orthogonalisation (breakdown reference)

                // Modified Gram-Schmidt: h_ij = w . v_i, w <- w - h_ij v_i.
                for (int i = 0; i <= j; ++i) {
                    h(i, j) = dot(w, V[i]);
                    axpy(-h(i, j), V[i], w);
                }
                const double hnext = norm2(w); // h_{j+1,j}
                h(j + 1, j) = hnext;

                // Apply the previous rotations to the new column, then build the one that zeroes h_{j+1,j}.
                for (int i = 0; i < j; ++i) {
                    const double t = cs[i] * h(i, j) + sn[i] * h(i + 1, j); // rotated entry i
                    h(i + 1, j) = -sn[i] * h(i, j) + cs[i] * h(i + 1, j);
                    h(i, j) = t;
                }
                const double den = std::hypot(h(j, j), h(j + 1, j)); // length of the pair to rotate
                if (den > 0.0) {
                    cs[j] = h(j, j) / den;
                    sn[j] = h(j + 1, j) / den;
                }
                else {
                    cs[j] = 1.0;
                    sn[j] = 0.0;
                }
                h(j, j) = cs[j] * h(j, j) + sn[j] * h(j + 1, j);
                h(j + 1, j) = 0.0;
                g[j + 1] = -sn[j] * g[j];
                g[j] = cs[j] * g[j];

                ++total;
                kdim = j + 1;

                // Stop the cycle on convergence (|g_{j+1}| = residual estimate) or happy breakdown.
                const bool breakdown = hnext <= 1e-14 * wnorm0;
                if (breakdown || std::fabs(g[j + 1]) <= rtol * bnorm) {
                    break;
                }
                for (size_t k = 0; k < n; ++k) {
                    V[j + 1][k] = w[k] / hnext; // v_{j+1} = w / h_{j+1,j}
                }
            }

            // y = argmin ||beta e_1 - H y||: back substitution on the rotated (upper triangular) H.
            for (int i = kdim - 1; i >= 0; --i) {
                double s = g[i]; // right-hand side minus the known terms
                for (int k = i + 1; k < kdim; ++k) {
                    s -= h(i, k) * y[k];
                }
                y[i] = (h(i, i) != 0.0) ? s / h(i, i) : 0.0;
            }

            // x = x + Minv (V y).
            std::fill(w.begin(), w.end(), 0.0);
            for (int i = 0; i < kdim; ++i) {
                axpy(y[i], V[i], w);
            }
            Minv(w, z);
            axpy(1.0, z, x);
        }
    }

}
