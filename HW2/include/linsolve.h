#pragma once

#include "blockmat.h"

#include <functional>
#include <vector>

namespace euler {

    /*
       Linear solvers for the implicit system  M dW = -R  (Roe_scheme_accel slides 32-40):
         - block matrix-vector product y = M x
         - inverse of the 4x4 diagonal blocks (point implicit, slide 37)
         - block symmetric Gauss-Seidel, forward + backward sweeps (slides 33-35)
         - restarted GMRES with right preconditioning (slide 39), usable matrix-free (JFNK, slide 40)
       Vectors are laid out like BlockMatrix: size NEQ * nic * njc, entry [p * NEQ + k], p = j * nic + i.
    */

    // Operator y = A(x) for GMRES: first argument is the input vector, second the output (gmres passes it already sized).
    using LinearOperator = std::function<void(const std::vector<double>&, std::vector<double>&)>;

    // ---- small vector helpers (plain loops over the whole vector) ----

    double dot(const std::vector<double>& a, const std::vector<double>& b); // sum a_k b_k
    double norm2(const std::vector<double>& a);                             // Euclidean norm sqrt(a . a)
    void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y); // y = y + alpha x

    // ---- block operator ----

    /*
       y = M x, row p (blockmat.h):
         y_p = D_p x_p + Wb_p x_west + Eb_p x_east + Sb_p x_south + Nb_p x_north
       with periodic west / east across the cut and no south (j = 0) / north (j = njc-1) term.
       y is resized to the size of x.
    */
    void block_matvec(const BlockMatrix& M, const std::vector<double>& x, std::vector<double>& y);

    // r = b - M x (the linear residual), r resized to the size of b.
    void block_residual(const BlockMatrix& M, const std::vector<double>& b, const std::vector<double>& x,
                        std::vector<double>& r);

    // ---- diagonal blocks (slide 33: "we need to invert a 4x4 matrix") ----

    /*
       Inverse of one 4x4 row-major block by Gauss-Jordan elimination with partial pivoting.
       Returns false (Ainv left undefined) if a pivot is not finite or below 1e-14 * max|A_ij|
       (singular to double precision; an all-zero block is singular too).
    */
    bool invert4(const double A[NEQ * NEQ], double Ainv[NEQ * NEQ]);

    /*
       Factorisation of the diagonal blocks, stored as the EXPLICIT INVERSE of every D_p:
       Dinv has the layout of M.D (size ncell * 16, block p at p * 16, row-major).
       Reason: the inverse is computed once per implicit step and then applied many times
       (every Gauss-Seidel visit of a cell, every sweep, every GMRES preconditioner call);
       applying it is one 4x4 matrix-vector product (16 multiply-adds, no branches, no pivot array),
       cheaper than a forward + backward LU substitution. The diagonal blocks of the implicit operator
       (V/dt I + dR/dW) are well conditioned, so the round-off of the explicit inverse is harmless.
       Throws std::runtime_error naming the cell if a block is singular.
    */
    void factor_diagonal(const BlockMatrix& M, std::vector<double>& Dinv);

    // Point-implicit / block-Jacobi solve x_p = D_p^-1 b_p for every cell (slide 37); x resized to the size of b.
    void diagonal_solve(const std::vector<double>& Dinv, const std::vector<double>& b, std::vector<double>& x);

    // ---- block symmetric Gauss-Seidel (slides 33-35) ----

    /*
       Block symmetric Gauss-Seidel on M x = b.
       x is the initial guess on input (typically zero) and the result on output.
       One sweep = forward pass over the cells in lexicographic order (i fastest, then j),
       followed by a backward pass in the reverse order. Each cell visit uses the latest values
       of its neighbours (slide 34):
         x_p* = D_p^-1 ( b_p - sum_{q neighbour of p} M_pq x_q )
         x_p  = x_p + omega (x_p* - x_p)                (omega < 1: under-relaxation, slide 33)
       The periodic neighbours across the cut (west of i = 0 is i = nic-1, east of i = nic-1 is i = 0)
       just use their current value. Dinv comes from factor_diagonal(M, Dinv).
       Slide 33: the system is usually not converged; 2 or 3 sweeps are used.
    */
    void block_sgs(const BlockMatrix& M, const std::vector<double>& Dinv, const std::vector<double>& b,
                   std::vector<double>& x, int sweeps, double omega = 1.0);

    // ---- restarted GMRES with right preconditioning (slide 39) ----

    struct GmresResult {
        int iterations = 0;        // total number of inner (Arnoldi) iterations = operator applications in the Krylov loop
        double rel_residual = 0.0; // true ||b - A x|| / ||b|| at return (0 if b = 0)
        bool converged = false;    // rel_residual <= rtol
    };

    /*
       GMRES(m) for A x = b with right preconditioning: solves A Minv u = b, x = x0 + Minv (V y),
       so the minimised residual is the true residual of the unpreconditioned system.
       Arnoldi with modified Gram-Schmidt, least-squares problem by Givens rotations (Saad & Schultz 1986).
         A(v, Av)    : applies the operator (matrix, or matrix-free J v for JFNK, slide 40)
         Minv(r, z)  : applies the preconditioner z ~ M^-1 r (e.g. a few SGS sweeps from zero);
                       it must be a FIXED linear operator (same result for the same r).
         x           : initial guess on input (resized to zeros if its size differs from b), solution on output
         restart     : Krylov basis size m (slide 39: 20 to 50); the basis (m+1 vectors) is allocated once per call
         max_iter    : maximum total number of inner iterations
         rtol        : stop when ||b - A x|| <= rtol ||b||
       b = 0 gives x = 0 at once. A happy breakdown (new Krylov vector of zero norm: the solution is in the
       subspace) ends the cycle early. Inside a cycle the Givens estimate |g_{j+1}| of the residual decides
       when to stop; convergence is then confirmed with the true residual b - A x recomputed at each restart
       (one extra A and one extra Minv per cycle).
    */
    GmresResult gmres(const LinearOperator& A, const LinearOperator& Minv, const std::vector<double>& b,
                      std::vector<double>& x, int restart, int max_iter, double rtol);

}
