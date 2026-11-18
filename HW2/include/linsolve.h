#pragma once

#include "blockmat.h"

#include <functional>
#include <vector>

namespace euler {



    // Operator y = A(x) for GMRES
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

    // ---- diagonal blocks (4x4 inverses) ----

    /*
       Inverse of one 4x4 row-major block by Gauss-Jordan elimination with partial pivoting.
       claude's work
    */
    bool invert4(const double A[NEQ * NEQ], double Ainv[NEQ * NEQ]);

    /*
       Factorisation of the diagonal blocks, stored as the EXPLICIT INVERSE of every D_p

    */
    void factor_diagonal(const BlockMatrix& M, std::vector<double>& Dinv);

    // Point-implicit / block-Jacobi solve x_p = D_p^-1 b_p for every cell; x resized to the size of b.
    void diagonal_solve(const std::vector<double>& Dinv, const std::vector<double>& b, std::vector<double>& x);

    // ---- block symmetric Gauss-Seidel

    void block_sgs(const BlockMatrix& M, const std::vector<double>& Dinv, const std::vector<double>& b,
                   std::vector<double>& x, int sweeps, double omega = 1.0);

    // ---- restarted GMRES with right preconditioning ----

    struct GmresResult {
        int iterations = 0;        // total number of inner (Arnoldi) iterations = operator applications in the Krylov loop
        double rel_residual = 0.0; // true ||b - A x|| / ||b|| at return (0 if b = 0)
        bool converged = false;    // rel_residual <= rtol
    };

    GmresResult gmres(const LinearOperator& A, const LinearOperator& Minv, const std::vector<double>& b,
                      std::vector<double>& x, int restart, int max_iter, double rtol);

}
