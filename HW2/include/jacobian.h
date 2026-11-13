#pragma once

#include "bc.h"
#include "blockmat.h"
#include "flow.h"
#include "mesh.h"
#include "roe.h"

#include <vector>

namespace euler {

    /*
       Jacobian of the physical convective flux Fc.n (per unit face length, as physical_flux)
       with respect to the conservative variables W = (rho, rho u, rho v, rho E)
       (Roe_scheme_accel slide 26, the matrices A(W_i) = dF(W_i)/dW_i; Blazek ch. 4, 2D form):
         A = dFc/dW, row-major, A[r * NEQ + c] = d Fc_r / d W_c
       q = (rho, u, v, p) primitive state, (nx, ny) unit normal.
       The Euler flux is homogeneous of degree 1: A(W) W = Fc(W).
    */
    void flux_jacobian(const double q[NEQ], double nx, double ny, double A[NEQ * NEQ]);

    /*
       Roe matrix |A~| = T~ |Lambda~| T~^-1 at the Roe average of qL and qR (slides 8-11, 26),
       with the same Harten-corrected eigenvalues |V~ - c~|, |V~|, |V~ + c~| as roe_flux().
       Built column by column: column c is the roe_flux dissipation formula applied to the
       unit conservative jump e_c, so that |A~| (W_R - W_L) is exactly the dissipation of roe_flux.
       Row-major, absA[r * NEQ + c].
    */
    void roe_abs_matrix(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double absA[NEQ * NEQ],
                        double harten_fraction = HARTEN_FRACTION);

    /*
       First-order face Jacobians of the Roe flux with |A~| frozen (slides 26-27,
       "we assume the matrices |A| do not depend on W"), per unit face length:
         F = 1/2 [Fc(W_L) + Fc(W_R)] - 1/2 |A~| (W_R - W_L)
         dF/dW_L = 1/2 (A(W_L) + |A~|)
         dF/dW_R = 1/2 (A(W_R) - |A~|)
    */
    void roe_face_jacobians(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double dFdL[NEQ * NEQ],
                            double dFdR[NEQ * NEQ]);

    /*
       Implicit operator of the FIRST-ORDER residual on the physical cells (slides 26-30):
         M = dR1/dW + diag(diag_add)
       R1 is compute_residual(Order::First): cell values as face states, Roe on the interior
       faces and on the i-cut, averaged-state flux on the wall / farfield faces.
         - interior and i-cut faces: frozen-|A~| Jacobians of roe_face_jacobians() (slide 27),
           added with the residual signs (R_left += F dS, R_right -= F dS);
         - wall (j = 0) and farfield (j = nj-1) faces: the flux depends only on the interior cell
           (the ghost is a function of it), so dF/dW_interior (finite differences of the exact
           boundary flux) is lumped on the diagonal block D (slide 30, chain rule).
       diag_add: one scalar per physical cell (size nic * njc, e.g. Omega / dt), added to the
       diagonal of D_p; an empty vector adds nothing.
       f.W and f.Q must be up to date with the ghosts filled (apply_bc + compute_primitives).
       M is resized to nic x njc (blockmat.h layout).
    */
    void assemble_first_order_jacobian(const Mesh& m, const Flow& f, const Freestream& fs, JminBc jmin,
                                       const std::vector<double>& diag_add, BlockMatrix& M);

}
