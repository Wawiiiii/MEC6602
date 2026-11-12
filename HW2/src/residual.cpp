#include "residual.h"

#include "muscl.h"
#include "roe.h"

#include <algorithm>

namespace euler {

    /*
       Left and right primitive states at the face between cells c0 and cp1 (slide 12 step 1).
       The four cells are in a line: cm1, c0 | face | cp1, cp2.
       First order : qL = Q(c0), qR = Q(cp1).
       Second order: MUSCL + Hemker-Koren on each primitive variable (slides 15-16),
                     eps = limiter_eps(area) of the cell the state belongs to (muscl.h).
    */
    static void reconstruct(const Mesh& m, const Flow& f, Order order, int cm1, int c0, int cp1, int cp2,
                            double qL[NEQ], double qR[NEQ])
    {
        if (order == Order::First) {
            for (int k = 0; k < NEQ; ++k) {
                qL[k] = f.Q[k][c0];
                qR[k] = f.Q[k][cp1];
            }
            return;
        }

        const double epsL = limiter_eps(m.vol[c0]);  // limiter eps of the left cell
        const double epsR = limiter_eps(m.vol[cp1]); // limiter eps of the right cell
        for (int k = 0; k < NEQ; ++k) {
            const std::vector<double>& q = f.Q[k]; // primitive variable k
            qL[k] = muscl_left(q[cm1], q[c0], q[cp1], epsL);
            qR[k] = muscl_right(q[c0], q[cp1], q[cp2], epsR);
        }
    }

    void i_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double F[NEQ])
    {
        const int cm1 = m.cell(i - 2, j); // second cell on the left
        const int c0 = m.cell(i - 1, j);  // left cell
        const int cp1 = m.cell(i, j);     // right cell
        const int cp2 = m.cell(i + 1, j); // second cell on the right
        const int fc = m.node(i, j);      // face storage index

        double qL[NEQ]; // left state at the face
        double qR[NEQ]; // right state at the face
        reconstruct(m, f, order, cm1, c0, cp1, cp2, qL, qR);
        roe_flux(qL, qR, m.nix[fc], m.niy[fc], F);
        for (int k = 0; k < NEQ; ++k) {
            F[k] *= m.si[fc];
        }
    }

    void j_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double F[NEQ])
    {
        const int cm1 = m.cell(i, j - 2); // second cell below (towards the wall)
        const int c0 = m.cell(i, j - 1);  // left cell
        const int cp1 = m.cell(i, j);     // right cell
        const int cp2 = m.cell(i, j + 1); // second cell above (towards the farfield)
        const int fc = m.node(i, j);      // face storage index

        /*
        Physical boundary faces (wall j = 0, farfield j = nj-1): flux of the averaged state,
        Fc(W01) with W01 = 1/2 (W_ghost + W_interior) in conservative variables
        (basic_algo slide 7, method #2, "to be used at physical boundaries"; slide 13).
        The ghost cell carries the boundary condition, so at the wall W01 has V.n = 0 and
        at the farfield W01 = Wb exactly (see apply_bc).
        */
        if (j == 0 || j == m.nj - 1) {
            double w01[NEQ]; // averaged conservative state at the face
            double q01[NEQ]; // same, primitive
            for (int k = 0; k < NEQ; ++k) {
                w01[k] = 0.5 * (f.W[k][c0] + f.W[k][cp1]);
            }
            cons_to_prim(w01, q01);
            physical_flux(q01, m.njx[fc], m.njy[fc], F);
            for (int k = 0; k < NEQ; ++k) {
                F[k] *= m.sj[fc];
            }
            return;
        }

        double qL[NEQ]; // left state at the face
        double qR[NEQ]; // right state at the face
        reconstruct(m, f, order, cm1, c0, cp1, cp2, qL, qR);
        roe_flux(qL, qR, m.njx[fc], m.njy[fc], F);
        for (int k = 0; k < NEQ; ++k) {
            F[k] *= m.sj[fc];
        }
    }

    void compute_residual(const Mesh& m, Flow& f, Order order)
    {
        compute_primitives(f);
        for (int k = 0; k < NEQ; ++k) {
            std::fill(f.R[k].begin(), f.R[k].end(), 0.0);
        }

        double F[NEQ]; // (Fc dS) of the current face

        // i-faces: what leaves the left cell (i-1, j) enters the right cell (i, j).
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.ni; ++i) {
                i_face_flux(m, f, order, i, j, F);
                const int cl = m.cell(i - 1, j); // left cell
                const int cr = m.cell(i, j);     // right cell
                for (int k = 0; k < NEQ; ++k) {
                    f.R[k][cl] += F[k];
                    f.R[k][cr] -= F[k];
                }
            }
        }

        // j-faces: left cell (i, j-1), right cell (i, j).
        for (int j = 0; j < m.nj; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                j_face_flux(m, f, order, i, j, F);
                const int cl = m.cell(i, j - 1); // left cell
                const int cr = m.cell(i, j);     // right cell
                for (int k = 0; k < NEQ; ++k) {
                    f.R[k][cl] += F[k];
                    f.R[k][cr] -= F[k];
                }
            }
        }
    }

}
