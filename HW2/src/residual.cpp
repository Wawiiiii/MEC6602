#include "residual.h"

#include "muscl.h"
#include "roe.h"

#include <algorithm>

namespace euler {

    /*
       Left and right primitive states at the face between cells `left` and `right`.
       The four cells lie on one grid line: far_left, left | face | right, far_right.
       First order: the cell values. Second order: MUSCL (muscl.h), with the limiter eps of the
       cell each state belongs to.
    */
    static void reconstruct(const Mesh& m, const Flow& f, Order order, int far_left, int left, int right,
                            int far_right, double q_left[NEQ], double q_right[NEQ])
    {
        if (order == Order::First) {
            for (int k = 0; k < NEQ; ++k) {
                q_left[k] = f.Q[k][left];
                q_right[k] = f.Q[k][right];
            }
            return;
        }

        const double eps_left = limiter_eps(m.vol[left]);   // limiter eps of the left cell
        const double eps_right = limiter_eps(m.vol[right]); // limiter eps of the right cell
        for (int k = 0; k < NEQ; ++k) {
            const std::vector<double>& q = f.Q[k]; // primitive variable k in every cell
            q_left[k] = muscl_left(q[far_left], q[left], q[right], eps_left);
            q_right[k] = muscl_right(q[left], q[right], q[far_right], eps_right);
        }
    }

    void i_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double flux[NEQ])
    {
        const int far_left = m.cell(i - 2, j);  // second cell on the left
        const int left = m.cell(i - 1, j);      // cell on the left of the face
        const int right = m.cell(i, j);         // cell on the right of the face
        const int far_right = m.cell(i + 1, j); // second cell on the right
        const int face = m.node(i, j);          // where the face data is stored

        double q_left[NEQ];  // left state at the face
        double q_right[NEQ]; // right state at the face
        reconstruct(m, f, order, far_left, left, right, far_right, q_left, q_right);
        roe_flux(q_left, q_right, m.nix[face], m.niy[face], flux);
        for (int k = 0; k < NEQ; ++k) {
            flux[k] *= m.si[face];
        }
    }

    void j_face_flux(const Mesh& m, const Flow& f, Order order, int i, int j, double flux[NEQ])
    {
        const int far_left = m.cell(i, j - 2);  // second cell below (towards the wall)
        const int left = m.cell(i, j - 1);      // cell below the face
        const int right = m.cell(i, j);         // cell above the face
        const int far_right = m.cell(i, j + 1); // second cell above (towards the farfield)
        const int face = m.node(i, j);          // where the face data is stored

        /*
        Wall (j = 0) and farfield (j = nj-1): physical flux of the averaged state
        1/2 (W_ghost + W_interior), averaged in conservative variables. The ghost cell holds the
        boundary condition, so at the wall the averaged state has V.n = 0 and at the farfield it is
        exactly the boundary state (bc.cpp).
        */
        if (j == 0 || j == m.nj - 1) {
            double w_face[NEQ]; // averaged state, conservative variables
            double q_face[NEQ]; // same state, primitive variables
            for (int k = 0; k < NEQ; ++k) {
                w_face[k] = 0.5 * (f.W[k][left] + f.W[k][right]);
            }
            cons_to_prim(w_face, q_face);
            physical_flux(q_face, m.njx[face], m.njy[face], flux);
            for (int k = 0; k < NEQ; ++k) {
                flux[k] *= m.sj[face];
            }
            return;
        }

        double q_left[NEQ];  // left state at the face
        double q_right[NEQ]; // right state at the face
        reconstruct(m, f, order, far_left, left, right, far_right, q_left, q_right);
        roe_flux(q_left, q_right, m.njx[face], m.njy[face], flux);
        for (int k = 0; k < NEQ; ++k) {
            flux[k] *= m.sj[face];
        }
    }

    void compute_residual(const Mesh& m, Flow& f, Order order)
    {
        compute_primitives(f);

        // Face flux arrays: one slot per mesh node (allocate_flow sizes them; done here too in case).
        const size_t n_slots = static_cast<size_t>(m.ni_tot) * m.nj_tot; // number of stored nodes
        if (f.Fi[0].size() != n_slots) {
            for (int k = 0; k < NEQ; ++k) {
                f.Fi[k].assign(n_slots, 0.0);
                f.Fj[k].assign(n_slots, 0.0);
            }
        }

        // Step 1: the flux of every face, computed once. One j row per thread.
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            double flux[NEQ]; // Fc dS of the current face
            for (int i = 0; i < m.ni; ++i) {
                i_face_flux(m, f, order, i, j, flux);
                const int face = m.node(i, j);
                for (int k = 0; k < NEQ; ++k) {
                    f.Fi[k][face] = flux[k];
                }
            }
        }
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.nj; ++j) {
            double flux[NEQ]; // Fc dS of the current face
            for (int i = 0; i < m.nic; ++i) {
                j_face_flux(m, f, order, i, j, flux);
                const int face = m.node(i, j);
                for (int k = 0; k < NEQ; ++k) {
                    f.Fj[k][face] = flux[k];
                }
            }
        }

        /*
        Step 2: each cell adds its 4 faces, with the normals pointing out of the cell:
          R(i, j) = F_i(i+1, j) - F_i(i, j) + F_j(i, j+1) - F_j(i, j)
        (what leaves the left cell of a face enters its right cell). Every cell always adds the same
        4 numbers in the same order, so the result does not depend on the number of threads.
        */
        for (int k = 0; k < NEQ; ++k) {
            std::fill(f.R[k].begin(), f.R[k].end(), 0.0);
        }
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int cell = m.cell(i, j);          // the cell
                const int west_south = m.node(i, j);    // its west i-face and its south j-face
                const int east = m.node(i + 1, j);      // its east i-face
                const int north = m.node(i, j + 1);     // its north j-face
                for (int k = 0; k < NEQ; ++k) {
                    f.R[k][cell] = f.Fi[k][east] - f.Fi[k][west_south] + f.Fj[k][north] - f.Fj[k][west_south];
                }
            }
        }
    }

}
