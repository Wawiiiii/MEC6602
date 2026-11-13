#include "jacobian.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace euler {

    constexpr int NB = NEQ * NEQ; // number of entries of one 4x4 block

    void flux_jacobian(const double q[NEQ], double nx, double ny, double A[NEQ * NEQ])
    {
        const double g1 = GAMMA - 1.0;                              // gamma - 1
        const double u = q[1];                                      // x velocity
        const double v = q[2];                                      // y velocity
        const double V = u * nx + v * ny;                           // contravariant velocity V = V.n
        const double q2 = u * u + v * v;                            // |V|^2
        const double phi = 0.5 * g1 * q2;                           // phi = (gamma-1) |V|^2 / 2
        const double H = GAMMA / g1 * q[3] / q[0] + 0.5 * q2;       // total enthalpy H

        /*
        dFc/dW with p = (gamma-1) (rho E - ((rho u)^2 + (rho v)^2) / (2 rho)),
        i.e. dp/dW = (gamma-1) (|V|^2/2, -u, -v, 1)  (Blazek, 2D Euler Jacobian):
          row 0: d(rho V)
          row 1: d(rho u V + nx p)
          row 2: d(rho v V + ny p)
          row 3: d(rho H V)
        */
        A[0] = 0.0;
        A[1] = nx;
        A[2] = ny;
        A[3] = 0.0;

        A[4] = nx * phi - u * V;
        A[5] = V + u * nx - g1 * u * nx;
        A[6] = u * ny - g1 * v * nx;
        A[7] = g1 * nx;

        A[8] = ny * phi - v * V;
        A[9] = v * nx - g1 * u * ny;
        A[10] = V + v * ny - g1 * v * ny;
        A[11] = g1 * ny;

        A[12] = V * (phi - H);
        A[13] = H * nx - g1 * u * V;
        A[14] = H * ny - g1 * v * V;
        A[15] = GAMMA * V;
    }

    void roe_abs_matrix(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double absA[NEQ * NEQ],
                        double harten_fraction)
    {
        const double g1 = GAMMA - 1.0; // gamma - 1

        // Total enthalpy on each side (same as roe_flux).
        const double HL = GAMMA / g1 * qL[3] / qL[0] + 0.5 * (qL[1] * qL[1] + qL[2] * qL[2]); // left total enthalpy
        const double HR = GAMMA / g1 * qR[3] / qR[0] + 0.5 * (qR[1] * qR[1] + qR[2] * qR[2]); // right total enthalpy

        // Roe average state (slide 9), copied from roe_flux.
        const double sL = std::sqrt(qL[0]);                       // sqrt of left density
        const double sR = std::sqrt(qR[0]);                       // sqrt of right density
        const double wL = sL / (sL + sR);                         // Roe weight of the left state
        const double wR = sR / (sL + sR);                         // Roe weight of the right state
        const double rho_t = sL * sR;                             // rho~
        const double u_t = wL * qL[1] + wR * qR[1];               // u~
        const double v_t = wL * qL[2] + wR * qR[2];               // v~
        const double H_t = wL * HL + wR * HR;                     // H~
        const double q2_t = u_t * u_t + v_t * v_t;                // q~^2
        const double c_t = std::sqrt(g1 * (H_t - 0.5 * q2_t));    // c~
        const double V_t = u_t * nx + v_t * ny;                   // V~ (contravariant)
        const double c2 = c_t * c_t;                              // c~^2

        // Harten-corrected eigenvalues (slide 11), same as roe_flux.
        const double delta = harten_fraction * c_t;          // entropy fix width
        const double l1 = harten(std::abs(V_t - c_t), delta); // |lambda_1|,     V~ - c~
        const double l234 = harten(std::abs(V_t), delta);     // |lambda_2,3,4|, V~
        const double l5 = harten(std::abs(V_t + c_t), delta); // |lambda_5|,     V~ + c~

        /*
        Column c of |A~| = |A~| e_c (e_c = unit conservative vector).
        The conservative jump e_c is turned into primitive jumps by linearizing at the Roe
        state (Roe's identities, exact for the real jumps W_R - W_L):
          drho = e0,  du = (e1 - u~ e0)/rho~,  dv = (e2 - v~ e0)/rho~,
          dp = (gamma-1) (e3 - u~ e1 - v~ e2 + q~^2/2 e0)
        then the dissipation formula of roe_flux (slide 8, |dF1| + |dF2,3,4| + |dF5|) is applied.
        */
        for (int c = 0; c < NEQ; ++c) {
            double e[NEQ] = {0.0, 0.0, 0.0, 0.0}; // unit conservative jump e_c
            e[c] = 1.0;

            const double drho = e[0];                                                 // density jump
            const double du = (e[1] - u_t * e[0]) / rho_t;                            // x velocity jump
            const double dv = (e[2] - v_t * e[0]) / rho_t;                            // y velocity jump
            const double dp = g1 * (e[3] - u_t * e[1] - v_t * e[2] + 0.5 * q2_t * e[0]); // pressure jump
            const double dV = du * nx + dv * ny;                                      // contravariant velocity jump

            const double s1 = l1 * (dp - rho_t * c_t * dV) / (2.0 * c2); // |lambda_1| x strength of wave 1
            const double s5 = l5 * (dp + rho_t * c_t * dV) / (2.0 * c2); // |lambda_5| x strength of wave 5
            const double s2 = l234 * (drho - dp / c2);                   // |lambda_2| x entropy wave strength
            const double s3 = l234 * rho_t;                              // |lambda_3,4| x rho~ (shear waves)

            absA[0 * NEQ + c] = s1 + s2 + s5;
            absA[1 * NEQ + c] = s1 * (u_t - c_t * nx) + s2 * u_t + s3 * (du - dV * nx) + s5 * (u_t + c_t * nx);
            absA[2 * NEQ + c] = s1 * (v_t - c_t * ny) + s2 * v_t + s3 * (dv - dV * ny) + s5 * (v_t + c_t * ny);
            absA[3 * NEQ + c] = s1 * (H_t - c_t * V_t) + s2 * 0.5 * q2_t + s3 * (u_t * du + v_t * dv - V_t * dV)
                              + s5 * (H_t + c_t * V_t);
        }
    }

    void roe_face_jacobians(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double dFdL[NEQ * NEQ],
                            double dFdR[NEQ * NEQ])
    {
        double AL[NB];   // A(W_L)
        double AR[NB];   // A(W_R)
        double absA[NB]; // |A~| at the Roe state, frozen (slide 26)
        flux_jacobian(qL, nx, ny, AL);
        flux_jacobian(qR, nx, ny, AR);
        roe_abs_matrix(qL, qR, nx, ny, absA);

        // Slide 27: dF/dW_L = 1/2 (A_L + |A~|), dF/dW_R = 1/2 (A_R - |A~|).
        for (int k = 0; k < NB; ++k) {
            dFdL[k] = 0.5 * (AL[k] + absA[k]);
            dFdR[k] = 0.5 * (AR[k] - absA[k]);
        }
    }

    /*
       Exact flux (per unit length) of a physical boundary face as a function of the interior
       cell only, rebuilt like apply_bc + j_face_flux do it:
         wall    : W_g = wall_ghost(W_d)
         farfield: W_g = 2 Wb - W_d, Wb = farfield_state(W_d, free stream, outward normal)
         flux    = Fc(1/2 (W_g + W_d)) through the stored face normal (fnx, fny).
       (onx, ony) = outward normal (farfield) or the wall normal (wall, either direction).
    */
    static void boundary_flux(const double wd[NEQ], bool wall, const Freestream& fs, double onx, double ony,
                              double fnx, double fny, double F[NEQ], const double* w1 = nullptr)
    {
        double qd[NEQ]; // interior cell, primitive
        double wg[NEQ]; // ghost cell, conservative
        cons_to_prim(wd, qd);
        if (wall) {
            double qg[NEQ]; // ghost, primitive
            wall_ghost(qd, onx, ony, qg);
            if (w1) {
                // Extrapolated wall pressure (bc.h / apply_bc): p_g = max(2 p_0 - p_1, 0.2 p_0).
                double q1[NEQ]; // second cell off the wall, primitive
                cons_to_prim(w1, q1);
                qg[3] = std::max(2.0 * qd[3] - q1[3], 0.2 * qd[3]);
            }
            prim_to_cons(qg, wg);
        } else {
            double qb[NEQ]; // boundary state, primitive
            double wb[NEQ]; // boundary state, conservative
            farfield_state(qd, fs.q, onx, ony, qb);
            prim_to_cons(qb, wb);
            for (int k = 0; k < NEQ; ++k) {
                wg[k] = 2.0 * wb[k] - wd[k]; // linear extrapolation to the ghost centre
            }
        }

        double w01[NEQ]; // averaged conservative state at the face
        double q01[NEQ]; // same, primitive
        for (int k = 0; k < NEQ; ++k) {
            w01[k] = 0.5 * (wg[k] + wd[k]);
        }
        cons_to_prim(w01, q01);
        physical_flux(q01, fnx, fny, F);
    }

    /*
       dF_boundary/dW_d by central finite differences of boundary_flux (slide 30: derivative of
       the BC times the flux derivative, done in one go), step h = 1e-7 (|W_c| + rho).
       J[r * NEQ + c] = d F_r / d W_d,c.
    */
    static void boundary_flux_jacobian(const double wd[NEQ], bool wall, const Freestream& fs, double onx, double ony,
                                       double fnx, double fny, double J[NEQ * NEQ])
    {
        for (int c = 0; c < NEQ; ++c) {
            const double h = 1e-7 * (std::abs(wd[c]) + std::abs(wd[0])); // finite-difference step
            double wp[NEQ];                                              // W_d + h e_c
            double wm[NEQ];                                              // W_d - h e_c
            for (int k = 0; k < NEQ; ++k) {
                wp[k] = wd[k];
                wm[k] = wd[k];
            }
            wp[c] += h;
            wm[c] -= h;
            const double h2 = wp[c] - wm[c]; // the step actually taken (round-off safe)

            double Fp[NEQ]; // flux at W_d + h e_c
            double Fm[NEQ]; // flux at W_d - h e_c
            boundary_flux(wp, wall, fs, onx, ony, fnx, fny, Fp);
            boundary_flux(wm, wall, fs, onx, ony, fnx, fny, Fm);
            for (int r = 0; r < NEQ; ++r) {
                J[r * NEQ + c] = (Fp[r] - Fm[r]) / h2;
            }
        }
    }

    /*
       Wall face with the extrapolated ghost pressure: the flux depends on the first two cells off
       the wall, W_0 (= wd) and W_1. J0 = dF/dW_0, J1 = dF/dW_1, central finite differences.
    */
    static void wall_flux_jacobians(const double wd[NEQ], const double w1[NEQ], const Freestream& fs, double fnx,
                                    double fny, double J0[NEQ * NEQ], double J1[NEQ * NEQ])
    {
        for (int which = 0; which < 2; ++which) {
            const double* base = which == 0 ? wd : w1; // the state being perturbed
            double* J = which == 0 ? J0 : J1;
            for (int c = 0; c < NEQ; ++c) {
                const double h = 1e-7 * (std::abs(base[c]) + std::abs(base[0])); // finite-difference step
                double wp[NEQ]; // perturbed +h
                double wm[NEQ]; // perturbed -h
                for (int k = 0; k < NEQ; ++k) {
                    wp[k] = base[k];
                    wm[k] = base[k];
                }
                wp[c] += h;
                wm[c] -= h;
                const double h2 = wp[c] - wm[c]; // the step actually taken
                double Fp[NEQ];
                double Fm[NEQ];
                if (which == 0) {
                    boundary_flux(wp, true, fs, fnx, fny, fnx, fny, Fp, w1);
                    boundary_flux(wm, true, fs, fnx, fny, fnx, fny, Fm, w1);
                } else {
                    boundary_flux(wd, true, fs, fnx, fny, fnx, fny, Fp, wp);
                    boundary_flux(wd, true, fs, fnx, fny, fnx, fny, Fm, wm);
                }
                for (int r = 0; r < NEQ; ++r) {
                    J[r * NEQ + c] = (Fp[r] - Fm[r]) / h2;
                }
            }
        }
    }

    // B += s * J for one 4x4 block (B points to the start of the block).
    static void add_block(double* B, const double J[NEQ * NEQ], double s)
    {
        for (int k = 0; k < NB; ++k) {
            B[k] += s * J[k];
        }
    }

    void assemble_first_order_jacobian(const Mesh& m, const Flow& f, const Freestream& fs, JminBc jmin,
                                       const std::vector<double>& diag_add, BlockMatrix& M)
    {
        const int ncell = m.nic * m.njc; // number of physical cells
        if (!diag_add.empty() && static_cast<int>(diag_add.size()) != ncell) {
            throw std::invalid_argument("assemble_first_order_jacobian: diag_add must have nic * njc entries");
        }
        M.resize(m.nic, m.njc);

        /*
        i-faces, i = 0..nic-1: face i has left cell west(i) and right cell i. Face 0 is the
        TE cut, it couples cell nic-1 (across the cut) with cell 0; the face i = ni-1 of the
        residual is the same physical face, so it is not visited again (each face once).
        Residual signs: R_left += F dS, R_right -= F dS (slide 26).
        */
        // OpenMP: one j row per thread iteration; a row's i-faces only touch blocks of that row.
#pragma omp parallel for schedule(static)
        for (int j = 0; j < m.njc; ++j) {
            double qL[NEQ];  // left state at the face (cell value, first order)
            double qR[NEQ];  // right state at the face
            double dFdL[NB]; // dF/dW_L per unit length
            double dFdR[NB]; // dF/dW_R per unit length
            for (int i = 0; i < m.nic; ++i) {
                const int iw = M.west(i);         // i of the left cell (periodic)
                const int cl = m.cell(iw, j);     // left cell storage index
                const int cr = m.cell(i, j);      // right cell storage index
                const int fc = m.node(i, j);      // face storage index
                const int pL = M.index(iw, j);    // left cell compact index
                const int pR = M.index(i, j);     // right cell compact index
                for (int k = 0; k < NEQ; ++k) {
                    qL[k] = f.Q[k][cl];
                    qR[k] = f.Q[k][cr];
                }
                roe_face_jacobians(qL, qR, m.nix[fc], m.niy[fc], dFdL, dFdR);
                const double S = m.si[fc]; // face length
                add_block(&M.D[pL * NB], dFdL, S);   // dR_L/dW_L += dF/dW_L dS
                add_block(&M.Eb[pL * NB], dFdR, S);  // dR_L/dW_R += dF/dW_R dS (R is east of L)
                add_block(&M.Wb[pR * NB], dFdL, -S); // dR_R/dW_L -= dF/dW_L dS (L is west of R)
                add_block(&M.D[pR * NB], dFdR, -S);  // dR_R/dW_R -= dF/dW_R dS
            }
        }

        /*
        Interior j-faces, j = 1..njc-1: left cell (i, j-1) = south, right cell (i, j) = north.
        OpenMP: one i column per thread iteration (a j-face touches two rows but one column).
        */
#pragma omp parallel for schedule(static)
        for (int i = 0; i < m.nic; ++i) {
            double qL[NEQ];  // left state at the face
            double qR[NEQ];  // right state at the face
            double dFdL[NB]; // dF/dW_L per unit length
            double dFdR[NB]; // dF/dW_R per unit length
            for (int j = 1; j < m.njc; ++j) {
                const int cl = m.cell(i, j - 1);   // left cell storage index
                const int cr = m.cell(i, j);       // right cell storage index
                const int fc = m.node(i, j);       // face storage index
                const int pL = M.index(i, j - 1);  // left cell compact index
                const int pR = M.index(i, j);      // right cell compact index
                for (int k = 0; k < NEQ; ++k) {
                    qL[k] = f.Q[k][cl];
                    qR[k] = f.Q[k][cr];
                }
                roe_face_jacobians(qL, qR, m.njx[fc], m.njy[fc], dFdL, dFdR);
                const double S = m.sj[fc]; // face length
                add_block(&M.D[pL * NB], dFdL, S);   // dR_L/dW_L
                add_block(&M.Nb[pL * NB], dFdR, S);  // dR_L/dW_R (R is north of L)
                add_block(&M.Sb[pR * NB], dFdL, -S); // dR_R/dW_L (L is south of R)
                add_block(&M.D[pR * NB], dFdR, -S);  // dR_R/dW_R
            }
        }

        /*
        Physical boundary faces (slide 30): the ghost is a function of the interior cell, so the
        whole face flux depends on the interior cell only and its derivative goes on D.
          j = 0     : interior cell on the RIGHT -> R_d -= F dS.
                      Wall: wall_ghost with the face normal; Farfield: outward normal -(njx, njy).
          j = nj-1  : farfield, interior cell on the LEFT -> R_d += F dS, outward normal +(njx, njy).
        */
        // OpenMP: one i column per thread iteration (each column has its own boundary cells).
#pragma omp parallel for schedule(static)
        for (int i = 0; i < m.nic; ++i) {
            double wd[NEQ]; // interior cell, conservative
            double J[NB];   // boundary flux Jacobian per unit length

            // jmin
            {
                const int cd = m.cell(i, 0);  // interior cell storage index
                const int fc = m.node(i, 0);  // face storage index
                const int p = M.index(i, 0);  // interior cell compact index
                const double fnx = m.njx[fc]; // face normal (into the fluid)
                const double fny = m.njy[fc];
                for (int k = 0; k < NEQ; ++k) {
                    wd[k] = f.W[k][cd];
                }
                if (jmin == JminBc::Wall && WALL_PRESSURE_EXTRAPOLATE && WALL_EXTRAP_ACTIVE && m.njc > 1) {
                    // Extrapolated wall pressure: the flux also depends on cell (i, 1) -> north block.
                    const int c1 = m.cell(i, 1); // second cell off the wall
                    double w1[NEQ];
                    double J1[NB];
                    for (int k = 0; k < NEQ; ++k) {
                        w1[k] = f.W[k][c1];
                    }
                    wall_flux_jacobians(wd, w1, fs, fnx, fny, J, J1);
                    add_block(&M.Nb[p * NB], J1, -m.sj[fc]);
                } else if (jmin == JminBc::Wall) {
                    boundary_flux_jacobian(wd, true, fs, fnx, fny, fnx, fny, J);
                } else {
                    boundary_flux_jacobian(wd, false, fs, -fnx, -fny, fnx, fny, J);
                }
                add_block(&M.D[p * NB], J, -m.sj[fc]);
            }

            // jmax, farfield
            {
                const int cd = m.cell(i, m.njc - 1);  // interior cell storage index
                const int fc = m.node(i, m.nj - 1);   // face storage index
                const int p = M.index(i, m.njc - 1);  // interior cell compact index
                const double fnx = m.njx[fc];         // face normal (out of the domain)
                const double fny = m.njy[fc];
                for (int k = 0; k < NEQ; ++k) {
                    wd[k] = f.W[k][cd];
                }
                boundary_flux_jacobian(wd, false, fs, fnx, fny, fnx, fny, J);
                add_block(&M.D[p * NB], J, m.sj[fc]);
            }
        }

        // Time term (slide 28): diag_add[p] (e.g. Omega_p / dt_p) on the diagonal of D_p.
        if (!diag_add.empty()) {
#pragma omp parallel for schedule(static)
            for (int p = 0; p < ncell; ++p) {
                for (int k = 0; k < NEQ; ++k) {
                    M.D[p * NB + k * NEQ + k] += diag_add[p];
                }
            }
        }
    }

}
