#include "bc.h"

#include "forces.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace euler {

    constexpr double PI = 3.14159265358979323846; // MSVC has no M_PI by default

    Freestream make_freestream(double mach, double alpha_deg)
    {
        Freestream fs;
        const double a = alpha_deg * PI / 180.0; // angle of attack in radians
        fs.mach = mach;
        fs.alpha_deg = alpha_deg;
        fs.q[0] = 1.0;               // rho_inf
        fs.q[1] = mach * std::cos(a); // u_inf
        fs.q[2] = mach * std::sin(a); // v_inf
        fs.q[3] = 1.0 / GAMMA;       // p_inf, so c_inf = sqrt(gamma p / rho) = 1
        return fs;
    }

    void init_freestream(Flow& f, const Freestream& fs)
    {
        double w[NEQ]; // conservative free stream
        prim_to_cons(fs.q, w);
        for (int k = 0; k < NEQ; ++k) {
            for (double& value : f.W[k]) {
                value = w[k];
            }
        }
    }

    void wall_ghost(const double qi[NEQ], double nx, double ny, double qg[NEQ])
    {
        const double Vn = qi[1] * nx + qi[2] * ny; // normal velocity of the interior cell
        qg[0] = qi[0];                             // rho_g = rho_i
        qg[1] = qi[1] - 2.0 * Vn * nx;             // velocity mirrored about the wall
        qg[2] = qi[2] - 2.0 * Vn * ny;
        qg[3] = qi[3];                             // p_g = p_i
    }

    void farfield_state(const double qd[NEQ], const double qa[NEQ], double nx, double ny, double qb[NEQ])
    {
        const double rho0 = qd[0];                         // rho0: interior density
        const double c0 = std::sqrt(GAMMA * qd[3] / qd[0]); // c0: interior speed of sound
        const double Vd = qd[1] * nx + qd[2] * ny;         // interior normal velocity (n outward)

        // Supersonic (normal Mach >= 1): all characteristics go the same way.
        if (std::abs(Vd) >= c0) {
            const double* src = Vd < 0.0 ? qa : qd; // inflow: free stream; outflow: interior
            for (int k = 0; k < NEQ; ++k) {
                qb[k] = src[k];
            }
            return;
        }

        const double rc = rho0 * c0; // acoustic impedance rho0 c0
        double pb;                   // boundary pressure
        if (Vd < 0.0) {
            // Subsonic inflow (slide 15).
            pb = 0.5 * (qa[3] + qd[3] - rc * (nx * (qa[1] - qd[1]) + ny * (qa[2] - qd[2])));
            qb[0] = qa[0] + (pb - qa[3]) / (c0 * c0);
            qb[1] = qa[1] - nx * (qa[3] - pb) / rc;
            qb[2] = qa[2] - ny * (qa[3] - pb) / rc;
        } else {
            // Subsonic outflow (slide 16).
            pb = qa[3];
            qb[0] = qd[0] + (pb - qd[3]) / (c0 * c0);
            qb[1] = qd[1] + nx * (qd[3] - pb) / rc;
            qb[2] = qd[2] + ny * (qd[3] - pb) / rc;
        }
        qb[3] = pb;
    }

    void vortex_freestream(const Freestream& fs, double cl, double x, double y, double qa[NEQ])
    {
        const double M = fs.mach;                                  // free-stream Mach number (|V_inf|, c_inf = 1)
        const double a = fs.alpha_deg * PI / 180.0;                // angle of attack
        const double gam = 0.5 * M * VORTEX_CHORD * cl;            // circulation Gamma = 1/2 |V_inf| c CL
        const double beta = std::sqrt(1.0 - M * M);                // Prandtl-Glauert factor
        const double dx = x - VORTEX_XREF;
        const double dy = y - VORTEX_YREF;
        const double d = std::sqrt(dx * dx + dy * dy);             // distance to the vortex
        const double th = std::atan2(dy, dx);                      // polar angle theta
        const double s = std::sin(th - a);
        const double fac = gam * beta / (2.0 * PI * d * (1.0 - M * M * s * s)); // vortex strength at (x, y)
        const double u = fs.q[1] + fac * std::sin(th);
        const double v = fs.q[2] - fac * std::cos(th);
        const double e = (GAMMA - 1.0) / GAMMA;                    // (gamma - 1) / gamma
        const double pe = std::pow(fs.q[3], e) + e * fs.q[0] / std::pow(fs.q[3], 1.0 / GAMMA) * 0.5 * (M * M - u * u - v * v);
        qa[3] = std::pow(pe, 1.0 / e);                             // isentropic, same total enthalpy
        qa[0] = fs.q[0] * std::pow(qa[3] / fs.q[3], 1.0 / GAMMA);
        qa[1] = u;
        qa[2] = v;
    }

    void apply_bc(const Mesh& m, Flow& f, const Freestream& fs, JminBc jmin)
    {
        /*
        Farfield ghost from the boundary state: W_g = 2 Wb - Wd (conservative), both layers.
        (nx, ny) must point out of the domain; qa = the outer state (free stream, or with the vortex).
        */
        auto farfield_ghost = [&](int ci, int cg1, int cg2, double nx, double ny, const double qa[NEQ]) {
            double wd[NEQ]; // interior cell, conservative
            double qd[NEQ]; // interior cell, primitive
            double qb[NEQ]; // boundary state, primitive
            double wb[NEQ]; // boundary state, conservative
            for (int k = 0; k < NEQ; ++k) {
                wd[k] = f.W[k][ci];
            }
            cons_to_prim(wd, qd);
            farfield_state(qd, qa, nx, ny, qb);
            prim_to_cons(qb, wb);
            for (int k = 0; k < NEQ; ++k) {
                const double wg = 2.0 * wb[k] - wd[k]; // linear extrapolation to the ghost centre
                f.W[k][cg1] = wg;
                f.W[k][cg2] = wg;
            }
        };

        // jmax, farfield (after the wall, whose ghosts the vortex correction's CL needs).
        auto farfield_jmax = [&](double cl, bool vortex) {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < m.nic; ++i) {
                const int fc = m.node(i, m.nj - 1); // farfield face storage index; its normal points out (+j)
                double qa[NEQ];                     // outer state of this face
                if (vortex) {
                    const double xf = 0.5 * (m.x[m.node(i, m.nj - 1)] + m.x[m.node(i + 1, m.nj - 1)]); // face midpoint
                    const double yf = 0.5 * (m.y[m.node(i, m.nj - 1)] + m.y[m.node(i + 1, m.nj - 1)]);
                    vortex_freestream(fs, cl, xf, yf, qa);
                } else {
                    for (int k = 0; k < NEQ; ++k) {
                        qa[k] = fs.q[k];
                    }
                }
                farfield_ghost(m.cell(i, m.njc - 1), m.cell(i, m.njc), m.cell(i, m.njc + 1), m.njx[fc], m.njy[fc], qa);
            }
        };
        const bool vortex = FARFIELD_VORTEX && fs.mach < 1.0 && jmin == JminBc::Wall; // vortex correction on
        if (!vortex) {
            farfield_jmax(0.0, false);
        }

        // jmin: every i column is independent (its own wall ghosts): parallel loop.
#pragma omp parallel for schedule(static)
        for (int i = 0; i < m.nic; ++i) {
            // jmin: the j-face normal at j = 0 points into the fluid.
            const int fw = m.node(i, 0); // wall face storage index
            if (jmin == JminBc::Farfield) {
                farfield_ghost(m.cell(i, 0), m.cell(i, -1), m.cell(i, -2), -m.njx[fw], -m.njy[fw], fs.q); // flip: out of the domain
                continue;
            }

            const int ci = m.cell(i, 0); // first interior cell
            double wi[NEQ];              // interior cell, conservative
            double qi[NEQ];              // interior cell, primitive
            double qg[NEQ];              // ghost, primitive
            double wg[NEQ];              // ghost, conservative
            for (int k = 0; k < NEQ; ++k) {
                wi[k] = f.W[k][ci];
            }
            cons_to_prim(wi, qi);
            wall_ghost(qi, m.njx[fw], m.njy[fw], qg);
            if (WALL_PRESSURE_EXTRAPOLATE && WALL_EXTRAP_ACTIVE) {
                // Study option (bc.h): ghost pressure extrapolated from cells (i, 0) and (i, 1).
                const int c1 = m.cell(i, 1); // second cell off the wall
                double w1[NEQ];              // its conservative state
                double q1[NEQ];              // its primitive state
                for (int k = 0; k < NEQ; ++k) {
                    w1[k] = f.W[k][c1];
                }
                cons_to_prim(w1, q1);
                /*
                Only the pressure changes: rho_g = rho_i and the mirrored velocity are kept, so the
                face average 1/2 (W_g + W_i) still has exactly zero normal momentum (no mass or energy
                through the wall). (A first version also changed rho_g isentropically: that leaked mass.)
                */
                qg[3] = std::max(2.0 * qi[3] - q1[3], 0.2 * qi[3]); // extrapolated, kept positive
            }
            prim_to_cons(qg, wg);
            for (int k = 0; k < NEQ; ++k) {
                f.W[k][m.cell(i, -1)] = wg[k]; // ghost 1 (slide: cell 1)
                f.W[k][m.cell(i, -2)] = wg[k]; // ghost 2 (slide: cell 0), same rho and p as cell 2
            }
        }

        // Vortex correction: current CL from the wall pressure (wall ghosts are now set), then the farfield.
        if (vortex) {
            std::vector<double> pw; // wall-face pressure
            wall_pressure(m, f, pw);
            const ForceRef ref{VORTEX_CHORD, VORTEX_XREF, VORTEX_YREF}; // reference of CL
            farfield_jmax(integrate_forces(m, pw, fs, ref).cl, true);
        }

        /*
        i-cut, done last: cell (-g, j) is cell (nic-g, j), cell (nic-1+g, j) is cell (g-1, j).
        The copy also runs over the ghost rows j < 0 and j >= njc, which were just filled above:
        this fills the corner ghost blocks with the wall / farfield ghosts of the cells across
        the cut (needed by the cell -> node interpolation of the output at the TE nodes).
        */
        for (int k = 0; k < NEQ; ++k) {
            for (int j = -NG; j < m.njc + NG; ++j) {
                for (int g = 1; g <= NG; ++g) {
                    f.W[k][m.cell(-g, j)] = f.W[k][m.cell(m.nic - g, j)];
                    f.W[k][m.cell(m.nic - 1 + g, j)] = f.W[k][m.cell(g - 1, j)];
                }
            }
        }
    }

}
