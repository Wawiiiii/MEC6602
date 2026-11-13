#include "forces.h"

#include <cmath>

namespace euler {

    constexpr double PI = 3.14159265358979323846; // MSVC has no M_PI by default

    void wall_pressure(const Mesh& m, const Flow& f, std::vector<double>& pw)
    {
        pw.resize(m.nic);
        for (int i = 0; i < m.nic; ++i) {
            double w01[NEQ]; // averaged conservative state of the wall face
            double q01[NEQ]; // same, primitive
            for (int k = 0; k < NEQ; ++k) {
                w01[k] = 0.5 * (f.W[k][m.cell(i, -1)] + f.W[k][m.cell(i, 0)]);
            }
            cons_to_prim(w01, q01);
            pw[i] = q01[3];
        }
    }

    void wall_pressure_extrapolated(const Mesh& m, const Flow& f, std::vector<double>& pw)
    {
        pw.resize(m.nic);
        for (int i = 0; i < m.nic; ++i) {
            const int a = m.node(i, 0);     // wall face, first node
            const int b = m.node(i + 1, 0); // wall face, second node
            const double xm = 0.5 * (m.x[a] + m.x[b]); // face midpoint
            const double ym = 0.5 * (m.y[a] + m.y[b]);
            const double nx = m.njx[a];     // wall normal into the fluid
            const double ny = m.njy[a];
            double p[2]; // pressure of cells (i, 0) and (i, 1)
            double d[2]; // their normal distances to the wall face
            for (int jj = 0; jj < 2; ++jj) {
                const int c = m.cell(i, jj);
                double w[NEQ];
                double q[NEQ];
                for (int k = 0; k < NEQ; ++k) {
                    w[k] = f.W[k][c];
                }
                cons_to_prim(w, q);
                p[jj] = q[3];
                const double xc = 0.25 * (m.x[m.node(i, jj)] + m.x[m.node(i + 1, jj)] + m.x[m.node(i + 1, jj + 1)] +
                                          m.x[m.node(i, jj + 1)]);
                const double yc = 0.25 * (m.y[m.node(i, jj)] + m.y[m.node(i + 1, jj)] + m.y[m.node(i + 1, jj + 1)] +
                                          m.y[m.node(i, jj + 1)]);
                d[jj] = (xc - xm) * nx + (yc - ym) * ny;
            }
            pw[i] = p[0] - d[0] * (p[1] - p[0]) / (d[1] - d[0]);
        }
    }

    Forces integrate_forces(const Mesh& m, const std::vector<double>& pw, const Freestream& fs, const ForceRef& ref)
    {
        Forces F;
        for (int i = 0; i < m.nic; ++i) {
            const int fc = m.node(i, 0);                          // wall face storage index
            const int n2 = m.node(i + 1, 0);                      // second node of the face
            const double xm = 0.5 * (m.x[fc] + m.x[n2]);          // face midpoint x
            const double ym = 0.5 * (m.y[fc] + m.y[n2]);          // face midpoint y
            const double dp = pw[i] - fs.q[3];                    // p_w - p_inf
            const double fx = -dp * m.njx[fc] * m.sj[fc];         // force of this face on the airfoil, x
            const double fy = -dp * m.njy[fc] * m.sj[fc];         // same, y
            F.fx += fx;
            F.fy += fy;
            F.mz += (xm - ref.xref) * fy - (ym - ref.yref) * fx;
        }
        const double a = fs.alpha_deg * PI / 180.0;  // angle of attack
        const double q = 0.5 * fs.mach * fs.mach;    // dynamic pressure (rho_inf = 1, |V_inf| = M)
        F.cl = (F.fy * std::cos(a) - F.fx * std::sin(a)) / (q * ref.chord);
        F.cd = (F.fx * std::cos(a) + F.fy * std::sin(a)) / (q * ref.chord);
        F.cm = -F.mz / (q * ref.chord * ref.chord);
        return F;
    }

    Forces compute_forces(const Mesh& m, const Flow& f, const Freestream& fs, const ForceRef& ref)
    {
        std::vector<double> pw; // wall pressure
        wall_pressure(m, f, pw);
        return integrate_forces(m, pw, fs, ref);
    }

}
