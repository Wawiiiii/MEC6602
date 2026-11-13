#include "output.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace euler {

    void write_tecplot_mesh(const Mesh& m, const std::string& filename)
    {
        std::FILE* f = std::fopen(filename.c_str(), "w");
        if (!f) {
            throw std::runtime_error("cannot open output file: " + filename);
        }

        // Header: I, J are the NODE counts; variable 3 is given per cell (ni-1)*(nj-1).
        std::fprintf(f, "TITLE = \"mesh %dx%d\"\n", m.ni, m.nj);
        std::fprintf(f, "VARIABLES = \"X\", \"Y\", \"Area\"\n");
        std::fprintf(f, "ZONE T=\"%dx%d\", I=%d, J=%d, DATAPACKING=BLOCK, ZONETYPE=ORDERED, "
                        "VARLOCATION=([3]=CELLCENTERED)\n",
                     m.ni, m.nj, m.ni, m.nj);

        // BLOCK packing: all X (i fastest), then all Y, then all Area. Physical entries only.
        for (const std::vector<double>* v : {&m.x, &m.y}) {
            for (int j = 0; j < m.nj; ++j) {
                for (int i = 0; i < m.ni; ++i) {
                    std::fprintf(f, "%.12e\n", (*v)[m.node(i, j)]);
                }
            }
        }
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                std::fprintf(f, "%.12e\n", m.vol[m.cell(i, j)]);
            }
        }

        std::fclose(f);
    }

    // Write the values of one node variable, 8 per line (BLOCK packing).
    static void write_block(std::FILE* f, const std::vector<double>& v)
    {
        for (size_t k = 0; k < v.size(); ++k) {
            std::fprintf(f, "%.9e%c", v[k], (k % 8 == 7 || k + 1 == v.size()) ? '\n' : ' ');
        }
    }

    void write_tecplot_field(const Mesh& m, const Flow& f, const Freestream& fs, const std::string& filename)
    {
        const size_t nn = static_cast<size_t>(m.ni) * m.nj;         // number of nodes
        const double q_inf = 0.5 * fs.mach * fs.mach;               // dynamic pressure
        const double s_inf = fs.q[3] / std::pow(fs.q[0], GAMMA);    // free-stream p / rho^gamma
        const char* names[9] = {"X", "Y", "Rho", "U", "V", "P", "Mach", "Cp", "Entropy"};
        std::vector<std::vector<double>> var(9, std::vector<double>(nn)); // node values, i fastest

        for (int j = 0; j < m.nj; ++j) {
            for (int i = 0; i < m.ni; ++i) {
                const size_t n = static_cast<size_t>(j) * m.ni + i; // output node index
                // The 4 cells around node (i, j): (i-1, j-1), (i, j-1), (i-1, j), (i, j).
                const int cells[4] = {m.cell(i - 1, j - 1), m.cell(i, j - 1), m.cell(i - 1, j), m.cell(i, j)};
                double w[NEQ] = {0.0, 0.0, 0.0, 0.0}; // averaged conservative state
                for (int c : cells) {
                    for (int k = 0; k < NEQ; ++k) {
                        w[k] += 0.25 * f.W[k][c];
                    }
                }
                double q[NEQ]; // primitive state at the node
                cons_to_prim(w, q);
                const double a = std::sqrt(GAMMA * q[3] / q[0]); // speed of sound
                var[0][n] = m.x[m.node(i, j)];
                var[1][n] = m.y[m.node(i, j)];
                var[2][n] = q[0];
                var[3][n] = q[1];
                var[4][n] = q[2];
                var[5][n] = q[3];
                var[6][n] = std::sqrt(q[1] * q[1] + q[2] * q[2]) / a;
                var[7][n] = (q[3] - fs.q[3]) / q_inf;
                var[8][n] = q[3] / std::pow(q[0], GAMMA) / s_inf - 1.0;
            }
        }

        std::FILE* out = std::fopen(filename.c_str(), "w");
        if (!out) {
            throw std::runtime_error("cannot open output file: " + filename);
        }
        std::fprintf(out, "TITLE = \"Euler 2D, M = %.4f, alpha = %.4f deg\"\n", fs.mach, fs.alpha_deg);
        std::fprintf(out, "VARIABLES =");
        for (int k = 0; k < 9; ++k) {
            std::fprintf(out, " \"%s\"", names[k]);
        }
        std::fprintf(out, "\nZONE T=\"%dx%d\", I=%d, J=%d, DATAPACKING=BLOCK, ZONETYPE=ORDERED\n", m.ni, m.nj, m.ni, m.nj);
        for (int k = 0; k < 9; ++k) {
            write_block(out, var[k]);
        }
        std::fclose(out);
    }

    void write_surface_cp(const Mesh& m, const std::vector<double>& pw, const Freestream& fs, const ForceRef& ref,
                          const Forces& F, const std::string& filename)
    {
        std::FILE* out = std::fopen(filename.c_str(), "w");
        if (!out) {
            throw std::runtime_error("cannot open output file: " + filename);
        }
        const double q_inf = 0.5 * fs.mach * fs.mach; // dynamic pressure
        std::fprintf(out, "# Surface pressure, wall face midpoints, order TE -> lower -> LE -> upper -> TE\n");
        std::fprintf(out, "# M = %.6f  alpha = %.6f deg  chord = %.6f  xref = %.6f\n", fs.mach, fs.alpha_deg, ref.chord,
                     ref.xref);
        std::fprintf(out, "# CL = %.10f  CD = %.10f  CM = %.10f\n", F.cl, F.cd, F.cm);
        std::fprintf(out, "TITLE = \"surface Cp\"\nVARIABLES = \"x/c\" \"y/c\" \"Cp\" \"p\"\n");
        std::fprintf(out, "ZONE T=\"wall\", I=%d, DATAPACKING=POINT\n", m.nic);
        for (int i = 0; i < m.nic; ++i) {
            const int a = m.node(i, 0);
            const int b = m.node(i + 1, 0);
            std::fprintf(out, "%.10e %.10e %.10e %.10e\n", 0.5 * (m.x[a] + m.x[b]) / ref.chord,
                         0.5 * (m.y[a] + m.y[b]) / ref.chord, (pw[i] - fs.q[3]) / q_inf, pw[i]);
        }
        std::fclose(out);
    }

}
