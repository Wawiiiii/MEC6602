#include "output.h"

#include <cmath>
#include <cstdio>
#include <ctime>
#include <stdexcept>

namespace euler {

    // Tecplot file of the mesh nodes and the cell areas.
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

    // Tecplot file of the flow at the nodes: each node gets the average of its 4 cells.
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

    // Cp at the wall face midpoints, with the force coefficients in the header.
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

    // Append one line per run to the summary file.
    void append_run_summary(const std::string& filename, const RunSummary& s, const ForceRef& ref)
    {
        std::FILE* out = std::fopen(filename.c_str(), "a"); // append: earlier runs are kept
        if (!out) {
            throw std::runtime_error("cannot open summary file: " + filename);
        }

        // Header only when the file is new (empty).
        std::fseek(out, 0, SEEK_END);
        if (std::ftell(out) == 0) {
            std::fprintf(out, "# HW2 Euler solver - run summary, one line per run (new runs are appended)\n");
            std::fprintf(out, "# CL, CD, CM : chord %.5f, moment about (%.4f, %.4f), CM nose-up positive\n", ref.chord,
                         ref.xref, ref.yref);
            std::fprintf(out, "# time_s     : wall-clock seconds of the solver iterations\n");
            std::fprintf(out, "# drop       : final / first density residual\n");
            std::fprintf(out, "# status     : converged (residual drop reached), cauchy (force criterion),\n");
            std::fprintf(out, "#              max_iter (iteration limit, NOT converged), diverged\n");
            std::fprintf(out, "#\n");
            std::fprintf(out, "# %-11s %6s %7s %17s %17s %17s %6s %9s  %-10s %9s  %-11s %5s  %s\n", "grid", "Mach",
                         "alpha", "CL", "CD", "CM", "iter", "time_s", "status", "drop", "scheme", "order", "date");
        }

        // Date and time of the end of the run.
        char date[32];
        const std::time_t now = std::time(nullptr);
        std::strftime(date, sizeof(date), "%Y-%m-%d %H:%M", std::localtime(&now));

        char grid[32];
        std::snprintf(grid, sizeof(grid), "%dx%d", s.ni, s.nj);

        // Forces in scientific notation so that small values (e.g. CD -> 0) keep their digits.
        std::fprintf(out, "  %-11s %6.3f %7.3f %17.10e %17.10e %17.10e %6d %9.2f  %-10s %9.2e  %-11s %5d  %s\n", grid,
                     s.mach, s.alpha_deg, s.forces.cl, s.forces.cd, s.forces.cm, s.iterations, s.seconds,
                     s.status.c_str(), s.drop, s.scheme.c_str(), s.order, date);
        std::fclose(out);
    }

}
