#include "mesh.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace euler {

    // Read a single-block 2D plot3D file into the physical nodes (ghost nodes stay NaN).
    Mesh read_plot3d(const std::string& filename)
    {
        std::ifstream in(filename);
        if (!in) {
            throw std::runtime_error("cannot open mesh file: " + filename);
        }

        // Header: number of blocks (always 1 here), then the node counts.
        int nblock = 0;
        if (!(in >> nblock) || nblock != 1) {
            throw std::runtime_error("expected a single-block plot3D file: " + filename);
        }

        Mesh m;
        if (!(in >> m.ni >> m.nj) || m.ni < 3 || m.nj < 3) {
            throw std::runtime_error("bad grid dimensions in " + filename);
        }

        // Derived sizes. Example 9x9 file: 9 nodes -> 8 cells; with ghosts 13 nodes, 12 cells.
        m.nic     = m.ni - 1;       // physical cells in i (around the airfoil)
        m.njc     = m.nj - 1;       // physical cells in j (wall to farfield)
        m.ni_tot  = m.ni + 2 * NG;  // stored nodes in i, ghosts included
        m.nj_tot  = m.nj + 2 * NG;  // stored nodes in j, ghosts included
        m.nic_tot = m.nic + 2 * NG; // stored cells in i, ghosts included
        m.njc_tot = m.njc + 2 * NG; // stored cells in j, ghosts included

    
        const double nan = std::numeric_limits<double>::quiet_NaN();
        m.x.assign(static_cast<size_t>(m.ni_tot) * m.nj_tot, nan);
        m.y.assign(m.x.size(), nan);

        // plot3D order: all x (i fastest), then all y.
        for (std::vector<double>* coord : {&m.x, &m.y}) {
            for (int j = 0; j < m.nj; ++j) {
                for (int i = 0; i < m.ni; ++i) {
                    if (!(in >> (*coord)[m.node(i, j)])) {
                        throw std::runtime_error("mesh file ended early: " + filename);
                    }
                }
            }
        }

        double extra = 0.0;
        if (in >> extra) {
            throw std::runtime_error("unexpected extra values in " + filename);
        }

        // O-grid cut at the trailing edge: the i = 0 and i = ni-1 lines must be the same points.
        for (int j = 0; j < m.nj; ++j) {
            const int a = m.node(0, j);
            const int b = m.node(m.ni - 1, j);
            if (std::abs(m.x[a] - m.x[b]) > 1e-10 || std::abs(m.y[a] - m.y[b]) > 1e-10) {
                throw std::runtime_error("i-min and i-max lines do not coincide in " + filename);
            }
        }

        return m;
    }

    // Cell areas, face lengths and unit face normals.
    void compute_metrics(Mesh& m)
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const std::vector<double>& x = m.x;
        const std::vector<double>& y = m.y;

        m.vol.assign(static_cast<size_t>(m.nic_tot) * m.njc_tot, nan);
        m.nix.assign(x.size(), nan);
        m.niy.assign(x.size(), nan);
        m.si.assign(x.size(), nan);
        m.njx.assign(x.size(), nan);
        m.njy.assign(x.size(), nan);
        m.sj.assign(x.size(), nan);

        /*
           Cell areas: Omega = 0.5 |a x b|, a and b the two diagonals.
           The signed area must be positive (corners counter-clockwise): the normals below rely on it.
        */
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int n00 = m.node(i, j);         // corner (i, j)
                const int n10 = m.node(i + 1, j);     // corner (i+1, j)
                const int n11 = m.node(i + 1, j + 1); // corner (i+1, j+1)
                const int n01 = m.node(i, j + 1);     // corner (i, j+1)

                const double ax = x[n11] - x[n00]; // diagonal a: (i,j) -> (i+1,j+1)
                const double ay = y[n11] - y[n00];
                const double bx = x[n01] - x[n10]; // diagonal b: (i+1,j) -> (i,j+1)
                const double by = y[n01] - y[n10];

                const double area = 0.5 * (ax * by - ay * bx); // signed area
                if (!(area > 0.0)) {
                    throw std::runtime_error("non-positive cell area at cell (" + std::to_string(i) +
                                             ", " + std::to_string(j) + ")");
                }
                m.vol[m.cell(i, j)] = area;
            }
        }

        /*
           Ghost-cell areas, for the MUSCL limiter eps (muscl.h) at boundary faces.
             i-cut: the ghost cells are the cells across the cut -> copy them.
             wall and farfield: mirror of the interior cell (ghost -1 <-> cell 0, ghost -2 <-> cell 1).
           The four corner ghost blocks are never used and stay NaN.
        */
        for (int j = 0; j < m.njc; ++j) {
            for (int g = 1; g <= NG; ++g) {
                m.vol[m.cell(-g, j)] = m.vol[m.cell(m.nic - g, j)];        // left of the cut = last cells
                m.vol[m.cell(m.nic - 1 + g, j)] = m.vol[m.cell(g - 1, j)]; // right of the cut = first cells
            }
        }
        for (int i = 0; i < m.nic; ++i) {
            for (int g = 1; g <= NG; ++g) {
                m.vol[m.cell(i, -g)] = m.vol[m.cell(i, g - 1)];                 // wall mirror
                m.vol[m.cell(i, m.njc - 1 + g)] = m.vol[m.cell(i, m.njc - g)]; // farfield mirror
            }
        }

        /*
           i-faces: edge t = node(i,j+1) - node(i,j). Normal (ty, -tx) / |t|,
           pointing from cell (i-1, j) to cell (i, j).
        */
        for (int j = 0; j < m.njc; ++j) {
            for (int i = 0; i < m.ni; ++i) {
                const int a = m.node(i, j);     // first node of the face (storage index)
                const int b = m.node(i, j + 1); // second node of the face
                const double tx = x[b] - x[a];
                const double ty = y[b] - y[a];
                const double len = std::hypot(tx, ty); // face length Delta S
                m.si[a] = len;
                m.nix[a] = ty / len;
                m.niy[a] = -tx / len;
            }
        }

        /*
           j-faces: edge t = node(i+1,j) - node(i,j). Normal (-ty, tx) / |t|,
           pointing from cell (i, j-1) to cell (i, j). At the wall (j = 0) it points into the fluid.
        */
        for (int j = 0; j < m.nj; ++j) {
            for (int i = 0; i < m.nic; ++i) {
                const int a = m.node(i, j);     // first node of the face (storage index)
                const int b = m.node(i + 1, j); // second node of the face
                const double tx = x[b] - x[a];
                const double ty = y[b] - y[a];
                const double len = std::hypot(tx, ty); // face length Delta S
                m.sj[a] = len;
                m.njx[a] = -ty / len;
                m.njy[a] = tx / len;
            }
        }
    }

}
