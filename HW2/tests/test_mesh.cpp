/*
   Step 1 check: read a plot3D mesh, verify the indexing, print the nodes.
   Usage: test_mesh <mesh.x>
*/

#include "mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mesh.x>\n", argv[0]);
        return 1;
    }

    try {
        const euler::Mesh m = euler::read_plot3d(argv[1]);
        const int NG = euler::NG;

        std::printf("nodes : %d x %d   (with ghosts %d x %d)\n", m.ni, m.nj, m.ni_tot, m.nj_tot);
        std::printf("cells : %d x %d   (with ghosts %d x %d)\n", m.nic, m.njc, m.nic_tot, m.njc_tot);

        // Index sanity: the corner ghost nodes/cells must map to the two ends of the array.
        std::printf("node(-NG, -NG)           = %d   (expect 0)\n", m.node(-NG, -NG));
        std::printf("node(ni-1+NG, nj-1+NG)   = %d   (expect %zu)\n",
                    m.node(m.ni - 1 + NG, m.nj - 1 + NG), m.x.size() - 1);
        std::printf("cell(-NG, -NG)           = %d   (expect 0)\n", m.cell(-NG, -NG));
        std::printf("cell(nic-1+NG, njc-1+NG) = %d   (expect %d)\n",
                    m.cell(m.nic - 1 + NG, m.njc - 1 + NG), m.nic_tot * m.njc_tot - 1);

        // Farfield (j = nj-1) radius range.
        double rmin = 1e300, rmax = 0.0;
        for (int i = 0; i < m.ni; ++i) {
            const int k = m.node(i, m.nj - 1);
            const double r = std::hypot(m.x[k], m.y[k]);
            rmin = std::min(rmin, r);
            rmax = std::max(rmax, r);
        }
        std::printf("farfield radius: %.3f .. %.3f\n", rmin, rmax);
        std::printf("O-grid cut closed (checked by the reader)\n");

        // Full node table for small meshes only: compare with the file / a sketch.
        if (m.ni <= 17) {
            std::printf("\n   j    i            x            y\n");
            for (int j = 0; j < m.nj; ++j) {
                for (int i = 0; i < m.ni; ++i) {
                    const int k = m.node(i, j);
                    std::printf("%4d %4d %12.6f %12.6f\n", j, i, m.x[k], m.y[k]);
                }
            }
        }
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
