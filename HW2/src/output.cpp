#include "output.h"

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

}
