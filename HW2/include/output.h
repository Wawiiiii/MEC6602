#pragma once

#include "mesh.h"

#include <string>

namespace euler {

    /*
       Write the mesh back with the cell areas (Euler2D slide 22, item 4),
       as a Tecplot ASCII file (format of Euler2D slide 20):
         X, Y at the nodes, Area at the cell centres, BLOCK packing.
       compute_metrics() must have been called first.
    */
    void write_tecplot_mesh(const Mesh& m, const std::string& filename);

}
