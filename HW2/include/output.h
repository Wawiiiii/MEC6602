#pragma once

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"

#include <string>
#include <vector>

namespace euler {

    /*
       Write the mesh back with the cell areas (Euler2D slide 22, item 4),
       as a Tecplot ASCII file (format of Euler2D slide 20):
         X, Y at the nodes, Area at the cell centres, BLOCK packing.
       compute_metrics() must have been called first.
    */
    void write_tecplot_mesh(const Mesh& m, const std::string& filename);

    /*
       Field solution, Tecplot ASCII, ordered zone, BLOCK packing, all variables at the nodes
       (Euler2D slide 20: "interpolate the variables from cell-centre to the node").
       Node value = average of the conservative variables of the 4 cells around the node, ghost
       cells included (wall: mirror ghost -> tangential velocity at the wall; farfield: average =
       the boundary state; i-cut: the cells across the cut), then converted to primitive.
       Variables: X, Y, Rho, U, V, P, Mach, Cp, Entropy = p/rho^gamma / (p/rho^gamma)_inf - 1
       (0 in exact subsonic Euler flow: a direct measure of the numerical error).
       apply_bc must have been called on f (all ghosts, corners included, up to date).
    */
    void write_tecplot_field(const Mesh& m, const Flow& f, const Freestream& fs, const std::string& filename);

    /*
       Surface pressure, one line per wall face (face midpoint), in grid order
       TE -> lower surface -> LE -> upper surface -> TE. Columns: x/c, y/c, Cp, p.
       Tecplot ASCII (POINT packing) with '#' comment lines, also readable by numpy.loadtxt.
       pw = wall_pressure(), Cp = (p_w - p_inf) / q_inf.
    */
    void write_surface_cp(const Mesh& m, const std::vector<double>& pw, const Freestream& fs, const ForceRef& ref,
                          const Forces& F, const std::string& filename);

}
