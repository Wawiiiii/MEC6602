#pragma once

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"

#include <string>
#include <vector>

namespace euler {

    void write_tecplot_mesh(const Mesh& m, const std::string& filename);


    void write_tecplot_field(const Mesh& m, const Flow& f, const Freestream& fs, const std::string& filename);


    void write_surface_cp(const Mesh& m, const std::vector<double>& pw, const Freestream& fs, const ForceRef& ref,
                          const Forces& F, const std::string& filename);

  
    struct RunSummary {
        int ni = 0;            // grid nodes in i
        int nj = 0;            // grid nodes in j
        double mach = 0.0;     // free-stream Mach number
        double alpha_deg = 0.0; // angle of attack, degrees
        int order = 2;         // spatial order
        std::string scheme;    // solver label, e.g. "bdf1-gmres", "rk5"
        std::string status;    // converged, cauchy, max_iter or diverged
        int iterations = 0;    // iterations (or multigrid cycles) of the run
        double seconds = 0.0;  // wall-clock time of the solver iterations
        double drop = 0.0;     // final density residual / first one
        Forces forces;         // final CL, CD, CM
    };


    void append_run_summary(const std::string& filename, const RunSummary& s, const ForceRef& ref);

}
