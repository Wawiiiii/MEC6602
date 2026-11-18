#pragma once

#include "bc.h"
#include "flow.h"
#include "mesh.h"

#include <vector>

namespace euler {

    // Reference values for the force and moment coefficients 
    struct ForceRef {
        double chord = 1.0; // reference chord c - mettre a 1 mais jameson la met un peu plus que 1 
        double xref = 0.25; // moment reference point x 
        double yref = 0.0;  // moment reference point y
    };

    // Forces on the airfoil and their coefficients.
    struct Forces {
        double fx = 0.0; // px - pression y 
        double fy = 0.0; // py - pression y
        double mz = 0.0; // moment anti-horaire --  +
        double cl = 0.0; // lift coefficient - perpendiculaire au free-stream (sign error?)
        double cd = 0.0; // drag coefficient - avec le freestream
        double cm = 0.0; // nose up + nose down - 
    };


    // 1/2 (W_ghost + W_interior) 
    void wall_pressure(const Mesh& m, const Flow& f, std::vector<double>& pw);

    // extrpolation 2D pour conserver l'ordre 
    void wall_pressure_extrapolated(const Mesh& m, const Flow& f, std::vector<double>& pw);

    // fomrules pour trouver les coeffs a partir de fy et fx 
    Forces integrate_forces(const Mesh& m, const std::vector<double>& pw, const Freestream& fs, const ForceRef& ref);

    // wall_pressure + integrate_forces.
    Forces compute_forces(const Mesh& m, const Flow& f, const Freestream& fs, const ForceRef& ref);

}
