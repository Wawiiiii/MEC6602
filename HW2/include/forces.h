#pragma once

#include "bc.h"
#include "flow.h"
#include "mesh.h"

#include <vector>

namespace euler {

    // Reference values for the force and moment coefficients (post.pdf slides 16-18).
    struct ForceRef {
        double chord = 1.0; // reference chord c (open item: 1.0 or the grid's 1.00893, DEVIATIONS.md)
        double xref = 0.25; // moment reference point x (quarter chord)
        double yref = 0.0;  // moment reference point y
    };

    // Forces on the airfoil and their coefficients.
    struct Forces {
        double fx = 0.0; // pressure force on the airfoil, x (grid axes)
        double fy = 0.0; // pressure force on the airfoil, y
        double mz = 0.0; // moment about (xref, yref), counter-clockwise positive
        double cl = 0.0; // lift coefficient, perpendicular to the free stream
        double cd = 0.0; // drag coefficient, along the free stream
        double cm = 0.0; // pitching moment coefficient about (xref, yref), nose-up positive
    };

    /*
       Pressure on every wall face i = 0..nic-1: the pressure of the averaged state
       1/2 (W_ghost + W_interior) of the face, i.e. exactly the pressure the wall flux uses
       (basic_algo slide 7, residual.cpp), so the forces are the momentum fluxes through the wall.
       apply_bc must have been called on f (ghosts up to date).
    */
    void wall_pressure(const Mesh& m, const Flow& f, std::vector<double>& pw);

    /*
       Wall pressure extrapolated linearly to the wall face from the first two cells off the wall
       (i, 0) and (i, 1), using their normal distances to the face (d = (x_cell - x_face) . n):
         p_w = p0 - d0 (p1 - p0) / (d1 - d0)
       Second-order accurate at the wall, whereas wall_pressure() is the first-cell value
       (half a cell away, first order). Used for the force integration test of the grid
       convergence study (DEVIATIONS.md B6). Cell centres = average of the 4 corners.
    */
    void wall_pressure_extrapolated(const Mesh& m, const Flow& f, std::vector<double>& pw);

    /*
       Integrate a wall pressure (one value per wall face) into forces and coefficients
       (post.pdf slides 16-18, integrals replaced by sums over the wall faces):
         F  = - sum (p_w - p_inf) n dS          (n = wall normal into the fluid)
         Mz = sum [ (x_m - xref) F_y,face - (y_m - yref) F_x,face ],  (x_m, y_m) = face midpoint
         CL = (F_y cos(alpha) - F_x sin(alpha)) / (q_inf c),  CD = (F_x cos(alpha) + F_y sin(alpha)) / (q_inf c)
         CM = -Mz / (q_inf c^2)   (nose-up positive),  q_inf = rho_inf V_inf^2 / 2 = M^2 / 2
       Subtracting p_inf changes nothing on a closed contour but reduces round-off.
    */
    Forces integrate_forces(const Mesh& m, const std::vector<double>& pw, const Freestream& fs, const ForceRef& ref);

    // wall_pressure + integrate_forces.
    Forces compute_forces(const Mesh& m, const Flow& f, const Freestream& fs, const ForceRef& ref);

}
