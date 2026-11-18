#pragma once

#include "flow.h"
#include "mesh.h"

namespace euler {

   // Make the freestream non dimentiojnal. 
    struct Freestream {
        double mach = 0.0;      // free-stream Mach number
        double alpha_deg = 0.0; // angle of attack in degrees
        double q[NEQ] = {};     // primitive state (rho, u, v, p)
    };

    Freestream make_freestream(double mach, double alpha_deg);

    // Fill every stored cell (ghosts included) with the free stream (cold start)
    void init_freestream(Flow& f, const Freestream& fs);

    void wall_ghost(const double qi[NEQ], double nx, double ny, double qg[NEQ]);

    /*
       Farfield boundary state Wb (nx, ny) = unit normal pointing OUT of the domain 
         supersonic inflow : Wb = Wa              supersonic outflow: Wb = Wd
         subsonic inflow   : pb = 1/2 {pa + pd - rho0 c0 [nx (ua - ud) + ny (va - vd)]}
                             rho_b = rho_a + (pb - pa)/c0^2
                             ub = ua - nx (pa - pb)/(rho0 c0),  vb = va - ny (pa - pb)/(rho0 c0)
         subsonic outflow  : pb = pa
                             rho_b = rho_d + (pb - pd)/c0^2
                             ub = ud + nx (pd - pb)/(rho0 c0),  vb = vd + ny (pd - pb)/(rho0 c0)
       Inflow / outflow: sign of the interior normal velocity V_d = ud nx + vd ny (< 0 = inflow).
    */
    void farfield_state(const double qd[NEQ], const double qa[NEQ], double nx, double ny, double qb[NEQ]);

    // yavai des bug ici, just pour sanity check sont pas utilises. 
    inline bool WALL_PRESSURE_EXTRAPOLATE = true;

    inline double WALL_EXTRAP_SWITCH = 1.0;

    inline bool WALL_EXTRAP_ACTIVE = true;

    // Fafield correction implemented but not used 
    inline bool FARFIELD_VORTEX = false;
    inline double VORTEX_CHORD = 1.0; // reference chord of CL in Gamma
    inline double VORTEX_XREF = 0.25; // vortex position x
    inline double VORTEX_YREF = 0.0;  // vortex position y

    // Free stream plus the point vortex of lift coefficient cl, at (x, y) (see FARFIELD_VORTEX).
    void vortex_freestream(const Freestream& fs, double cl, double x, double y, double qa[NEQ]);

    // What sits at jmin: the airfoil wall, or a farfield.
    enum class JminBc {
        Wall,
        Farfield
    };

    /*
       Fill the ghost cells 
         i-cut   : copy of the cells across the cut - i wrapping around the airfoil
         wall    : velocity reflected, density and pressure from the first cells
         farfield: W_g = 2 Wb - Wd, so the face average is exactly the boundary state Wb
       Only f.W is updated.
    */
    void apply_bc(const Mesh& m, Flow& f, const Freestream& fs, JminBc jmin = JminBc::Wall);

}
