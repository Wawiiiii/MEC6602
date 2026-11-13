#pragma once

#include "flow.h"
#include "mesh.h"

namespace euler {

    /*
       Free stream, non-dimensional (the code is non-dimensional, post.pdf slides 13-15):
         rho_inf = 1, p_inf = 1/gamma (so c_inf = 1), |V_inf| = M,
         u_inf = M cos(alpha), v_inf = M sin(alpha).
    */
    struct Freestream {
        double mach = 0.0;      // free-stream Mach number
        double alpha_deg = 0.0; // angle of attack in degrees
        double q[NEQ] = {};     // primitive state (rho, u, v, p)
    };

    Freestream make_freestream(double mach, double alpha_deg);

    // Fill every stored cell (ghosts included) with the free stream (cold start, post.pdf slide 8).
    void init_freestream(Flow& f, const Freestream& fs);

    /*
       Slip wall ghost state (Euler2D_basic_algo slide 14): V.n = 0 at the wall.
         rho_g = rho_i,  p_g = p_i,  v_g = v_i - 2 (v_i . n) n
       qi: first interior cell, (nx, ny): wall unit normal (either direction).
    */
    void wall_ghost(const double qi[NEQ], double nx, double ny, double qg[NEQ]);

    /*
       Farfield boundary state Wb (slides 15-16). (nx, ny) = unit normal pointing OUT of the
       domain (the textbook convention the slide formulas use), qd = interior cell (d),
       qa = free stream (a), rho0 and c0 taken from the interior cell.
         supersonic inflow : Wb = Wa              supersonic outflow: Wb = Wd
         subsonic inflow   : pb = 1/2 {pa + pd - rho0 c0 [nx (ua - ud) + ny (va - vd)]}
                             rho_b = rho_a + (pb - pa)/c0^2
                             ub = ua - nx (pa - pb)/(rho0 c0),  vb = va - ny (pa - pb)/(rho0 c0)
         subsonic outflow  : pb = pa
                             rho_b = rho_d + (pb - pd)/c0^2
                             ub = ud + nx (pd - pb)/(rho0 c0),  vb = vd + ny (pd - pb)/(rho0 c0)
       Inflow / outflow: sign of the interior normal velocity V_d = ud nx + vd ny (< 0 = inflow).
       Supersonic / subsonic: |V_d| >= c0 or not, i.e. the normal Mach number, which decides how
       many characteristics (V - c, V, V + c along n) enter the domain.
    */
    void farfield_state(const double qd[NEQ], const double qa[NEQ], double nx, double ny, double qb[NEQ]);

    /*
       Wall ghost pressure (case-file key wall_pressure):
         false ("mirror", slide 14): p_g = p_i (P1 = P2), i.e. dp/dn = 0 at the wall;
         true ("extrapolated", DEFAULT): p_g = 2 p_0 - p_1 from the first two cells off the wall
               (linear extrapolation, so the wall-face average is the pressure extrapolated to the
               wall); rho_g = rho_0 and the mirrored velocity are kept, so no mass crosses the wall.
       On a curved wall dp/dn = rho V^2 / R != 0, so the slide's dp/dn = 0 is only first-order accurate:
       with it CM converges at order ~1 and CL at ~1.5; with the extrapolation CL, CD, CM all converge
       at order >= 2 (DEVIATIONS.md A20, grid study 129 / 257 / 513 at M = 0.5).
    */
    inline bool WALL_PRESSURE_EXTRAPOLATE = true;

    /*
       Start-up of the extrapolated wall pressure (only used when WALL_PRESSURE_EXTRAPOLATE is true):
       the solvers use the mirror ghost (slide 14) until the density residual has dropped by
       WALL_EXTRAP_SWITCH (1 = extrapolated from the first iteration, default), then switch to the
       extrapolated ghost for good. Optional robustness knob; the converged solution only depends on
       the final BC. (The supersonic cold-start problem found on 257x257 was the starting CFL, fixed by
       cfl0 = 1, not the wall BC.)
       WALL_EXTRAP_ACTIVE is the switch state, set by run_explicit / run_implicit; apply_bc reads it.
    */
    inline double WALL_EXTRAP_SWITCH = 1.0;
    inline bool WALL_EXTRAP_ACTIVE = true;

    // What sits at jmin: the airfoil wall, or (slide 20 sanity check only) a farfield.
    enum class JminBc {
        Wall,
        Farfield
    };

    /*
       Fill all the ghost cells used by the residual from the physical cells (slides 13-16):
         i-cut (imin/imax): ghost cells = the cells across the cut (connectivity).
         wall (jmin)      : both ghost layers = wall_ghost() of the first interior cell
                            (slide 14: rho, p of cell 2 in both ghosts, velocity reflected).
         farfield (jmax)  : Wb from farfield_state(), then extrapolated to the ghost centre
                            in conservative variables, W_g = 2 Wb - Wd (slide 16, last line),
                            so that the boundary-face average 1/2 (W_g + W_d) is exactly Wb.
                            The second farfield ghost layer = the first one.
       The i-cut copy is done last and also covers the ghost rows, so the corner ghost blocks
       (ghost in i AND in j) hold the wall / farfield ghosts of the cells across the cut.
       Updates f.W only; compute_residual() derives f.Q from it.
    */
    void apply_bc(const Mesh& m, Flow& f, const Freestream& fs, JminBc jmin = JminBc::Wall);

}
