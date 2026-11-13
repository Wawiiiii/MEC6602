#pragma once

#include <cmath>

namespace euler {

    /*
       Constant K of the limiter eps = (K h)^3. Default 5 (DEVIATIONS.md A14). It can be changed from
       the case file (key limiter_k) for studies: a very large K switches the limiter off (pure
       kappa = 1/3 MUSCL), only acceptable for smooth (subsonic) flow.
    */
    inline double LIMITER_K = 5.0;

    /*
       Limiter eps of a cell of area Omega: eps = (K h)^3, h = sqrt(Omega) the cell size
       (Venkatakrishnan's scaling, J. Comput. Phys. 118, 1995).
       The slide 16 suggestion eps = Omega^2 = h^4 was implemented first: with it the 2nd-order
       residual stalls at a drop of 1e-3 on the 65x65 and 129x129 grids (limit cycle of the limiter in
       the wall row near the leading edge), exactly the limiter-induced stall Venkatakrishnan's eps
       was designed to remove. With (5 h)^3 the residual goes to machine zero and the spurious drag
       at M = 0.5 is lower (129x129: CD 0.0024 vs 0.0031 stalled). Comparison in DEVIATIONS.md A14.
    */
    inline double limiter_eps(double area)
    {
        const double kh = LIMITER_K * std::sqrt(area); // K h
        return kh * kh * kh;
    }

    /*
       MUSCL reconstruction with the Hemker-Koren limiter, kappa^ = 1/3
       (Roe_scheme_accel slides 13-16). Applied to each primitive variable separately.

       Face I+1/2 between cells I and I+1 (slide 15):
         UL = U_I     + 1/2 delta_L
         UR = U_{I+1} - 1/2 delta_R

       Smooth limited slope of slide 16:
         delta(a, b) = [ (2 a^2 + eps) b + (b^2 + 2 eps) a ] / (2 a^2 + 2 b^2 - a b + 3 eps)

       a is the difference ACROSS the face, b the difference on the FAR side:
         left state : a = U_{I+1} - U_I,  b = U_I - U_{I-1}      (a_L = D+U_I, b_L = D-U_I as on slide 16)
         right state: a = U_{I+1} - U_I,  b = U_{I+2} - U_{I+1}

       Note on the right state: slide 16 lists a_R = D+U_{I+1}, b_R = D-U_{I+1}. That labelling only
       matters for Hemker-Koren (van Albada is symmetric in a, b). Writing the slide 15 form
       UR = U_{I+1} - 1/2 Psi(r_R) (U_{I+2} - U_{I+1}), r_R = (U_{I+1} - U_I)/(U_{I+2} - U_{I+1}),
       with Psi(r) = 1/2 [(1 + kappa^) r + (1 - kappa^)] psi(r), psi(r) = 3r/(2r^2 - r + 2),
       gives exactly delta(a, b) with a = U_{I+1} - U_I and b = U_{I+2} - U_{I+1}, i.e. the
       mirror image of the left state. test_roe checks this against slide 15 and by symmetry.

       eps: small number (slide 16), given by limiter_eps(). It must be > 0 when a = b = 0.

       eps keeps the limiter from clipping smooth extrema (there a, b = O(h^2)), while at a
       discontinuity a^2 = O(1) >> eps and the limiter acts fully. The formula is NOT switched off
       when a b < 0 (no TVD "psi = 0 for r <= 0"): that switch was tried (2026-11-12) and made the
       residual error on a smooth field 6x larger (DEVIATIONS.md A13).
       Cost: at a sharp extremum the face value can overshoot (data 0.5, 1, 0 -> 1.125);
       on monotone data the face value always stays between its two neighbour cells (test_roe B5).
    */
    inline double hk_slope(double a, double b, double eps)
    {
        return ((2.0 * a * a + eps) * b + (b * b + 2.0 * eps) * a) / (2.0 * a * a + 2.0 * b * b - a * b + 3.0 * eps);
    }

    // Left state at face I+1/2 from cells I-1, I, I+1. eps = (area of cell I)^2.
    inline double muscl_left(double Um1, double U0, double Up1, double eps)
    {
        return U0 + 0.5 * hk_slope(Up1 - U0, U0 - Um1, eps);
    }

    // Right state at face I+1/2 from cells I, I+1, I+2. eps = (area of cell I+1)^2.
    inline double muscl_right(double U0, double Up1, double Up2, double eps)
    {
        return Up1 - 0.5 * hk_slope(Up1 - U0, Up2 - Up1, eps);
    }

}
