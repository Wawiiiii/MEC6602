#pragma once

#include <cmath>

namespace euler {

    /*
       MUSCL reconstruction of the face states 

       Face I+1/2 between cells I and I+1:
         left state  U_L = U_I     + 1/2 delta_L,   from cells I-1, I, I+1
         right state U_R = U_{I+1} - 1/2 delta_R,   from cells I, I+1, I+2
       delta is the limited slope. It is written with two differences:
         a = the jump ACROSS the face        (U_{I+1} - U_I for both states)
         b = the jump on the FAR side        (U_I - U_{I-1} for U_L, U_{I+2} - U_{I+1} for U_R)
       The right state is the mirror image of the left one. 
       way round; for Hemker-Koren, which is not symmetric in a and b, that order would not give the
       kappa^ = 1/3 scheme, so the mirrored form is used.)

       Each primitive variable (rho, u, v, p) is reconstructed on its own.
    */

    // Settings read from the case file 
    inline double LIMITER_K = 5.0;       // K in eps = (K h)^3 -> plus ca c'est rand plus moins le limiteur influence so un off switch 
    inline bool VAN_ALBADA = false;      // false: Hemker-Koren (kappa^ = 1/3), true: Van Albada (kappa^ = 0)
    inline bool MUSCL_KAPPA_ON = false;  // true: unlimited kappa scheme instead (study option, smooth flow only)
    inline double MUSCL_KAPPA = 1.0 / 3.0; // kappa^ of the unlimited scheme

    // kh^3 ... maybe a mistake?
    inline double limiter_eps(double area)
    {
        const double k_times_h = LIMITER_K * std::sqrt(area); // K h
        return k_times_h * k_times_h * k_times_h;
    }

    // Unlimited kappa scheme: delta = 1/2 [(1 - kappa^) b + (1 + kappa^) a].
    inline double kappa_slope(double a, double b)
    {
        return 0.5 * ((1.0 - MUSCL_KAPPA) * b + (1.0 + MUSCL_KAPPA) * a);
    }

    // Hemker-Koren limited slope, kappa^ = 1/3.
    inline double hk_slope(double a, double b, double eps)
    {
        return ((2.0 * a * a + eps) * b + (b * b + 2.0 * eps) * a) / (2.0 * a * a + 2.0 * b * b - a * b + 3.0 * eps);
    }

    // Van Albada limited slope, kappa^ = 0; symmetric in a and b.
    inline double va_slope(double a, double b, double eps)
    {
        return (a * (b * b + eps) + b * (a * a + eps)) / (a * a + b * b + 2.0 * eps);
    }

    // Slope of the selected scheme.
    inline double muscl_slope(double a, double b, double eps)
    {
        if (MUSCL_KAPPA_ON) {
            return kappa_slope(a, b);
        }
        if (VAN_ALBADA) {
            return va_slope(a, b, eps);
        }
        return hk_slope(a, b, eps);
    }

    // Left state of face I+1/2. eps: limiter eps of cell I.
    inline double muscl_left(double u_far_left, double u_left, double u_right, double eps)
    {
        const double across = u_right - u_left;     // a: jump across the face
        const double far_side = u_left - u_far_left; // b: jump on the far side
        return u_left + 0.5 * muscl_slope(across, far_side, eps);
    }

    // Right state of face I+1/2. eps: limiter eps of cell I+1.
    inline double muscl_right(double u_left, double u_right, double u_far_right, double eps)
    {
        const double across = u_right - u_left;       // a: jump across the face
        const double far_side = u_far_right - u_right; // b: jump on the far side
        return u_right - 0.5 * muscl_slope(across, far_side, eps);
    }

}
