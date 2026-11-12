#pragma once

#include "flow.h"

namespace euler {

    constexpr double HARTEN_FRACTION = 0.1; // entropy fix width delta = c~ / 10 (Roe_scheme_accel slide 11)

    /*
       Physical convective flux through a face of unit normal (nx, ny), per unit face length
       (Euler2D slide 7, 2D):
         Fc = ( rho V,  rho u V + nx p,  rho v V + ny p,  rho H V ),   V = u nx + v ny
       with rho H = rho E + p. q = (rho, u, v, p).
    */
    void physical_flux(const double q[NEQ], double nx, double ny, double F[NEQ]);

    /*
       Harten entropy correction of one eigenvalue (slide 11):
         |lambda|                          if |lambda| >  delta
         (lambda^2 + delta^2) / (2 delta)  if |lambda| <= delta
       lambda_abs = |lambda|. delta = 0 switches the correction off.
    */
    double harten(double lambda_abs, double delta);

    /*
       Roe flux at one face, per unit face length (slides 7-12):
         F = 1/2 [ Fc(WR) + Fc(WL) - |A_Roe| (WR - WL) ]
       qL, qR: primitive states (rho, u, v, p) on each side (cell values at first order,
       MUSCL values at second order). (nx, ny): unit normal from L to R.
       harten_fraction: delta = harten_fraction * c~ (slide 11); 0 turns the fix off
       (only used by the tests).
    */
    void roe_flux(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double F[NEQ],
                  double harten_fraction = HARTEN_FRACTION);

}
