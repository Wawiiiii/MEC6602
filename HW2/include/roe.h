#pragma once

#include "flow.h"

namespace euler {

    constexpr double HARTEN_FRACTION = 0.1; // entropy fix 

    void physical_flux(const double q[NEQ], double nx, double ny, double F[NEQ]);

    double harten(double lambda_abs, double delta);

    void roe_flux(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double F[NEQ],
                  double harten_fraction = HARTEN_FRACTION);

}
