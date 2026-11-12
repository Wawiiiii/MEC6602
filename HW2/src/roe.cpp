#include "roe.h"

#include <cmath>

namespace euler {

    void physical_flux(const double q[NEQ], double nx, double ny, double F[NEQ])
    {
        const double rho = q[0];                                                  // density
        const double u = q[1];                                                    // x velocity
        const double v = q[2];                                                    // y velocity
        const double p = q[3];                                                    // pressure
        const double V = u * nx + v * ny;                                         // contravariant velocity V = V.n
        const double rhoH = GAMMA / (GAMMA - 1.0) * p + 0.5 * rho * (u * u + v * v); // rho H = rho E + p

        F[0] = rho * V;
        F[1] = rho * u * V + nx * p;
        F[2] = rho * v * V + ny * p;
        F[3] = rhoH * V;
    }

    double harten(double lambda_abs, double delta)
    {
        // ">=" so that delta = 0 (fix off) never divides by zero; both branches give delta at |lambda| = delta.
        if (lambda_abs >= delta) {
            return lambda_abs;
        }
        return (lambda_abs * lambda_abs + delta * delta) / (2.0 * delta);
    }

    void roe_flux(const double qL[NEQ], const double qR[NEQ], double nx, double ny, double F[NEQ],
                  double harten_fraction)
    {
        // Left and right states (step 1 of slide 12 is done by the caller).
        const double rhoL = qL[0]; // left density
        const double uL = qL[1];   // left x velocity
        const double vL = qL[2];   // left y velocity
        const double pL = qL[3];   // left pressure
        const double rhoR = qR[0]; // right density
        const double uR = qR[1];   // right x velocity
        const double vR = qR[2];   // right y velocity
        const double pR = qR[3];   // right pressure

        // Total enthalpy H = gamma/(gamma-1) p/rho + (u^2 + v^2)/2 on each side.
        const double HL = GAMMA / (GAMMA - 1.0) * pL / rhoL + 0.5 * (uL * uL + vL * vL); // left total enthalpy
        const double HR = GAMMA / (GAMMA - 1.0) * pR / rhoR + 0.5 * (uR * uR + vR * vR); // right total enthalpy

        /*
        Step 2: Roe average state (slide 9). Weights sqrt(rho_L) and sqrt(rho_R):
          rho~ = sqrt(rhoL rhoR),  u~ = (uL sqrt(rhoL) + uR sqrt(rhoR)) / (sqrt(rhoL) + sqrt(rhoR)),
          same for v~ and H~,  c~ = sqrt((gamma-1)(H~ - q~^2/2)),  V~ = u~ nx + v~ ny
        */
        const double sL = std::sqrt(rhoL);                       // sqrt of left density
        const double sR = std::sqrt(rhoR);                       // sqrt of right density
        const double wL = sL / (sL + sR);                        // Roe weight of the left state
        const double wR = sR / (sL + sR);                        // Roe weight of the right state
        const double rho_t = sL * sR;                            // rho~
        const double u_t = wL * uL + wR * uR;                    // u~
        const double v_t = wL * vL + wR * vR;                    // v~
        const double H_t = wL * HL + wR * HR;                    // H~
        const double q2_t = u_t * u_t + v_t * v_t;               // q~^2
        const double c_t = std::sqrt((GAMMA - 1.0) * (H_t - 0.5 * q2_t)); // c~
        const double V_t = u_t * nx + v_t * ny;                  // V~ (contravariant)

        // Jumps Delta() = ()_R - ()_L (slide 9).
        const double drho = rhoR - rhoL;    // density jump
        const double du = uR - uL;          // x velocity jump
        const double dv = vR - vL;          // y velocity jump
        const double dp = pR - pL;          // pressure jump
        const double dV = du * nx + dv * ny; // contravariant velocity jump

        // Step 3: eigenvalues V~ - c~, V~, V~ + c~ with the Harten entropy correction (slide 11).
        const double delta = harten_fraction * c_t;            // entropy fix width
        const double l1 = harten(std::abs(V_t - c_t), delta);   // |lambda_1|,   acoustic wave V~ - c~
        const double l234 = harten(std::abs(V_t), delta);       // |lambda_2,3,4|, entropy and shear waves V~
        const double l5 = harten(std::abs(V_t + c_t), delta);   // |lambda_5|,   acoustic wave V~ + c~

        /*
        Step 4: |A_Roe| (WR - WL) = |dF1| + |dF2,3,4| + |dF5| (slide 8).
          |dF1|     = |V~-c~| (dp - rho~ c~ dV)/(2 c~^2) * (1, u~ - c~ nx, v~ - c~ ny, H~ - c~ V~)
          |dF5|     = |V~+c~| (dp + rho~ c~ dV)/(2 c~^2) * (1, u~ + c~ nx, v~ + c~ ny, H~ + c~ V~)
          |dF2,3,4| = |V~| { (drho - dp/c~^2) * (1, u~, v~, q~^2/2)
                            + rho~ * (0, du - dV nx, dv - dV ny, u~ du + v~ dv - V~ dV) }
        */
        const double c2 = c_t * c_t;                                // c~^2
        const double s1 = l1 * (dp - rho_t * c_t * dV) / (2.0 * c2); // |lambda_1| x strength of wave 1
        const double s5 = l5 * (dp + rho_t * c_t * dV) / (2.0 * c2); // |lambda_5| x strength of wave 5
        const double s2 = l234 * (drho - dp / c2);                   // |lambda_2| x entropy wave strength
        const double s3 = l234 * rho_t;                              // |lambda_3,4| x rho~ (shear waves)

        double D[NEQ]; // |A_Roe| (WR - WL), the upwind dissipation
        D[0] = s1 + s2 + s5;
        D[1] = s1 * (u_t - c_t * nx) + s2 * u_t + s3 * (du - dV * nx) + s5 * (u_t + c_t * nx);
        D[2] = s1 * (v_t - c_t * ny) + s2 * v_t + s3 * (dv - dV * ny) + s5 * (v_t + c_t * ny);
        D[3] = s1 * (H_t - c_t * V_t) + s2 * 0.5 * q2_t + s3 * (u_t * du + v_t * dv - V_t * dV)
             + s5 * (H_t + c_t * V_t);

        // Step 5: F = 1/2 [Fc(WR) + Fc(WL) - |A_Roe| dW] (slide 8). The caller multiplies by dS.
        double FL[NEQ]; // physical flux of the left state
        double FR[NEQ]; // physical flux of the right state
        physical_flux(qL, nx, ny, FL);
        physical_flux(qR, nx, ny, FR);
        for (int k = 0; k < NEQ; ++k) {
            F[k] = 0.5 * (FL[k] + FR[k] - D[k]);
        }
    }

}
