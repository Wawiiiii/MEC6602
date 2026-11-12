#pragma once

#include "mesh.h"

#include <array>
#include <vector>

namespace euler {

    constexpr int NEQ = 4;        // number of equations in 2D: mass, x-momentum, y-momentum, energy
    constexpr double GAMMA = 1.4; // ratio of specific heats (air)

    /*
       Flow variables, one value per cell (cell-centred finite volume, Euler2D slides 8-9).
       Stored like the mesh: one flat vector per variable, index Mesh::cell(i, j),
       ghost cells included. Everything starts as NaN, so a ghost cell that the
       boundary conditions forgot to fill shows up as NaN in the residual.
    */
    struct Flow {
        std::array<std::vector<double>, NEQ> W; // conservative variables: rho, rho*u, rho*v, rho*E
        std::array<std::vector<double>, NEQ> Q; // primitive variables: rho, u, v, p (computed from W)
        std::array<std::vector<double>, NEQ> R; // residual: sum over the cell faces of Fc*dS, normal pointing out
    };

    // Allocate W, Q and R for every stored cell of the mesh (ghosts included), all NaN.
    void allocate_flow(const Mesh& m, Flow& f);

    // Fill Q from W in every stored cell (ghosts included; a NaN cell simply gives NaN).
    void compute_primitives(Flow& f);

    /*
       Primitive -> conservative (Euler2D slide 7):
         q = (rho, u, v, p)  ->  w = (rho, rho u, rho v, rho E),
         rho E = p / (gamma - 1) + rho (u^2 + v^2) / 2
    */
    inline void prim_to_cons(const double q[NEQ], double w[NEQ])
    {
        w[0] = q[0];
        w[1] = q[0] * q[1];
        w[2] = q[0] * q[2];
        w[3] = q[3] / (GAMMA - 1.0) + 0.5 * q[0] * (q[1] * q[1] + q[2] * q[2]);
    }

    // Conservative -> primitive: p = (gamma - 1) (rho E - rho (u^2 + v^2) / 2).
    inline void cons_to_prim(const double w[NEQ], double q[NEQ])
    {
        const double rho = w[0];     // density
        const double u = w[1] / rho; // x velocity
        const double v = w[2] / rho; // y velocity
        q[0] = rho;
        q[1] = u;
        q[2] = v;
        q[3] = (GAMMA - 1.0) * (w[3] - 0.5 * rho * (u * u + v * v));
    }

}
