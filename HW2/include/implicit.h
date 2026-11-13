#pragma once

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"
#include "residual.h"
#include "solver.h"

#include <string>

namespace euler {

    /*
       Implicit time schemes for the steady state (Roe_scheme_accel slides 25-40). With dt -> infinity
       all of them become Newton's method on R(W) = 0 (slide 38). Pseudo-time, local time steps.
         BDF1 : implicit Euler (slide 25):  (Omega/dt + J) dW = -R(W^n)
         BDF2 : second-order backward difference (added at the user's request, not in the slides):
                (3/2 Omega/dt + J) dW = -R(W^n) + Omega/(2 dt) dW^(n-1)     (first iteration = BDF1)
         ROW  : Rosenbrock-Wanner ROS2, 2 stages, L-stable, gamma = 1 + 1/sqrt(2) (added at the
                user's request, not in the slides), same matrix M = Omega/(gamma dt) + J for both:
                M k1 = -R(W^n) / gamma
                M k2 = -R(W^n + k1) / gamma - 2 Omega/(gamma dt) k1
                W^(n+1) = W^n + 3/2 k1 + 1/2 k2
       J = dR/dW of the full (second-order) residual.
    */
    enum class TimeScheme {
        BDF1,
        BDF2,
        ROW
    };

    /*
       How each linear system is solved:
         GMRES : Jacobian-free Newton-Krylov (slides 39-40): GMRES on the exact operator
                 c Omega/dt + J, J v ~ (R(W + eps v) - R(W)) / eps, right-preconditioned by
                 block SGS sweeps on the first-order Roe Jacobian (slides 26-35).
         SGS   : block SGS sweeps only, on the first-order operator (slide 33: approximate Jacobian,
                 a few forward/backward sweeps, the update uses the 2nd-order residual on the right).
    */
    enum class LinearSolver {
        GMRES,
        SGS
    };

    // Settings of one implicit run.
    struct ImplicitSettings {
        Order order = Order::Second;            // spatial order of the residual
        TimeScheme scheme = TimeScheme::BDF1;   // time scheme
        LinearSolver linsolve = LinearSolver::GMRES; // linear solver
        JminBc jmin = JminBc::Wall;             // jmin boundary
        double cfl0 = 1.0;                      // starting CFL (10 fails on the supersonic 257x257 cold start)
        double cfl_max = 1e8;                   // largest CFL (CFL ramp of slide 40)
        int gmres_restart = 30;                 // GMRES restart length (slide 39: 20-50)
        int gmres_max_iter = 60;                // GMRES iterations per linear solve
        double gmres_rtol = 1e-3;               // GMRES relative tolerance (inexact Newton)
        int sgs_sweeps = 2;                     // SGS sweeps (preconditioner, or the solver itself)
        double sgs_omega = 1.0;                 // SGS relaxation factor
        double max_change = 0.5;                // update limiter: |d rho|/rho and |d p|/p at most this
        int max_iter = 500;                     // nonlinear iteration limit
        double tol_drop = 1e-12;                // stop when ||R_rho|| has dropped by this factor
        int cauchy_window = 0;                  // Cauchy criterion window (0 = off)
        double cauchy_tol = 1e-6;               // Cauchy tolerance
        ForceRef ref;                           // force reference
        int print_every = 1;                    // print every N iterations
        std::string history_file;               // history file (empty = none)
    };

    /*
       Implicit iterations from the current f.W (cold start = free stream).
       Each iteration: apply_bc, R(W^n), forces, convergence check; CFL ramp
       sigma^n = sigma^0 ||R^0|| / ||R^n|| clamped to [cfl0, cfl_max] (slide 40); local dt;
       first-order Jacobian + c Omega/dt on the diagonal; linear solve(s); update with the
       positivity limiter (the update is scaled down if it would change rho or p by more than
       max_change in some cell).
    */
    RunResult run_implicit(const Mesh& m, Flow& f, const Freestream& fs, const ImplicitSettings& s);

}
