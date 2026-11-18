#pragma once

#include "bc.h"
#include "fixedcl.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"
#include "residual.h"
#include "solver.h"

#include <string>

namespace euler {

    // on va just prendre bdf1 - je voulais experimenter un peu :)
    enum class TimeScheme {
        BDF1,
        BDF2,
        ROW
    };

    /*
       How each linear system is solved:
         GMRES : Jacobian-free Newton-Krylov: GMRES on the exact operator
                 c Omega/dt + J, J v ~ (R(W + eps v) - R(W)) / eps, right-preconditioned by
                 block SGS sweeps on the first-order Roe Jacobian.
         SGS   : block SGS sweeps only
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
        double cfl_max = 1e8;                   // largest CFL of the CFL ramp
        int gmres_restart = 30;                 // GMRES restart length (typically 20-50)
        int gmres_max_iter = 60;                // GMRES iterations per linear solve
        double gmres_rtol = 1e-3;               // GMRES relative tolerance (inexact Newton)
        double gmres_fail = 0.5;                // linear solve failed if its relative residual is above this (CFL backtracking)
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
        FixedClSettings fixed_cl;               // fixed-CL mode, secant on alpha (fixedcl.h)
    };

    RunResult run_implicit(const Mesh& m, Flow& f, const Freestream& fs, const ImplicitSettings& s);

}
