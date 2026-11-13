#pragma once

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"
#include "residual.h"

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

namespace euler {

    /*
       Explicit multistage scheme (Euler2D_basic_algo slide 18):
         W(0) = W^n
         W(k) = W(0) - alpha_k dt/Omega R(W(k-1)),   k = 1 .. m
         W^(n+1) = W(m)
       The stage coefficients and the CFL come from the slide 18 table (upwind schemes).
       Explicit Euler (slide 17) is the 1-stage case alpha_1 = 1, CFL < 1.
    */
    struct RkScheme {
        int stages = 1;                            // number of stages m
        double alpha[5] = {1.0, 0.0, 0.0, 0.0, 0.0}; // stage coefficients alpha_1 .. alpha_m
        double cfl = 0.9;                          // CFL number sigma of the table
    };

    /*
       Coefficients of slide 18:
         spatial = Order::First  -> "1st-order scheme" column, stages 3, 4 or 5 (CFL 1.5, 2.0, 2.5)
         spatial = Order::Second -> "2nd-order scheme" column, stages 3, 4 or 5 (CFL 0.69, 0.92, 1.15)
       stages = 1 -> explicit Euler (slide 17), CFL 0.9.
       Throws for any other number of stages.
    */
    RkScheme rk_scheme(Order spatial, int stages);

    /*
       Local time step of every physical cell (slide 17, Roe_scheme_accel slide 20):
         dt = sigma Omega / (Lambda_I + Lambda_J),
         Lambda_I = (|V . n_I| + c) dS_I,  Lambda_J = (|V . n_J| + c) dS_J,
       with n and dS at the cell centre = average of the two faces bounding the cell in that
       direction (basic_algo slide 11). Uses f.Q (call after compute_residual / compute_primitives).
       dt has the size of the cell arrays; ghost slots are left untouched.
    */
    void local_time_step(const Mesh& m, const Flow& f, double cfl, std::vector<double>& dt);

    /*
       L2 norm of the density residual (post.pdf slide 4), on the physical cells:
         ||R_rho|| = sqrt( (1/N) sum (R_rho / Omega)^2 )
       R / Omega is d(rho)/dt, the rate of change of the density in the cell.
    */
    double density_residual_norm(const Mesh& m, const Flow& f);

    // Settings of one explicit run.
    struct ExplicitSettings {
        Order order = Order::Second;     // spatial order of the residual
        RkScheme rk;                     // time scheme (rk_scheme(...))
        JminBc jmin = JminBc::Wall;      // jmin boundary (Farfield only for the slide 20 check)
        int max_iter = 100000;           // iteration limit
        double tol_drop = 1e-12;         // stop when ||R_rho|| / ||R_rho||_0 < tol_drop
        int cauchy_window = 0;           // force (Cauchy) criterion window N in iterations (0 = off)
        double cauchy_tol = 1e-6;        // Cauchy criterion: max - min of CL, CD, CM over the window < tol
        ForceRef ref;                    // chord and moment reference point for CL, CD, CM
        int print_every = 0;             // print one line every N iterations (0 = never)
        std::string history_file;        // residual / force history file (empty = none)
    };

    // What a run did (explicit or implicit).
    struct RunResult {
        int iterations = 0;            // iterations done
        double res0 = 0.0;             // ||R_rho|| at the first iteration
        double res = 0.0;              // ||R_rho|| at the last iteration
        bool converged = false;        // residual drop reached
        bool cauchy = false;           // force criterion reached
        bool diverged = false;         // NaN or Inf in the residual
        double seconds = 0.0;          // wall-clock time of the iterations
        std::vector<double> history;   // ||R_rho|| of every iteration
        std::vector<double> cl;        // CL of every iteration
        std::vector<double> cd;        // CD of every iteration
        std::vector<double> cm;        // CM of every iteration
        Forces forces;                 // forces of the last iteration
    };
    using ExplicitResult = RunResult;

    /*
       Convergence monitoring (post.pdf slides 4-5), shared by the explicit and implicit solvers:
         - L2 norm of the density residual: stop when it has dropped by tol_drop;
         - Cauchy criterion on the forces: stop when, over the last cauchy_window iterations,
           max - min of CL, CD and CM are all below cauchy_tol (|F_m - F_n| < eps for all m, n);
         - NaN / Inf in the residual: stop, "diverged".
       Also writes the history file (iteration, ||R_rho||, drop, CL, CD, CM, seconds) and prints.
    */
    class Monitor {
    public:
        Monitor(int max_iter, double tol_drop, int cauchy_window, double cauchy_tol, int print_every,
                const std::string& history_file);
        ~Monitor();
        Monitor(const Monitor&) = delete;
        Monitor& operator=(const Monitor&) = delete;

        // Record iteration it; returns true when the run must stop (flags set in r).
        bool record(int it, double norm, const Forces& F, RunResult& r);

        // Seconds since the monitor was created.
        double seconds() const;

    private:
        int max_iter_;                             // iteration limit
        double tol_drop_;                          // residual drop to reach
        int window_;                               // Cauchy window (0 = off)
        double cauchy_tol_;                        // Cauchy tolerance
        int print_every_;                          // print period (0 = never)
        std::FILE* hist_ = nullptr;                // history file
        std::chrono::steady_clock::time_point t0_; // start time
    };

    /*
       Explicit multistage iterations with local time stepping (basic_algo slides 17-18), from the
       current f.W (ghosts are refilled by apply_bc at every stage, slide 9 step 4d).
       One iteration:
         1. apply_bc, R(W^n), dt from W^n (dt is frozen over the stages)
         2. stages k = 1 .. m: W(k) = W(0) - alpha_k dt/Omega R(W(k-1)), with apply_bc before each R
       The residual norm of the iteration is the one of W^n (stage 1).
    */
    RunResult run_explicit(const Mesh& m, Flow& f, const Freestream& fs, const ExplicitSettings& s);

}
