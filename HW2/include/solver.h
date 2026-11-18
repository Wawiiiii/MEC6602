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

  // partie 1 - just pour verifier que tout est bon
    struct RkScheme {
        int stages = 1;                            
        double alpha[5] = {1.0, 0.0, 0.0, 0.0, 0.0}; 
        double cfl = 0.9;                          
    };


    RkScheme rk_scheme(Order spatial, int stages);
 
    
    void local_time_step(const Mesh& m, const Flow& f, double cfl, std::vector<double>& dt);

  
    double density_residual_norm(const Mesh& m, const Flow& f);

    struct ExplicitSettings {
        Order order = Order::Second;     // spatial order of the residual
        RkScheme rk;                     // time scheme (rk_scheme(...))
        JminBc jmin = JminBc::Wall;      // jmin boundary (Farfield only for the uniform-flow test)
        int max_iter = 100000;           // iteration limit
        double tol_drop = 1e-12;         // stop when ||R_rho|| / ||R_rho||_0 < tol_drop
        int cauchy_window = 0;           // force (Cauchy) criterion window N in iterations (0 = off)
        double cauchy_tol = 1e-6;        // Cauchy criterion: max - min of CL, CD, CM over the window < tol
        ForceRef ref;                    // chord and moment reference point for CL, CD, CM
        int print_every = 0;             // print one line every N iterations (0 = never)
        std::string history_file;        // residual / force history file (empty = none)
    };

    // status de la sim -> mostly sanity check on tests
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
        double alpha_deg = 0.0;        // angle of attack at the end 
        int alpha_changes = 0;         // alpha changes
    };
    using ExplicitResult = RunResult;


    // Claude a fait un monitor avec un fallback pour des simulations lancees avec des mauvaises entrees. 
    // par exemple si pas assez de GMRES restart et stall 
    class Monitor {
    public:
        Monitor(int max_iter, double tol_drop, int cauchy_window, double cauchy_tol, int print_every,
                const std::string& history_file);
        ~Monitor();
        Monitor(const Monitor&) = delete;
        Monitor& operator=(const Monitor&) = delete;

        // Record iteration it; returns true when the run must stop 
        bool record(int it, double norm, const Forces& F, RunResult& r);

        // Seconds since the monitor was created.
        double seconds() const;

        void set_extra_columns(const std::string& header);
        bool row_open() const { return row_open_; } 
        void end_row();                              

    private:
        std::string extra_header_;     
        bool header_printed_ = false;  // column header already printed
        bool row_open_ = false;        // last row printed without its end of line
        int max_iter_;                             // iteration limit
        double tol_drop_;                          // residual drop to reach
        int window_;                               // Cauchy window (0 = off)
        double cauchy_tol_;                        // Cauchy tolerance
        int print_every_;                          // print period (0 = never)
        std::FILE* hist_ = nullptr;                // history file
        std::chrono::steady_clock::time_point t0_; // start time
    };

 
    RunResult run_explicit(const Mesh& m, Flow& f, const Freestream& fs, const ExplicitSettings& s);

}
