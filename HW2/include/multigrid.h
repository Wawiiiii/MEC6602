#pragma once

#include "bc.h"
#include "flow.h"
#include "forces.h"
#include "mesh.h"
#include "residual.h"
#include "solver.h"

#include <array>
#include <string>
#include <vector>

// Ne pas utiliser, experimentation. j'ai pas reussi a bien la link ca bug
namespace euler {
    // Smoother of the multigrid cycle.
    enum class MgSmoother {
        RK,
        SGS
    };

    struct MultigridSettings {
        MgSmoother smoother = MgSmoother::RK; // smoother on every level
        double sgs_cfl = 10.0;               // CFL of the SGS smoother (implicit Euler step)
        int sgs_sweeps = 2;                  // SGS sweeps per smoothing step
        int levels = 0;                      // number of grid levels 
        int gamma = 2;                       // 1 = V-cycle, 2 = W-cycle
        int n_pre = 1;                       // smoothing iterations before going down
        int n_post = 1;                      // smoothing iterations after coming back up
        int n_coarse = 2;                    // iterations on the coarsest level
        int fmg_cycles = 0;                  // full multigrid start
        Order fine_order = Order::Second;    // spatial order on the finest level
        Order coarse_order = Order::First;   // spatial order on the coarse levels
        RkScheme fine_rk;                   
        RkScheme coarse_rk;                 
        double max_change = 0.5;            
        double relax = 0.0;                  // relaxation of the correction into a 2nd-order level
        JminBc jmin = JminBc::Wall;          // jmin boundary
        int max_iter = 2000;                 // number of cycles
        double tol_drop = 1e-12;             // stop when the fine ||R_rho|| has dropped by this
        int cauchy_window = 0;               // Cauchy criterion window, in cycles 
        double cauchy_tol = 1e-6;            // Cauchy tolerance
        ForceRef ref;                        
        int print_every = 10;                
        std::string history_file;           
    };

    // Coarse mesh made of every other node of m 
    Mesh coarsen_mesh(const Mesh& m);

    // Can m be coarsened once more? - c'est devision en 4 ici
    bool can_coarsen(const Mesh& m);

 
    void restrict_volume_average(const Mesh& mf, const Mesh& mc, const std::vector<double>& fine,
                                 std::vector<double>& coarse);

    void restrict_sum(const Mesh& mf, const Mesh& mc, const std::vector<double>& fine, std::vector<double>& coarse);
    void prolong_bilinear(const Mesh& mf, const Mesh& mc, const std::vector<double>& coarse, std::vector<double>& fine);

 
    RunResult run_multigrid(const Mesh& m, Flow& f, const Freestream& fs, const MultigridSettings& s);

    // Relaxation factor of the correction prolonged into a level of the given order (see relax above).
    double correction_relax(const MultigridSettings& s, Order level_order);

    // Number of levels run_multigrid would use for this mesh and settings.
    int multigrid_levels(const Mesh& m, const MultigridSettings& s);

  
    std::vector<Mesh> coarse_hierarchy(const Mesh& m, int n_levels);

    // Fine cells = value of their coarse cell, for every variable.
    void inject_piecewise_constant(const Mesh& mf, const Mesh& mc, const Flow& coarse, Flow& fine);

}
