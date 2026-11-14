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

namespace euler {

    /*
       Geometric FAS multigrid (Roe_scheme_accel slides 41-53), explicit multistage smoother.

       Coarse grids (slide 41): every other node of the finer grid in i and j, i.e. a coarse cell
       (I, J) is made of the four fine cells (2I, 2J), (2I+1, 2J), (2I, 2J+1), (2I+1, 2J+1).
       The grids are (2^n + 1)^2 nodes, so they coarsen down to 9x9 (8x8 cells).

       One cycle on level l (slide 44: 1 W -> 2 R -> 3 time steps on 2h -> 4 correction):
         1. pre-smoothing: n_pre multistage iterations on R_l(W) + Q_l
         2. restriction (slide 49):
              W_{l+1}^(0) = sum Omega_k W_k / sum Omega_k          (volume-weighted average)
              Q_{l+1} = sum_k (R_l + Q_l)_k - R_{l+1}(W_{l+1}^(0))  (forcing function, slides 46, 51)
         3. gamma cycles on level l+1 (gamma = 1: V-cycle, gamma = 2: W-cycle, slides 52-53);
            on the coarsest level: n_coarse iterations instead
         4. prolongation of the correction (slides 48, 50):
              W_l += relax I (W_{l+1} - W_{l+1}^(0)),  bilinear weights 9/16, 3/16, 3/16, 1/16 of the
              nearest coarse cells (slide 50); relax < 1 only into a 2nd-order level (see below)
         5. post-smoothing: n_post iterations
       Smoother on every level ("any temporal scheme can be used", slide 47):
         RK  : the multistage scheme of slide 18 with local time steps,
               W(k) = W(0) - alpha_k dt/Omega (R + Q)   (classic explicit multigrid)
         SGS : one implicit Euler step (slide 25) with block SGS sweeps on the first-order Jacobian
               of the level (slides 26-35), (Omega/dt + J1) dW = -(R + Q)   (implicit multigrid,
               Cagnone et al. 2011, in the slide references); CFL of the level = sgs_cfl.
       The fine level uses the second-
       order residual; the coarse levels use the first-order residual (cheaper and more robust; the
       converged fine solution does not depend on the coarse operator, slide 46: if R_h = 0 the
       forcing cancels and there is no correction).
       Full multigrid start (slide 42): optionally converge on the coarse grids first and inject
       the solution into the next finer grid (fine cells = their coarse cell) as initial condition.

       Two stability settings (DEVIATIONS.md B18, B19, measured by power iteration on the cycle
       linearised about the converged solution, NACA0012 65x65):
         relax   : the correction prolonged into a SECOND-order level is multiplied by relax
                   (W_l += relax I (W_{l+1} - W_{l+1}^(0))). Automatic value (relax <= 0): 0.5 with
                   the RK smoother, 1 with SGS; corrections into first-order levels are never relaxed.
                   Why: at the trailing edge the 2nd-order fine operator is twice as stiff as the
                   1st-order coarse one for the mode "pressure jump across the TE" (measured
                   <I A_h e, A_H I e> / |A_H I e|^2 = 1.98 there), so the coarse grid corrects that
                   mode about twice too much; the RK smoother does not damp it (one RK5 step even
                   amplifies it transiently, x1.29, x1.41 after two), so the cycle had an eigenvalue
                   of -1.7 (V, 3 levels) to -1.95 (V, 4 levels): period-2 limit cycle = the stall.
                   relax = 0.5 = 1 / 1.98 removes the overshoot (eigenvalues 0.90-0.92).
                   The implicit SGS smoother damps that mode itself and needs no relaxation.
         sgs_cfl : 10 (was 20). The SGS step on the 2nd-order residual with the 1st-order Jacobian
                   is a defect correction; at M 0.8 the limiter makes an odd-even mode along the
                   wall at the leading edge about 2x stiffer in R than in J1, and the step overshoots
                   it for CFL > ~15 (eigenvalue of one step -1.05 at CFL 20, -1.44 at CFL 40,
                   +0.96 at 15 and 10): period-2 orbit of the smoother itself, even on one grid.
    */
    // Smoother of the multigrid cycle.
    enum class MgSmoother {
        RK,
        SGS
    };

    struct MultigridSettings {
        MgSmoother smoother = MgSmoother::RK; // smoother on every level
        double sgs_cfl = 10.0;               // CFL of the SGS smoother (implicit Euler step), see below
        int sgs_sweeps = 2;                  // SGS sweeps per smoothing step
        int levels = 0;                      // number of grid levels (0 = as many as possible)
        int gamma = 2;                       // 1 = V-cycle, 2 = W-cycle
        int n_pre = 1;                       // smoothing iterations before going down
        int n_post = 1;                      // smoothing iterations after coming back up
        int n_coarse = 2;                    // iterations on the coarsest level
        int fmg_cycles = 0;                  // full multigrid start: cycles on each coarse level (0 = off)
        Order fine_order = Order::Second;    // spatial order on the finest level
        Order coarse_order = Order::First;   // spatial order on the coarse levels
        RkScheme fine_rk;                    // multistage scheme on the finest level
        RkScheme coarse_rk;                  // multistage scheme on the coarse levels
        double max_change = 0.5;             // limiter of the prolonged correction (rho, p)
        double relax = 0.0;                  // relaxation of the correction into a 2nd-order level (<= 0: auto)
        JminBc jmin = JminBc::Wall;          // jmin boundary
        int max_iter = 2000;                 // number of cycles
        double tol_drop = 1e-12;             // stop when the fine ||R_rho|| has dropped by this
        int cauchy_window = 0;               // Cauchy criterion window, in cycles (0 = off)
        double cauchy_tol = 1e-6;            // Cauchy tolerance
        ForceRef ref;                        // force reference
        int print_every = 10;                // print every N cycles
        std::string history_file;            // history file (one line per cycle)
    };

    // Coarse mesh made of every other node of m (in i and j). Throws if m cannot be coarsened.
    Mesh coarsen_mesh(const Mesh& m);

    // Can m be coarsened once more (even cell counts, at least 8 x 4 cells left)?
    bool can_coarsen(const Mesh& m);

    /*
       Restriction and prolongation helpers (also used by the tests). fine / coarse are the cell
       arrays of the two levels (ghost slots included, only physical cells are used / written).
         restrict_volume_average : coarse = sum Omega_k fine_k / sum Omega_k  (slide 49)
         restrict_sum            : coarse = sum fine_k                         (slide 49, residuals)
         prolong_bilinear        : fine = bilinear interpolation of coarse    (slide 50)
       The prolongation uses the coarse values of the neighbours, periodic across the TE cut and
       constant (copy of the first row) across the wall and the farfield.
    */
    void restrict_volume_average(const Mesh& mf, const Mesh& mc, const std::vector<double>& fine,
                                 std::vector<double>& coarse);
    void restrict_sum(const Mesh& mf, const Mesh& mc, const std::vector<double>& fine, std::vector<double>& coarse);
    void prolong_bilinear(const Mesh& mf, const Mesh& mc, const std::vector<double>& coarse, std::vector<double>& fine);

    /*
       Multigrid iterations from the current f.W on mesh m (finest level). Returns the history of the
       finest level, one entry per cycle (RunResult.iterations = number of cycles).
    */
    RunResult run_multigrid(const Mesh& m, Flow& f, const Freestream& fs, const MultigridSettings& s);

    // Relaxation factor of the correction prolonged into a level of the given order (see relax above).
    double correction_relax(const MultigridSettings& s, Order level_order);

    // Number of levels run_multigrid would use for this mesh and settings.
    int multigrid_levels(const Mesh& m, const MultigridSettings& s);

    /*
       Full multigrid start for ANY solver (slide 42): returns the list of coarse meshes of m,
       coarsest first (n_levels - 1 of them), so that the caller can solve on each one, inject the
       solution into the next finer one (inject_piecewise_constant) and finally start the fine run.
    */
    std::vector<Mesh> coarse_hierarchy(const Mesh& m, int n_levels);

    // Fine cells = value of their coarse cell (slide 42, "the easiest way"), for every variable.
    void inject_piecewise_constant(const Mesh& mf, const Mesh& mc, const Flow& coarse, Flow& fine);

}
