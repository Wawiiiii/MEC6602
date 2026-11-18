#pragma once

namespace euler {

    struct FixedClSettings {
        bool on = false;          // fixed-CL mode (off = fixed alpha, the normal run)
        double cl_target = 0.4;   // target lift coefficient CL_t
        double dalpha0 = 0.1;     // first perturbation of alpha, degrees
        double drop = 1e-2;       // residual drop since the last alpha change before the next change
        double cl_tol = 1e-6;     // |CL - CL_t| accepted at the end
        double max_dalpha = 2.0;  // largest alpha change in one step, degrees
    };

    class FixedClController {
    public:
        FixedClController(const FixedClSettings& s, double alpha0_deg);

        /*
        Call once per nonlinear iteration with the density residual norm and CL of the current state.
        flow_converged: the residual has reached the solver's stopping tolerance.
        Returns true when alpha was changed
        */
        bool update(double norm, double cl, bool flow_converged);

        double alpha() const { return alpha_; }       // current angle of attack, degrees
        bool on_target(double cl) const;               // |CL - CL_t| <= cl_tol
        int changes() const { return changes_; }       // number of alpha changes so far

    private:
        FixedClSettings s_;       // settings
        double alpha_;            // current alpha, degrees
        double alpha_old_ = 0.0;  // alpha of the previous change, degrees
        double cl_old_ = 0.0;     // CL at the previous change
        bool have_old_ = false;   // false until the first perturbation
        double ref_ = -1.0;       // largest residual since the last change (-1: none yet)
        int changes_ = 0;         // alpha changes so far
    };

}
