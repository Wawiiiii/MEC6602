#include "fixedcl.h"

#include <algorithm>
#include <cmath>

namespace euler {

    FixedClController::FixedClController(const FixedClSettings& s, double alpha0_deg) : s_(s), alpha_(alpha0_deg)
    {
    }

    bool FixedClController::on_target(double cl) const
    {
        return std::abs(cl - s_.cl_target) <= s_.cl_tol;
    }

    // Fixed-CL step: change alpha once the flow is converged enough and CL is off target.
    bool FixedClController::update(double norm, double cl, bool flow_converged)
    {
        if (!s_.on) {
            return false;
        }
        // Largest residual since the last change: reference of the "converged enough" test.
        ref_ = std::max(ref_, norm);
        const bool ready = flow_converged || norm <= s_.drop * ref_; // converged enough
        if (!ready || on_target(cl)) {
            return false; // keep iterating the flow at this alpha
        }

        const double err = s_.cl_target - cl; // CL_t - CL
        double step = 0.0;                    // alpha change, degrees
        if (!have_old_ || cl == cl_old_) {
            step = std::copysign(s_.dalpha0, err); // first step: small perturbation towards the target
        } else {
            // Secant: dalpha = (CL_t - CL) / (dCL/dalpha); CL too low -> alpha up
            step = err * (alpha_ - alpha_old_) / (cl - cl_old_);
        }
        step = std::clamp(step, -s_.max_dalpha, s_.max_dalpha); // at most max_dalpha per change

        // Store alpha_old, CL_old, then alpha = alpha + dalpha.
        alpha_old_ = alpha_;
        cl_old_ = cl;
        have_old_ = true;
        alpha_ += step;
        ref_ = -1.0; // new reference from the next (jumped) residual on
        ++changes_;
        return true;
    }

}
