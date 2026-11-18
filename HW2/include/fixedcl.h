#pragma once

namespace euler {

    /*
    Lift-coefficient targeting, iso-CL computation (post.pdf slide 9, Part 5 of the homework).

    The angle of attack becomes an unknown: the flow is iterated at the current alpha and, every time
    the flow is "converged" enough, alpha is corrected by a secant step on CL(alpha):
      1-3. start at the given alpha, iterate until the residual has dropped by `drop`
      4-5. store alpha_old, CL_old, perturb alpha by a small amount (`dalpha0`)
      6.   iterate again until the residual has dropped by `drop`
      7.   secant step  dalpha = (CL_t - CL) (alpha - alpha_old) / (CL - CL_old)
           store alpha_old, CL_old, alpha = alpha + dalpha
      8.   continue: the flow and alpha converge together; the new alpha enters through the farfield
           BC (the free stream direction). Slide: "do not fully converge before changing alpha".

    Sign of the secant step: the slide writes (CL - CL_t) in the numerator. With alpha = alpha + dalpha
    that steps AWAY from the target (CL below target, CL increasing with alpha -> dalpha < 0), so the
    numerator is used as (CL_t - CL), the standard secant step (DEVIATIONS.md A29).

    "Converged enough" (the slide does not give a number): the residual norm is `drop` times the
    largest residual seen since the last alpha change (the residual jumps after each change, slide 10).
    Once |CL - CL_t| <= cl_tol at such a point, alpha is left alone and the flow converges normally.
    The steps are clamped to +-max_dalpha (protects the first secant step when CL barely moved).
    */
    struct FixedClSettings {
        bool on = false;          // fixed-CL mode (off = fixed alpha, the normal run)
        double cl_target = 0.4;   // target lift coefficient CL_t
        double dalpha0 = 0.1;     // first perturbation of alpha, degrees (slide 9 step 5)
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
        Returns true when alpha was changed (the caller rebuilds the free stream from alpha()).
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
        bool have_old_ = false;   // false until the first perturbation (slide 9 steps 4-5)
        double ref_ = -1.0;       // largest residual since the last change (-1: none yet)
        int changes_ = 0;         // alpha changes so far
    };

}
