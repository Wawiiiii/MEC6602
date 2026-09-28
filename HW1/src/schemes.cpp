#include "schemes.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

Mesh make_mesh(int n, double xMin, double xMax) { 

    Mesh mesh; 
    mesh.n = n; 
    mesh.x = Eigen::VectorXd::LinSpaced(n, xMin, xMax); 
    mesh.dx = (xMax - xMin) / (n - 1); 
    return mesh; 

}

TimeParams make_time_params(double CFL, double dx, double c, double t_final) { 

    TimeParams tp; 
    double dt = CFL * dx / std::abs(c); 
    tp.nSteps = static_cast<int>(std::ceil(t_final / dt)); 
    tp.dt = t_final / tp.nSteps; 
    tp.CFL = c * tp.dt / dx;
    return tp;
}

InitialCondition make_initial_condition(const Mesh& mesh, double xStart, double xEnd) {
    
    InitialCondition ic;
    ic.u0 = Eigen::VectorXd::Zero(mesh.n);
    ic.u0 = ((mesh.x.array() >= xStart) && (mesh.x.array() <= xEnd)).select(1.0, ic.u0.array());
    return ic;
}

InitialCondition make_gaussian_initial_condition(const Mesh& mesh, double x0, double sigma) {
    InitialCondition ic;
    ic.u0 = (-((mesh.x.array() - x0) / sigma).square()).exp();
    return ic;
}

Eigen::VectorXd explicit_backward(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::VectorXd u = u0;
    Eigen::Index m = u.size() - 1;

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd u_new = u;
        u_new(0) = 0.0;
        u_new.tail(m) = (1.0 - CFL) * u.tail(m) + CFL * u.head(m);
        u = u_new;
    }

    return u;
}

Eigen::VectorXd explicit_forward(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::VectorXd u = u0;
    Eigen::Index m = u.size() - 1;

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd u_new = u;
        u_new.head(m) = (1.0 + CFL) * u.head(m) - CFL * u.tail(m);
        u = u_new;
    }

    return u;
}

Eigen::VectorXd leap_frog(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::Index m = u0.size() - 2;

    Eigen::VectorXd u_previous = u0;
    Eigen::VectorXd u_current = explicit_backward(u0, CFL, 1);

    for (int step = 1; step < nSteps; ++step) {
        Eigen::VectorXd u_next = u_current;
        u_next.segment(1, m) = u_previous.segment(1, m) - CFL * (u_current.tail(m) - u_current.head(m));
        u_previous = u_current;
        u_current = u_next;
    }

    return u_current;
}

Eigen::VectorXd lax_wendroff(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::Index m = u0.size() - 2;
    Eigen::VectorXd u_current = u0;

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd u_next = u_current;

        u_next.segment(1, m) = u_current.segment(1, m) - 0.5 * CFL * (u_current.tail(m) - u_current.head(m)) 
        + 0.5 * CFL * CFL * (u_current.tail(m) - 2 * u_current.segment(1, m) + u_current.head(m));

        u_current = u_next;
    }

    return u_current;
}

Eigen::VectorXd lax(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::Index m = u0.size() - 2;
    Eigen::VectorXd u_current = u0;

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd u_next = u_current;
        u_next.segment(1, m) = 0.5 * (u_current.tail(m) + u_current.head(m)) - 0.5 * CFL * (u_current.tail(m) - u_current.head(m));
        u_current = u_next;
    }

    return u_current;
}

Eigen::VectorXd scheme_2space_4time(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::Index m = u0.size() - 4;
    Eigen::VectorXd u_current = u0;

    const double c1 = CFL / 2.0;
    const double c2 = CFL * CFL / 2.0;
    const double c3 = CFL * CFL * CFL / 12.0;
    const double c4 = CFL * CFL * CFL * CFL / 24.0;

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd u_next = u_current;

        u_next.segment(2, m) = u_current.segment(2, m)
            - c1 * (u_current.segment(3, m) - u_current.segment(1, m))
            + c2 * (u_current.segment(3, m)
                    - 2.0 * u_current.segment(2, m)
                    + u_current.segment(1, m))
            - c3 * (u_current.tail(m)
                    - 2.0 * u_current.segment(3, m)
                    + 2.0 * u_current.segment(1, m)
                    - u_current.head(m))
            + c4 * (u_current.head(m)
                    - 4.0 * u_current.segment(1, m)
                    + 6.0 * u_current.segment(2, m)
                    - 4.0 * u_current.segment(3, m)
                    + u_current.tail(m));

        u_current = u_next;
    }

    return u_current;
}

Eigen::VectorXd scheme_4space_2time(const Eigen::VectorXd& u0, double CFL, int nSteps) {

    Eigen::Index m = u0.size() - 4;
    Eigen::VectorXd u_current = u0;

    const double c1 = CFL / 12.0;
    const double c2 = CFL * CFL / 24.0;

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd u_next = u_current;

        u_next.segment(2, m) = u_current.segment(2, m)
            - c1 * (u_current.head(m)
                    - 8.0 * u_current.segment(1, m)
                    + 8.0 * u_current.segment(3, m)
                    - u_current.tail(m))
            + c2 * (-u_current.tail(m)
                    + 16.0 * u_current.segment(3, m)
                    - 30.0 * u_current.segment(2, m)
                    + 16.0 * u_current.segment(1, m)
                    - u_current.head(m));

        u_current = u_next;
    }

    return u_current;
}

Eigen::VectorXd scheme_theta(const Eigen::VectorXd& u0, double CFL,
                             int nSteps, double theta) {

    Eigen::Index n = u0.size();
    Eigen::Index m = n - 2;
    Eigen::VectorXd u_current = u0;

    const double c_implicit = theta * CFL / 2.0;
    const double c_explicit = (1.0 - theta) * CFL / 2.0;

    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n);
    A.diagonal(1).tail(m).setConstant(c_implicit);
    A.diagonal(-1).head(m).setConstant(-c_implicit);

    Eigen::PartialPivLU<Eigen::MatrixXd> solver(A);

    for (int step = 0; step < nSteps; ++step) {
        Eigen::VectorXd b = u_current;

        b.segment(1, m) = u_current.segment(1, m)
            - c_explicit * (u_current.tail(m) - u_current.head(m));

        u_current = solver.solve(b);
    }

    return u_current;
}


double nozzle_area(double x)
{
    return 1.398 + 0.347 * std::tanh(0.8*x - 4);
}

double der_nozzle_area(double x)
{
    const double sech = 1.0 / std::cosh(0.8*x - 4);
    return 0.2776 * sech * sech;
}



void apply_boundary_conditions(
    Eigen::MatrixXd& Q,
    const Eigen::VectorXd& A,
    double gamma,
    double R,
    double T_in,
    double P_in,
    double Mach_in,
    OutletType outlet_type,
    double back_pressure_ratio)
{
    // ----- Inlet supersonic -----
    double c_in = std::sqrt(gamma * R * T_in);
    double u_in = Mach_in * c_in;
    double rho_in = P_in / (R * T_in);
    double e_in = R * T_in / (gamma - 1.0);
    double E_in = e_in + 0.5 * u_in * u_in;

    Q(0, 0) = rho_in * A(0);
    Q(1, 0) = rho_in * u_in * A(0);
    Q(2, 0) = rho_in * E_in * A(0);

    // ----- Outlet -----
    int N = Q.cols() - 1;

    if (outlet_type == OutletType::Supersonic)
    {
        // Extrapolation primitive
        Q.col(N) = Q.col(N - 1) * A(N) / A(N - 1);
    }
    else
    {
        // intérieur
        double rho_i = Q(0, N - 1) / A(N - 1);
        double u_i   = Q(1, N - 1) / Q(0, N - 1);
        double E_i   = Q(2, N - 1) / Q(0, N - 1);

        double P_i =
            (gamma - 1.0) * rho_i *
            (E_i - 0.5 * u_i * u_i);

        double c_i = std::sqrt(gamma * P_i / rho_i);

        // The slides impose the downstream static pressure at a subsonic exit.
        // Here back_pressure_ratio means P_B / P_in (static inlet pressure).
        const double P0_in = P_in * std::pow(
            1.0 + 0.5 * (gamma - 1.0) * Mach_in * Mach_in,
            gamma / (gamma - 1.0));
        const double P_L = back_pressure_ratio * P_in;
        if (!(P_L > 0.0) || !std::isfinite(P_L) || P_L >= P0_in)
            throw std::invalid_argument(
                "Subsonic-outlet back pressure must satisfy 0 < P_B < inlet total pressure");

        // Hypothèse isentropique à la sortie
        double rho_L =
            rho_i * std::pow(P_L / P_i, 1.0 / gamma);

        double c_L =
            std::sqrt(gamma * P_L / rho_L);

        // Invariant venant de l'intérieur
        double R1 = u_i + 2.0 * c_i / (gamma - 1.0);

        // Reconstruction
        double u_L =
            R1 - 2.0 * c_L / (gamma - 1.0);

        double e_L = P_L / ((gamma - 1.0) * rho_L);
        double E_L = e_L + 0.5 * u_L * u_L;

        Q(0, N) = rho_L * A(N);
        Q(1, N) = rho_L * u_L * A(N);
        Q(2, N) = rho_L * E_L * A(N);
    }
}



Eigen::MatrixXd euler1d_mackcormack(double CFL, double u, double dx, double Mach, double convergence, OutletType outlet_type, double back_pressure_ratio)
{
    if (CFL <= 0.0 || CFL >= 1.0 ||
        dx <= 0.0 ||
        Mach <= 0.0 ||
        convergence <= 0.0 ||
        dx > 5.0)
    {
        throw std::invalid_argument(
            "Require 0 < CFL < 1, dx > 0, Mach > 0 and convergence > 0"
        );
    }


    (void)u; // Kept in the public signature; inlet velocity is set by Mach below.
    const double gamma = 1.4;
    const double R = 287.0;
    const double T = 300.0;
    const double P = 101325.0;
    const double rho = P / (R * T);
    const double e = R * T / (gamma - 1.0);
    const double c = std::sqrt(gamma * R * T);
    const double u_in = Mach * c;
    const double E = e + 0.5 * u_in * u_in;

    int n = static_cast<int>(10.0 / dx) + 1;
    const double grid_dx = 10.0 / (n - 1);
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0, 10);
    Eigen::VectorXd A = x.unaryExpr([](double xi) { return nozzle_area(xi); });
    Eigen::VectorXd dAdx = x.unaryExpr([](double xi) { return der_nozzle_area(xi); });

    Eigen::MatrixXd Q_prev(3, n);
    Q_prev.row(0) = (rho * A).transpose();
    Q_prev.row(1) = (rho * u_in * A).transpose();
    Q_prev.row(2) = (rho * E * A).transpose();

    auto flux_and_source = [&](const Eigen::MatrixXd& Q, Eigen::MatrixXd& flux,
                               Eigen::MatrixXd& source) {
        double max_speed = 0.0;
        for (int i = 0; i < n; ++i) {
            const double density = Q(0, i) / A(i);
            const double velocity = Q(1, i) / Q(0, i);
            const double pressure = (gamma - 1.0) *
                (Q(2, i) / A(i) - 0.5 * density * velocity * velocity);
            if (!std::isfinite(pressure) || pressure <= 0.0 || density <= 0.0)
                throw std::runtime_error("MacCormack step produced a nonphysical state");
            max_speed = std::max(max_speed, std::abs(velocity) + std::sqrt(gamma * pressure / density));
            flux(0, i) = Q(1, i);
            flux(1, i) = Q(1, i) * velocity + pressure * A(i);
            flux(2, i) = velocity * (Q(2, i) + pressure * A(i));
            source(0, i) = 0.0;
            source(1, i) = pressure * dAdx(i);
            source(2, i) = 0.0;
        }
        return max_speed;
    };

    Eigen::MatrixXd flux_prev(3, n), source_prev(3, n);
    Eigen::MatrixXd flux_pred(3, n), source_pred(3, n);
    for (int step = 0; step < 100000; ++step) {
        const double max_speed = flux_and_source(Q_prev, flux_prev, source_prev);
        const double dt = CFL * grid_dx / max_speed;

        Eigen::MatrixXd Q_pred = Q_prev;
        for (int i = 1; i < n - 1; ++i)
            Q_pred.col(i) = Q_prev.col(i) - dt / grid_dx *
                (flux_prev.col(i) - flux_prev.col(i-1)) + dt * source_prev.col(i);

        apply_boundary_conditions(Q_pred, A, gamma, R, T, P, Mach, outlet_type, back_pressure_ratio);

        flux_and_source(Q_pred, flux_pred, source_pred);
        Eigen::MatrixXd Q_new = Q_prev;
        for (int i = 1; i < n - 1; ++i)
            Q_new.col(i) = 0.5 * (Q_prev.col(i) + Q_pred.col(i) - dt / grid_dx *
                (flux_pred.col(i + 1) - flux_pred.col(i)) + dt * source_pred.col(i));

        apply_boundary_conditions(Q_new, A, gamma, R, T, P, Mach, outlet_type, back_pressure_ratio);

        const double residual = (Q_new - Q_prev).norm() / Q_prev.norm();
        if (!std::isfinite(residual))
            throw std::runtime_error("MacCormack residual is not finite");
        if (residual < convergence)
            return Q_new;
        Q_prev.swap(Q_new);
    }
    throw std::runtime_error("MacCormack solver did not converge within 100000 steps");
}

Eigen::MatrixXd euler1d_implicit(double CFL, double u, double dx, double Mach,
                                 double convergence, OutletType outlet_type,
                                 double back_pressure_ratio)
{
    if (CFL <= 0.0 || CFL >= 1.0 || dx <= 0.0 || Mach <= 0.0 ||
        convergence <= 0.0 || dx > 5.0)
        throw std::invalid_argument("Require 0 < CFL < 1, dx > 0, Mach > 0 and convergence > 0");
    (void)u;

    constexpr double gamma = 1.4, R = 287.0, T = 300.0, P = 101325.0;
    constexpr double theta = 1.0;
    constexpr double eps_explicit = 0.02; // epsilon_e < 0.125 (course slides)
    constexpr double eps_implicit = 0.05; // epsilon_i > 2 epsilon_e
    const double c = std::sqrt(gamma * R * T);
    const double velocity = Mach * c;
    const double rho = P / (R * T);
    const double energy = R * T / (gamma - 1.0) + 0.5 * velocity * velocity;
    const int n = static_cast<int>(10.0 / dx) + 1;
    const double grid_dx = 10.0 / (n - 1);
    const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0.0, 10.0);
    const Eigen::VectorXd area = x.unaryExpr([](double xi) { return nozzle_area(xi); });
    const Eigen::VectorXd darea = x.unaryExpr([](double xi) { return der_nozzle_area(xi); });

    Eigen::MatrixXd Q(3, n);
    for (int i = 0; i < n; ++i) Q.col(i) << rho * area(i), rho * velocity * area(i), rho * energy * area(i);
    apply_boundary_conditions(Q, area, gamma, R, T, P, Mach, outlet_type, back_pressure_ratio);

    using Block = Eigen::Matrix3d;
    const int ni = n - 2;
    auto flux_source_jacobian = [&](int i, Eigen::Vector3d& flux, Eigen::Vector3d& source, Block& jac) {
        const double r = Q(0,i) / area(i), v = Q(1,i) / Q(0,i), E = Q(2,i) / Q(0,i);
        const double p = (gamma - 1.0) * r * (E - 0.5*v*v);
        if (!(r > 0.0 && p > 0.0) || !std::isfinite(p))
            throw std::runtime_error("Implicit Beam-Warming step produced a nonphysical state");
        const double H = (Q(2,i) + p*area(i)) / Q(0,i);
        flux << Q(1,i), Q(1,i)*v + p*area(i), v*(Q(2,i) + p*area(i));
        source << 0.0, p*darea(i), 0.0;
        jac << 0.0, 1.0, 0.0,
               0.5*(gamma-3.0)*v*v, (3.0-gamma)*v, gamma-1.0,
               v*(0.5*(gamma-1.0)*v*v-H), H-(gamma-1.0)*v*v, gamma*v;
    };

    for (int step = 0; step < 100000; ++step) {
        double max_speed = 0.0;
        for (int i = 0; i < n; ++i) {
            const double r = Q(0,i)/area(i), v = Q(1,i)/Q(0,i);
            const double p = (gamma-1.0)*r*(Q(2,i)/Q(0,i)-0.5*v*v);
            if (r <= 0.0 || p <= 0.0 || !std::isfinite(p)) throw std::runtime_error("Nonphysical Euler state");
            max_speed = std::max(max_speed, std::abs(v)+std::sqrt(gamma*p/r));
        }
        const double dt = CFL * grid_dx / max_speed;
        std::vector<Block> lower(ni, Block::Zero()), diagonal(ni), upper(ni, Block::Zero());
        std::vector<Eigen::Vector3d> rhs(ni);
        double rhs_norm_squared = 0.0;
        for (int j = 0; j < ni; ++j) {
            const int i = j+1;
            Eigen::Vector3d Fm, Fi, Fp, Sm, Si, Sp;
            Block Jm, Ji, Jp;
            flux_source_jacobian(i-1,Fm,Sm,Jm); flux_source_jacobian(i,Fi,Si,Ji); flux_source_jacobian(i+1,Fp,Sp,Jp);
            const Eigen::Vector3d residual = (Fp-Fm)/(2.0*grid_dx) - Si;
            Eigen::Vector3d fourth_difference = Eigen::Vector3d::Zero();
            if (i >= 2 && i + 2 < n)
                fourth_difference = Q.col(i-2) - 4.0*Q.col(i-1) + 6.0*Q.col(i)
                    - 4.0*Q.col(i+1) + Q.col(i+2);
            // Convert epsilon_e to a rate using the cell crossing time dx/lambda,
            // then multiply by dt as required by the course formula.
            rhs[j] = -dt*residual
                - eps_explicit*(dt*max_speed/grid_dx)*fourth_difference;
            rhs_norm_squared += rhs[j].squaredNorm();
            diagonal[j] = Block::Identity() + 2.0*eps_implicit*Block::Identity();
            lower[j] = -eps_implicit*Block::Identity(); upper[j] = -eps_implicit*Block::Identity();
            upper[j] += theta*dt/(2.0*grid_dx)*Ji;
            lower[j] -= theta*dt/(2.0*grid_dx)*Ji;
        }
        // The correction can be artificially small after positivity damping.
        // Test the discrete Beam-Warming equation itself before declaring convergence.
        if (std::sqrt(rhs_norm_squared) / Q.norm() < convergence) return Q;
        // Fixed boundary increments are zero; remove their couplings from the system.
        lower[0].setZero(); upper[ni-1].setZero();
        std::vector<Block> cprime(ni);
        std::vector<Eigen::Vector3d> dprime(ni), correction(ni);
        Eigen::PartialPivLU<Block> first_solver(diagonal[0]);
        cprime[0] = first_solver.solve(upper[0]);
        dprime[0] = first_solver.solve(rhs[0]);
        for (int j = 1; j < ni; ++j) {
            const Block pivot = diagonal[j] - lower[j]*cprime[j-1];
            Eigen::PartialPivLU<Block> solver(pivot);
            if (j + 1 < ni)
                cprime[j] = solver.solve(upper[j]);
            else
                cprime[j].setZero();
            dprime[j] = solver.solve(rhs[j]-lower[j]*dprime[j-1]);
        }
        correction[ni-1] = dprime[ni-1];
        for (int j = ni-2; j >= 0; --j) correction[j] = dprime[j]-cprime[j]*correction[j+1];
        Eigen::MatrixXd Qnew;
        double relaxation = 1.0;
        bool physical = false;
        Eigen::Index bad_cell = -1;
        double bad_density = 0.0, bad_pressure = 0.0;
        for (int attempt = 0; attempt < 60; ++attempt) {
            Qnew = Q;
            for (int j = 0; j < ni; ++j) Qnew.col(j+1) += relaxation*correction[j];
            apply_boundary_conditions(Qnew, area, gamma, R, T, P, Mach, outlet_type, back_pressure_ratio);
            physical = true;
            for (int i = 1; i < n; ++i) {
                const double r = Qnew(0,i)/area(i), v = Qnew(1,i)/Qnew(0,i);
                const double p = (gamma-1.0)*r*(Qnew(2,i)/Qnew(0,i)-0.5*v*v);
                if (!(r > 0.0 && p > 0.0) || !std::isfinite(r) || !std::isfinite(v) || !std::isfinite(p)) {
                    physical = false;
                    bad_cell = i;
                    bad_density = r;
                    bad_pressure = p;
                    break;
                }
            }
            if (physical) break;
            relaxation *= 0.5;
        }
        if (!physical) {
            throw std::runtime_error("Beam-Warming correction could not preserve a physical state at iteration "
                + std::to_string(step) + ", cell " + std::to_string(bad_cell)
                + " (rho=" + std::to_string(bad_density)
                + ", p=" + std::to_string(bad_pressure) + ")");
        }
        const double update = (Qnew-Q).norm()/Q.norm();
        if (!std::isfinite(update)) throw std::runtime_error("Implicit update is not finite");
        Q.swap(Qnew);
    }
    throw std::runtime_error("Implicit Beam-Warming solver did not converge within 100000 steps");
}
