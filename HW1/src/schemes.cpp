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

        // Pression imposée
        double P_L = back_pressure_ratio * P_in;

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
