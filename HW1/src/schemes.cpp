#include "schemes.h"
#include <cmath>

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


Eigen::VectorXd nozzle_area(Eigen::VectorXd x);
{
    return 1.398 + 0.347 * std::tanh(0.8*x - 4);
}

Eigen::VectorXd der_nozzle_area(Eigen::VectorXd x);
{
    return 0.2776 * (1/std::cosh(0.8*x - 4)) * (1/std::cosh(0.8*x - 4));
}

Eigen::MatrixXd euler1d_mackcormack(double CFL, double u, double dx, double Mach=1.25, double convergence = 1e-6)
{
    const double gamma = 1.4;
    const double R = 287.0;
    const double T = 300.0;
    const double P = 101325.0;
    const double rho = P / (R * T);
    const double e = R * T / (gamma - 1.0);
    const double E = e + 0.5 * u_in * u_in;
    const double H = E + P / rho;
    const double c = std::sqrt(gamma * R * T);
    const double u_in = Mach * c;

    bool conv = false;
    double residual = 1.0;
    int n = static_cast<int>(10.0 / dx) + 1;
    Eigen::MatrixXd Q_prev(3, n);
    Eigen::MatrixXd Q_new(3, n);
    Eigen::MatrixXd E_prev(3,n);
    Eigen::MatrixXd E_new(3,n);
    Eigen::MatrixXd S_prev(3,n);
    Eigen::MatrixXd S_new(3,n);

   
    Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, 0, 10);

    Eigen::VectorXd A = nozzle_area(x);
    Eigen::VectorXd dAdx = der_nozzle_area(x);

    Q_prev(0,Eigen::all) = (rho*A).transpose();
    Q_prev(1,Eigen::all) = (rho*u_in*A).transpose();
    Q_prev(2, Eigen::all) = (rho * E * A).transpose();

    E_prev(0,Eigen::all) = (rho*u_in*A).transpose();
    E_prev(1,Eigen::all) = ((rho*u_in*u_in + P)*A).transpose();
    E_prev(2,Eigen::all) = (rho*u_in*H*A).transpose();
    
    S_prev(0,Eigen::all) = Eigen::VectorXd::Zero(n);
    S_prev(1,Eigen::all) = (P*dAdx).transpose();
    S_prev(2,Eigen::all) = Eigen::VectorXd::Zero(n);

    while (!conv)
    {



        residual = (Q_new - Q_prev).norm();

        if (residual < convergence)
            conv = true;

        Q_prev = Q_new;
    }

    return Q_new;
}