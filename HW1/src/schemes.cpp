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