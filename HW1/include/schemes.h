#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <vector>

struct Mesh
{
    Eigen::VectorXd x;
    double dx;
    int n;
};

struct TimeParams
{
    double dt;
    double CFL;
    int nSteps;
};

Mesh make_mesh(int n, double xMin, double xMax);

double nozzle_area(double x);

double der_nozzle_area(double x);

TimeParams make_time_params(double CFL, double dx, double c, double t_final);

struct InitialCondition
{
    Eigen::VectorXd u0;
};

InitialCondition make_initial_condition(const Mesh &mesh, double xStart, double xEnd);

// Smooth (infinitely differentiable) pulse, for convergence studies where a
// discontinuous initial condition would cap the observed order of accuracy.
InitialCondition make_gaussian_initial_condition(const Mesh &mesh, double x0, double sigma);

Eigen::VectorXd explicit_backward(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd explicit_forward(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd leap_frog(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd lax_wendroff(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd lax(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd scheme_2space_4time(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd scheme_4space_2time(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd scheme_theta(const Eigen::VectorXd &u0, double CFL, int nSteps, double theta = 0.5);

Eigen::VectorXd euler1d_mackcormack(double CFL, double u, double dx, double Mach=1.25, double convergence = 1e-6);