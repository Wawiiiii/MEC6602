#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <vector>

// Uniform one-dimensional mesh and its spacing.
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

enum class OutletType
{
    Supersonic,
    Subsonic
};

// Cumulative solver time and duration of one iteration, in seconds.
struct IterationTiming
{
    double elapsed_seconds;
    double iteration_seconds;
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

// Linear advection schemes returning the state after nSteps.
Eigen::VectorXd explicit_backward(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd explicit_forward(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd leap_frog(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd lax_wendroff(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd lax(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd scheme_2space_4time(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd scheme_4space_2time(const Eigen::VectorXd &u0, double CFL, int nSteps);
Eigen::VectorXd scheme_theta(const Eigen::VectorXd &u0, double CFL, int nSteps, double theta = 0.5);


// back_pressure_ratio is P_B / P_in, with P_in the static inlet pressure.
void apply_boundary_conditions(Eigen::MatrixXd& Q, const Eigen::VectorXd& A, double gamma, double R, double T_in, double P_in, double Mach_in, OutletType outlet_type, double back_pressure_ratio = 1.9);
// Steady nozzle solvers returning rows rho*A, rho*u*A and rho*E*A.
Eigen::MatrixXd euler1d_mackcormack(double CFL, double u, double dx, double Mach,
                                    double convergence, OutletType outlet_type,
                                    double back_pressure_ratio = 1.9,
                                    std::vector<double>* residual_history = nullptr,
                                    int max_iterations = 100000,
                                    std::vector<IterationTiming>* timing_history = nullptr);
// residual_history stores R^n = ||Q^(n+1)-Q^n||_2 / ||Q^n||_2 after the
// accepted correction and boundary conditions. Convergence still uses the
// pre-correction relative discrete equation residual.
Eigen::MatrixXd euler1d_implicit(double CFL, double u, double dx, double Mach,
                                 double convergence, OutletType outlet_type,
                                 double back_pressure_ratio = 1.9,
                                 std::vector<double>* residual_history = nullptr,
                                 int max_iterations = 100000,
                                 std::vector<IterationTiming>* timing_history = nullptr);
