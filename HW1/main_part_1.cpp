#include "schemes.h"
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

constexpr double PI = 3.14159265358979323846;

int main() {
    int n = 1201;
    double xMin = 0.0;
    double xMax = PI;
    double c = 1;

    Mesh mesh = make_mesh(n, xMin, xMax);
    InitialCondition ic = make_initial_condition(mesh, 0.5, 1.0);

    double x_mid_0 = (0.5 + 1.0) / 2.0;
    double x_mid_final = 2.5;
    double t_final = (x_mid_final - x_mid_0) / c;

    double CFL_values[] = {0.1, 0.2,0,25, 0.3, 0.4, 0.5, 0.6, 0.7, 0,75, 0.8, 0.9, 1.0, 1.25};

    std::map<std::string, std::function<Eigen::VectorXd(const Eigen::VectorXd&, double, int)>> schemes = {
        {"explicit_backward", explicit_backward},
        {"explicit_forward", explicit_forward},
        {"leap_frog", leap_frog},
        {"lax_wendroff", lax_wendroff},
        {"lax", lax},
    };

    std::filesystem::create_directories("results");

    for (const auto& [name, scheme] : schemes) {
        for (double CFL : CFL_values) {
            TimeParams tp = make_time_params(CFL, mesh.dx, c, t_final);

            Eigen::VectorXd u = scheme(ic.u0, tp.CFL, tp.nSteps);

            std::ostringstream cfl_text;
            cfl_text << std::fixed << std::setprecision(2) << CFL;

            std::string filename = "results/" + name + "_CFL_" + cfl_text.str() + ".dat";
            std::ofstream file(filename);
            file << "# x u\n";
            for (int i = 0; i < mesh.n; ++i) {
                file << mesh.x(i) << " " << u(i) << "\n";
            }
            file.close();

            std::cout << name << ": CFL demande = " << CFL
                      << ", CFL reel = " << tp.CFL
                      << ", dt = " << tp.dt
                      << ", nSteps = " << tp.nSteps
                      << " -> " << filename << std::endl;
        }
    }

    return 0;
}
