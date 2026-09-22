#include "schemes.h"
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

constexpr double PI = 3.14159265358979323846;

struct SimulationConfig {
    int n = 1201;
    double c = 0.5;
    std::vector<double> CFL_values = {0.25, 0.5, 0.75, 1.0, 1.25};
};

SimulationConfig read_simulation_config(const std::string& filename) {
    SimulationConfig config;
    std::ifstream file(filename);
    if (!file) {
        std::cout << "Warning: could not open " << filename
                  << ", using default parameters." << std::endl;
        return config;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        std::size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) {
            continue;
        }

        std::istringstream key_stream(line.substr(0, eq_pos));
        std::string key;
        key_stream >> key;
        if (key.empty()) {
            continue;
        }

        std::istringstream value_stream(line.substr(eq_pos + 1));

        if (key == "n") {
            value_stream >> config.n;
        } else if (key == "c") {
            value_stream >> config.c;
        } else if (key == "CFL") {
            config.CFL_values.clear();
            double value;
            while (value_stream >> value) {
                config.CFL_values.push_back(value);
            }
        }
    }

    return config;
}

int main() {
    SimulationConfig config = read_simulation_config("cfl_study.txt");

    int n = config.n;
    double xMin = 0.0;
    double xMax = PI;
    double c = config.c;

    Mesh mesh = make_mesh(n, xMin, xMax);
    InitialCondition ic = make_initial_condition(mesh, 0.5, 1.0);

    double x_mid_0 = (0.5 + 1.0) / 2.0;
    double x_mid_final = 2.5;
    double t_final = (x_mid_final - x_mid_0) / c;

    const std::vector<double>& CFL_values = config.CFL_values;

    std::map<std::string, std::function<Eigen::VectorXd(const Eigen::VectorXd&, double, int)>> schemes = {
        {"explicit_backward", explicit_backward},
        {"explicit_forward", explicit_forward},
        {"leap_frog", leap_frog},
        {"lax_wendroff", lax_wendroff},
        {"lax", lax},
        {"crank_nicolson", [](const Eigen::VectorXd& u0, double CFL, int nSteps) {
            return scheme_theta(u0, CFL, nSteps, 0.5);
        }},
        {"scheme_2space_4time", scheme_2space_4time},
        {"scheme_4space_2time", scheme_4space_2time},
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
