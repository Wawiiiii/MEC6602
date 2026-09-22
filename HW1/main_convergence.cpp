#include "schemes.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

constexpr double PI = 3.14159265358979323846;

struct ConvergenceConfig {
    double c = 0.5;

    // initial condition used for the error measurement:
    //   "pulse"    -> the discontinuous top-hat from the homework (caps the
    //                 observed order well below the scheme's nominal order)
    //   "gaussian" -> a smooth bump, needed to actually observe the nominal
    //                 spatial order of accuracy
    std::string ic = "gaussian";
    double gaussian_x0 = 0.75;
    double gaussian_sigma = 0.1;

    // fixed dt, sweep over n (mesh points) to isolate spatial error
    double dt = 0.0005;
    std::vector<int> n_values = {51, 101, 201, 401, 801, 1601};
};

ConvergenceConfig read_convergence_config(const std::string& filename) {
    ConvergenceConfig config;
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

        if (key == "ic") {
            value_stream >> config.ic;
        } else if (key == "gaussian_x0") {
            value_stream >> config.gaussian_x0;
        } else if (key == "gaussian_sigma") {
            value_stream >> config.gaussian_sigma;
        } else if (key == "c") {
            value_stream >> config.c;
        } else if (key == "dt") {
            value_stream >> config.dt;
        } else if (key == "n_values") {
            config.n_values.clear();
            int value;
            while (value_stream >> value) {
                config.n_values.push_back(value);
            }
        }
    }

    return config;
}

double l2_error(const Eigen::VectorXd& u, const Eigen::VectorXd& u_exact, double dx) {
    return std::sqrt(dx * (u - u_exact).array().square().sum());
}

// Builds the initial condition and the exact solution shifted to time t,
// according to config.ic ("pulse" or "gaussian").
InitialCondition make_ic(const ConvergenceConfig& config, const Mesh& mesh) {
    if (config.ic == "gaussian") {
        return make_gaussian_initial_condition(mesh, config.gaussian_x0, config.gaussian_sigma);
    }
    return make_initial_condition(mesh, 0.5, 1.0);
}

InitialCondition make_exact(const ConvergenceConfig& config, const Mesh& mesh, double c, double t) {
    if (config.ic == "gaussian") {
        return make_gaussian_initial_condition(mesh, config.gaussian_x0 + c * t, config.gaussian_sigma);
    }
    return make_initial_condition(mesh, 0.5 + c * t, 1.0 + c * t);
}

// Reports one scheme's (dx, error, order) table, both to stdout and to a .dat file.
// std::ostream prints NaN/Inf differently per platform (e.g. MSVC writes
// "-nan(ind)", which numpy.loadtxt can't parse) -- normalize to the
// portable "nan"/"inf"/"-inf" text that every reader expects.
std::string format_double(double value) {
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "inf" : "-inf";
    }
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << value;
    return oss.str();
}

void report_convergence(const std::vector<double>& dx_values,
                         const std::map<std::string, std::vector<double>>& errors) {
    for (const auto& [name, error_values] : errors) {
        std::vector<double> order(error_values.size(), std::numeric_limits<double>::quiet_NaN());
        for (std::size_t i = 1; i < error_values.size(); ++i) {
            order[i] = std::log(error_values[i - 1] / error_values[i]) /
                       std::log(dx_values[i - 1] / dx_values[i]);
        }

        std::cout << "\n=== " << name << " (space study) ===\n";
        std::cout << std::left
                   << std::setw(16) << "dx"
                   << std::setw(16) << "L2 error"
                   << std::setw(10) << "order" << "\n";

        std::string filename = "results_convergence/space_" + name + ".dat";
        std::ofstream file(filename);
        file << "# dx error order\n";

        for (std::size_t i = 0; i < dx_values.size(); ++i) {
            std::string order_text = "-";
            std::string order_file_text = "nan";
            if (std::isnan(order[i])) {
                // both stay "-" / "nan"
            } else if (std::isinf(order[i])) {
                order_text = order_file_text = (order[i] > 0 ? "inf" : "-inf");
            } else {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << order[i];
                order_text = oss.str();
                order_file_text = order_text;
            }

            std::string dx_text = format_double(dx_values[i]);
            std::string error_text = format_double(error_values[i]);

            std::cout << std::left
                       << std::setw(16) << dx_text
                       << std::setw(16) << error_text
                       << std::setw(10) << order_text << "\n";

            file << dx_text << " " << error_text << " " << order_file_text << "\n";
        }

        file.close();
    }
    std::cout << std::endl;
}

int main() {
    ConvergenceConfig config = read_convergence_config("convergence_study.txt");

    double c = config.c;
    double x_mid_0 = (config.ic == "gaussian") ? config.gaussian_x0 : (0.5 + 1.0) / 2.0;
    double x_mid_final = 2.5;
    double t_final = (x_mid_final - x_mid_0) / c;

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

    std::filesystem::create_directories("results_convergence");

    std::vector<double> dx_values;
    std::map<std::string, std::vector<double>> errors;

    for (int n : config.n_values) {
        Mesh mesh = make_mesh(n, 0.0, PI);
        InitialCondition ic = make_ic(config, mesh);

        double dt = config.dt;
        int nSteps = static_cast<int>(std::round(t_final / dt));
        double t_actual = nSteps * dt;
        double CFL = c * dt / mesh.dx;

        InitialCondition exact_ic = make_exact(config, mesh, c, t_actual);

        if (CFL >= 1.0) {
            std::cout << "Warning: CFL = " << CFL << " >= 1 at n = " << n
                      << " (dx = " << mesh.dx << ", dt = " << dt
                      << ") -- conditionally-stable schemes will diverge here.\n";
        }

        dx_values.push_back(mesh.dx);

        for (const auto& [name, scheme] : schemes) {
            Eigen::VectorXd u = scheme(ic.u0, CFL, nSteps);
            double error = l2_error(u, exact_ic.u0, mesh.dx);
            errors[name].push_back(error);
        }
    }

    report_convergence(dx_values, errors);

    return 0;
}
