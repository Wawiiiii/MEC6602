#include "schemes.h"
#include "euler_study.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

void write_euler_results(const std::filesystem::path& filename,
                         const Eigen::MatrixXd& Q)
{
    // Same domain and gas constants as euler1d_mackcormack.
    constexpr double gamma = 1.4;
    constexpr double R = 287.0;
    const Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(Q.cols(), 0.0, 10.0);

    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open " + filename.string());
    }
    file << "# x A rho u p T Mach rhoA rhouA rhoEA\n"
         << "# Units: m m^2 kg/m^3 m/s Pa K - kg/m kg/s J/m\n"
         << std::scientific << std::setprecision(16);

    for (Eigen::Index i = 0; i < Q.cols(); ++i) {
        const double area = nozzle_area(x(i));
        const double rho = Q(0, i) / area;
        const double velocity = Q(1, i) / Q(0, i);
        const double pressure = (gamma - 1.0) *
            (Q(2, i) / area - 0.5 * rho * velocity * velocity);
        const double temperature = pressure / (rho * R);
        const double mach = velocity / std::sqrt(gamma * pressure / rho);

        file << x(i) << ' ' << area << ' ' << rho << ' ' << velocity << ' '
             << pressure << ' ' << temperature << ' ' << mach << ' '
             << Q(0, i) << ' ' << Q(1, i) << ' ' << Q(2, i) << '\n';
    }
    file.close();
    if (!file) {
        throw std::runtime_error("Could not write " + filename.string());
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::filesystem::path hw1_dir = HW1_SOURCE_DIR;
        const auto config = euler_study::read_config(hw1_dir / "euler1d_input.txt");
        const auto run_options = euler_study::parse_run_options(
            argc, argv, config.cfl_maccormack, 1.25, config.parallel_workers);
        const double dx = 10.0 / (config.n - 1);
        constexpr double u = 1.0; // The solver sets inlet velocity from Mach_in.
        const std::filesystem::path results_dir = hw1_dir / "Euler_1D_results";
        const auto scheme_dir = results_dir / "MacCormack";
        const auto profiles_dir = scheme_dir / "profiles";
        const auto residuals_dir = scheme_dir / "residuals";
        std::filesystem::create_directories(profiles_dir);
        std::filesystem::create_directories(residuals_dir);

        struct Case { double cfl; OutletType outlet; std::string outlet_name; std::string stem; };
        std::vector<Case> cases;
        for (double cfl : run_options.cfl_values) {
            for (OutletType outlet : {OutletType::Supersonic, OutletType::Subsonic}) {
                const std::string outlet_name = outlet == OutletType::Supersonic
                    ? "supersonic" : "subsonic";
                cases.push_back({cfl, outlet, outlet_name,
                    outlet_name + "_MacCormack_CFL_" + euler_study::cfl_label(cfl)});
            }
        }
        std::vector<std::string> messages(cases.size());
        std::vector<int> failed(cases.size(), 0);
        std::cout << "Running with " << euler_study::effective_worker_count(
            cases.size(), run_options.parallel_workers) << " thread(s).\n";
        euler_study::parallel_for(cases.size(), run_options.parallel_workers,
            [&](std::size_t index) {
                const auto& run = cases[index];
                const auto profile_path = profiles_dir / (run.stem + ".dat");
                const auto residual_path = residuals_dir / ("residual_" + run.stem + ".dat");
                std::filesystem::remove(profile_path);
                std::filesystem::remove(residual_path);
                std::vector<double> residuals;
                std::vector<IterationTiming> timings;
                bool solver_converged = false;
                std::ostringstream message;
                message << "MacCormack, " << run.outlet_name << " outlet, CFL = "
                        << run.cfl << ": ";
                try {
                    const Eigen::MatrixXd Q = euler1d_mackcormack(
                        run.cfl, u, dx, config.mach_in,
                        config.tolerance_maccormack, run.outlet,
                        config.back_pressure_ratio, &residuals,
                        config.max_iterations, &timings);
                    solver_converged = true;
                    write_euler_results(profile_path, Q);
                    message << "profile saved";
                } catch (const std::exception& e) {
                    std::error_code ignored;
                    std::filesystem::remove(profile_path, ignored);
                    message << "failed: " << e.what();
                    failed[index] = 1;
                }
                if (!timings.empty()) {
                    euler_study::write_iteration_history(residual_path, residuals, timings,
                                                         "relative update R^n", solver_converged);
                    message << ", iteration history saved";
                }
                messages[index] = message.str();
            });
        bool any_failed = false;
        for (std::size_t i = 0; i < cases.size(); ++i) {
            std::cout << messages[i] << '\n';
            any_failed = any_failed || failed[i] != 0;
        }
        return any_failed ? 1 : 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
