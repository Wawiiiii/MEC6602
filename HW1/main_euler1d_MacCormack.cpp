#include "schemes.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

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
    std::cout << "Results saved to " << filename << '\n';
}

int main()
{
    try
    {
        const double CFL = 0.5;
        const double dx = 0.05;
        const double Mach = 1.25;
        const double convergence = 1e-8;
        const double u = 1;
        // Like results/, this directory is relative to the working directory.
        const std::filesystem::path results_dir = "Euler_1D_results";
        std::filesystem::create_directories(results_dir);

        // -----------------------------
        // Cas 1 : sortie supersonique
        // -----------------------------
        std::cout << "Running supersonic outlet case...\n";

        Eigen::MatrixXd Q_supersonic =
            euler1d_mackcormack(
                CFL,
                u,
                dx,
                Mach,
                convergence,
                OutletType::Supersonic
            );

        std::cout << "Supersonic case converged.\n";
        write_euler_results(results_dir / "supersonic_MacCormack.dat", Q_supersonic);
        std::cout << "Q size = "
                  << Q_supersonic.rows()
                  << " x "
                  << Q_supersonic.cols()
                  << "\n\n";

        // -----------------------------
        // Cas 2 : sortie subsonique
        // -----------------------------
        std::cout << "Running subsonic outlet case...\n";

        const double back_pressure_ratio = 1.9;

        Eigen::MatrixXd Q_subsonic =
            euler1d_mackcormack(
                CFL,
                u,
                dx,
                Mach,
                convergence,
                OutletType::Subsonic,
                back_pressure_ratio
            );

        std::cout << "Subsonic case converged.\n";
        write_euler_results(results_dir / "subsonic_MacCormack.dat", Q_subsonic);
        std::cout << "Q size = "
                  << Q_subsonic.rows()
                  << " x "
                  << Q_subsonic.cols()
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
