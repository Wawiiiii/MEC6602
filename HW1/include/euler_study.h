#pragma once

#include <cmath>
#include <algorithm>
#include <atomic>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace euler_study {

inline std::string cfl_label(double cfl) {
    std::ostringstream text;
    text << std::fixed << std::setprecision(6) << cfl;
    std::string label = text.str();
    while (label.back() == '0') label.pop_back();
    if (label.back() == '.') label.pop_back();
    return label;
}

struct Config {
    int n = 201;
    int max_iterations = 100000;
    int parallel_workers = 0;
    double mach_in = 1.25;
    double back_pressure_ratio = 1.9;
    double tolerance_maccormack = 1e-8;
    double tolerance_implicit = 1e-8;
    std::vector<double> cfl_maccormack;
    std::vector<double> cfl_implicit;
};

template <typename T>
inline T parse_scalar(const std::string& text, const std::string& key) {
    std::istringstream stream(text);
    T value;
    if (!(stream >> value)) throw std::invalid_argument("Invalid value for " + key);
    stream >> std::ws;
    if (!stream.eof()) throw std::invalid_argument("Invalid value for " + key);
    return value;
}

inline std::vector<double> parse_cfl_list(const std::string& text,
                                           const std::string& key, double max_cfl) {
    std::istringstream stream(text);
    std::vector<double> values;
    std::set<std::string> labels;
    double cfl;
    while (stream >> cfl) {
        if (!std::isfinite(cfl) || cfl <= 0.0 || cfl > max_cfl)
            throw std::invalid_argument("CFL values for " + key + " must satisfy 0 < CFL <= "
                                        + std::to_string(max_cfl));
        if (!labels.insert(cfl_label(cfl)).second)
            throw std::invalid_argument("Duplicate CFL filename for " + key);
        values.push_back(cfl);
    }
    if (!stream.eof() || values.empty())
        throw std::invalid_argument("Invalid or empty CFL list for " + key);
    return values;
}

inline Config read_config(const std::filesystem::path& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Could not open " + filename.string());

    Config config;
    std::set<std::string> keys;
    std::string line;
    while (std::getline(file, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        const auto equal = line.find('=');
        if (equal == std::string::npos) {
            if (line.find_first_not_of(" \t\r") != std::string::npos)
                throw std::invalid_argument("Expected key = value in " + filename.string());
            continue;
        }
        std::istringstream key_stream(line.substr(0, equal));
        std::string key;
        key_stream >> key;
        if (key.empty() || !keys.insert(key).second)
            throw std::invalid_argument("Missing or duplicate key in " + filename.string());
        const std::string value = line.substr(equal + 1);
        if (key == "n") config.n = parse_scalar<int>(value, key);
        else if (key == "max_iterations") config.max_iterations = parse_scalar<int>(value, key);
        else if (key == "parallel_workers") config.parallel_workers = parse_scalar<int>(value, key);
        else if (key == "Mach_in") config.mach_in = parse_scalar<double>(value, key);
        else if (key == "back_pressure_ratio")
            config.back_pressure_ratio = parse_scalar<double>(value, key);
        else if (key == "tolerance_MacCormack")
            config.tolerance_maccormack = parse_scalar<double>(value, key);
        else if (key == "tolerance_Implicit")
            config.tolerance_implicit = parse_scalar<double>(value, key);
        else if (key == "CFL_MacCormack")
            config.cfl_maccormack = parse_cfl_list(value, key, 1.25);
        else if (key == "CFL_Implicit")
            config.cfl_implicit = parse_cfl_list(value, key,
                                                  std::numeric_limits<double>::infinity());
        else throw std::invalid_argument("Unknown Euler input key: " + key);
    }
    if (keys.size() != 9 || config.cfl_maccormack.empty() || config.cfl_implicit.empty())
        throw std::invalid_argument("Euler input requires n, max_iterations, Mach_in, "
                                    "parallel_workers, back_pressure_ratio, two tolerances "
                                    "and two CFL lists");
    if (config.n < 5 || config.max_iterations < 1 || config.parallel_workers < 0
            || !std::isfinite(config.mach_in)
            || config.mach_in <= 1.0 || !std::isfinite(config.back_pressure_ratio)
            || config.back_pressure_ratio <= 0.0
            || !std::isfinite(config.tolerance_maccormack)
            || config.tolerance_maccormack <= 0.0
            || !std::isfinite(config.tolerance_implicit)
            || config.tolerance_implicit <= 0.0)
        throw std::invalid_argument("Invalid Euler mesh, inlet, pressure or tolerance");
    return config;
}

// Run independent cases concurrently. Each case retains its own sequential
// arithmetic, so parallel scheduling does not change numerical results.
inline std::size_t effective_worker_count(std::size_t count, int requested_workers) {
    if (count == 0) return 0;
    const unsigned int available = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t desired = requested_workers == 0
        ? static_cast<std::size_t>(available)
        : static_cast<std::size_t>(requested_workers);
    return std::min(count, std::max<std::size_t>(1, desired));
}

template <typename Function>
inline void parallel_for(std::size_t count, int requested_workers, Function function) {
    const std::size_t worker_count = effective_worker_count(count, requested_workers);
    if (worker_count == 0) return;
    if (worker_count == 1) {
        for (std::size_t index = 0; index < count; ++index) function(index);
        return;
    }

    std::atomic<std::size_t> next{0};
    std::atomic<bool> stop{false};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    const auto worker = [&]() {
        while (!stop.load(std::memory_order_relaxed)) {
            const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
            if (index >= count) return;
            try {
                function(index);
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!first_error) first_error = std::current_exception();
                }
                stop.store(true, std::memory_order_relaxed);
                return;
            }
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (std::size_t i = 0; i < worker_count; ++i) workers.emplace_back(worker);
    for (auto& thread : workers) thread.join();
    if (first_error) std::rethrow_exception(first_error);
}

struct RunOptions {
    std::vector<double> cfl_values;
    int parallel_workers = 0;
};

inline RunOptions parse_run_options(int argc, char** argv,
                                    const std::vector<double>& configured_cfls,
                                    double max_cfl, int configured_workers) {
    RunOptions options{configured_cfls, configured_workers};
    bool cfl_override = false;
    bool workers_override = false;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--cfl" && !cfl_override && i + 1 < argc) {
            const double cfl = parse_scalar<double>(argv[++i], "--cfl");
            if (!std::isfinite(cfl) || cfl <= 0.0 || cfl > max_cfl)
                throw std::invalid_argument("Invalid Euler CFL override");
            options.cfl_values = {cfl};
            cfl_override = true;
        } else if (argument == "--threads" && !workers_override && i + 1 < argc) {
            options.parallel_workers = parse_scalar<int>(argv[++i], "--threads");
            if (options.parallel_workers < 1)
                throw std::invalid_argument("--threads must be a positive integer");
            workers_override = true;
        } else {
            throw std::invalid_argument(
                "Usage: executable [--cfl value] [--threads positive_integer]");
        }
    }
    return options;
}

inline void write_residuals(const std::filesystem::path& filename,
                            const std::vector<double>& residuals,
                            const char* metric) {
    std::ofstream file(filename);
    if (!file) throw std::runtime_error("Could not open " + filename.string());
    file << "# iteration " << metric << '\n' << std::scientific << std::setprecision(16);
    for (std::size_t i = 0; i < residuals.size(); ++i)
        file << i + 1 << ' ' << residuals[i] << '\n';
    file.close();
    if (!file) throw std::runtime_error("Could not write " + filename.string());
}

} // namespace euler_study
