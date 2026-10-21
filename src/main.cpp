// Command-line entry point: parses arguments, runs the requested solver and
// hands the result to the report / diagram writers.
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "powerflow/case_file.hpp"
#include "powerflow/contingency.hpp"
#include "powerflow/dc_opf.hpp"
#include "powerflow/network.hpp"
#include "powerflow/power_flow.hpp"
#include "powerflow/report.hpp"
#include "powerflow/svg_renderer.hpp"
#include "powerflow/synthetic.hpp"

namespace {

constexpr const char* kVersion = "1.2.0";

struct CliOptions {
    std::string case_path;
    std::string csv_directory;
    std::string json_path;
    std::string svg_path;
    std::string method = "nr";
    std::string linear = "sparse";
    double tolerance = 1e-8;
    int max_iterations = -1;  // -1 picks the per-method default
    int benchmark_buses = 0;
    int contingency_worst = 10;
    bool contingency = false;
    bool contingency_generators = false;
    bool opf = false;
    bool enforce_q_limits = false;
    bool show_tables = true;
    bool show_help = false;
    bool show_version = false;
};

void print_usage(std::ostream& out) {
    out << "powerflow " << kVersion << " - AC power-flow solver\n"
        << "\n"
        << "usage:\n"
        << "  powerflow <case-file> [options]\n"
        << "  powerflow --benchmark <buses>     scaling benchmark on a synthetic grid\n"
        << "\n"
        << "arguments:\n"
        << "  <case-file>              MATPOWER .m case file (mpc.bus / mpc.gen / mpc.branch)\n"
        << "\n"
        << "options:\n"
        << "  -h, --help               show this message\n"
        << "  --version                show the version\n"
        << "  --method <nr|gs|both>    solver selection, default nr\n"
        << "  --linear <sparse|dense>  Jacobian factorisation, default sparse\n"
        << "  --benchmark <buses>      solve a synthetic grid with both factorisations\n"
        << "  --tol <value>            power-mismatch tolerance in pu, default 1e-8\n"
        << "  --max-iter <count>       iteration limit, default 30 (nr) or 2000 (gs)\n"
        << "  --enforce-q-limits       turn PV buses that hit Q limits into PQ buses\n"
        << "  --csv <directory>        write bus.csv, branch.csv, gen.csv, summary.csv\n"
        << "  --json <file>            write the full result as JSON\n"
        << "  --svg <file>             write a single-line diagram as SVG\n"
        << "  --no-tables              print the summary only, skip the tables\n"
        << "\n"
        << "exit codes: 0 converged, 1 usage or data error, 2 not converged\n";
}

void write_text_file(const std::string& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("cannot open '" + path + "' for writing");
    }
    out << content;
    if (!out) {
        throw std::runtime_error("failed while writing '" + path + "'");
    }
}

CliOptions parse_arguments(int argc, char** argv) {
    CliOptions options;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        const auto value_of = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("option ") + name + " needs a value");
            }
            return argv[++i];
        };

        if (argument == "-h" || argument == "--help") {
            options.show_help = true;
        } else if (argument == "--version") {
            options.show_version = true;
        } else if (argument == "--method") {
            options.method = value_of("--method");
            if (options.method != "nr" && options.method != "gs" && options.method != "both") {
                throw std::runtime_error("--method must be one of nr, gs, both");
            }
        } else if (argument == "--linear") {
            options.linear = value_of("--linear");
            if (options.linear != "sparse" && options.linear != "dense") {
                throw std::runtime_error("--linear must be sparse or dense");
            }
        } else if (argument == "--n1") {
            options.contingency = true;
        } else if (argument == "--opf") {
            options.opf = true;
        } else if (argument == "--n1-generators") {
            options.contingency = true;
            options.contingency_generators = true;
        } else if (argument == "--n1-worst") {
            options.contingency_worst = std::stoi(value_of("--n1-worst"));
            if (options.contingency_worst <= 0) {
                throw std::runtime_error("--n1-worst must be positive");
            }
        } else if (argument == "--benchmark") {
            options.benchmark_buses = std::stoi(value_of("--benchmark"));
            if (options.benchmark_buses < 4) {
                throw std::runtime_error("--benchmark needs at least 4 buses");
            }
        } else if (argument == "--tol") {
            options.tolerance = std::stod(value_of("--tol"));
            if (!(options.tolerance > 0.0)) {
                throw std::runtime_error("--tol must be positive");
            }
        } else if (argument == "--max-iter") {
            options.max_iterations = std::stoi(value_of("--max-iter"));
            if (options.max_iterations <= 0) {
                throw std::runtime_error("--max-iter must be positive");
            }
        } else if (argument == "--enforce-q-limits") {
            options.enforce_q_limits = true;
        } else if (argument == "--csv") {
            options.csv_directory = value_of("--csv");
        } else if (argument == "--json") {
            options.json_path = value_of("--json");
        } else if (argument == "--svg") {
            options.svg_path = value_of("--svg");
        } else if (argument == "--no-tables") {
            options.show_tables = false;
        } else if (!argument.empty() && argument[0] == '-' && argument != "-") {
            throw std::runtime_error("unknown option '" + argument + "'");
        } else {
            positional.push_back(argument);
        }
    }

    if (!positional.empty()) {
        options.case_path = positional.front();
    }
    if (positional.size() > 1) {
        throw std::runtime_error("only one case file can be solved at a time");
    }
    return options;
}

std::string format_comparison(const pf::PowerFlowResult& newton, const pf::PowerFlowResult& gauss) {
    std::ostringstream out;
    out << "method comparison\n";
    out << std::left << std::setw(16) << "solver" << std::right << std::setw(12) << "iterations"
        << std::setw(16) << "final mismatch" << '\n';
    out << std::string(44, '-') << '\n';
    out << std::left << std::setw(16) << pf::to_string(newton.method) << std::right
        << std::setw(12) << newton.status.iterations << std::setw(16)
        << newton.status.final_mismatch << '\n';
    out << std::left << std::setw(16) << pf::to_string(gauss.method) << std::right << std::setw(12)
        << gauss.status.iterations << std::setw(16) << gauss.status.final_mismatch << '\n';
    return out.str();
}

std::string format_seconds(double seconds) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << seconds << " s";
    return out.str();
}

// Solves a synthetic grid of the requested size twice: once with the sparse
// factorisation and, when the problem is still small enough to be runnable, once
// with the dense one. The point is to show both that they agree and how the cost
// scales.
int run_benchmark(int bus_count) {
    const pf::CaseData data = pf::make_synthetic_case(bus_count);
    const pf::Network network(data);

    pf::SolverOptions options;
    options.max_iterations = 40;

    const std::size_t n = network.non_slack_buses().size() + network.pq_buses().size();
    std::cout << "synthetic grid : " << bus_count << " buses, " << network.branch_count()
              << " branches\n";

    options.linear_solver = pf::LinearSolver::Sparse;
    const auto sparse_start = std::chrono::steady_clock::now();
    const pf::PowerFlowResult sparse = pf::solve_newton_raphson(network, options);
    const double sparse_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sparse_start).count();

    if (!sparse.status.converged) {
        std::cout << "sparse LU      : NOT CONVERGED (" << sparse.status.message << ")\n";
        return 2;
    }

    const double dense_entries = static_cast<double>(n) * static_cast<double>(n);
    std::cout << "Jacobian       : " << n << " x " << n << ", " << sparse.jacobian_nonzeros
              << " nonzeros ("
              << std::fixed << std::setprecision(2)
              << 100.0 * static_cast<double>(sparse.jacobian_nonzeros) / dense_entries
              << " % of dense)\n";
    std::cout << "sparse LU      : " << sparse.status.iterations << " iterations, "
              << format_seconds(sparse_seconds) << ", factors " << sparse.factor_nonzeros
              << " nnz (includes the one-time fill-reducing ordering)\n";

    // A dense factorisation needs n^2 doubles and n^3/3 flops; past a few hundred
    // unknowns that is no longer a useful comparison, it is just slow.
    constexpr std::size_t kDenseLimit = 900;
    if (n <= kDenseLimit) {
        options.linear_solver = pf::LinearSolver::Dense;
        const auto dense_start = std::chrono::steady_clock::now();
        const pf::PowerFlowResult dense = pf::solve_newton_raphson(network, options);
        const double dense_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - dense_start).count();

        double max_voltage_difference = 0.0;
        for (std::size_t i = 0; i < network.bus_count(); ++i) {
            max_voltage_difference =
                std::max(max_voltage_difference, std::abs(sparse.vm[i] - dense.vm[i]));
        }

        std::cout << "dense LU       : " << dense.status.iterations << " iterations, "
                  << format_seconds(dense_seconds) << "\n";
        std::cout << "agreement      : max |dVm| = " << std::scientific
                  << std::setprecision(2) << max_voltage_difference << " pu, losses differ by "
                  << std::abs(sparse.p_loss_mw - dense.p_loss_mw) << " MW\n";
        if (dense_seconds > 0.0) {
            std::cout << "speed-up       : " << std::fixed << std::setprecision(2)
                      << dense_seconds / sparse_seconds << "x\n";
        }
    } else {
        const double megabytes = dense_entries * 8.0 / (1024.0 * 1024.0);
        std::cout << "dense LU       : skipped (n = " << n << " needs " << std::fixed
                  << std::setprecision(1) << megabytes
                  << " MB for the matrix alone, and O(n^3) work)\n";
    }
    return 0;
}

int run(const CliOptions& options) {
    const pf::CaseData case_data = pf::load_case_file(options.case_path);
    const pf::Network network(case_data);

    pf::SolverOptions solver_options;
    solver_options.tolerance = options.tolerance;
    solver_options.enforce_q_limits = options.enforce_q_limits;
    solver_options.linear_solver =
        options.linear == "dense" ? pf::LinearSolver::Dense : pf::LinearSolver::Sparse;

    if (options.contingency) {
        pf::ContingencyOptions contingency_options;
        contingency_options.solver = solver_options;
        contingency_options.include_generators = options.contingency_generators;
        contingency_options.worst_count = static_cast<std::size_t>(options.contingency_worst);

        const pf::ContingencyReport report =
            pf::run_contingency_scan(network, contingency_options);
        std::cout << pf::format_contingency_report(report, contingency_options);

        // 0: every outage is survivable, 3: an outage breaks a limit,
        // 2: an outage could not even be solved.
        if (report.solved < report.elements) {
            return 2;
        }
        return report.secure == report.elements ? 0 : 3;
    }

    if (options.opf) {
        const pf::DcOpfResult money = pf::solve_dc_opf(network);
        std::cout << pf::format_dc_opf_report(network, money);

        // The network-free dispatch is a lower bound on the cost: any flow limit
        // can only make things more expensive.
        const pf::DcOpfResult unconstrained = pf::solve_economic_dispatch(network);
        std::cout << "\nreference\n";
        std::cout << "  economic dispatch (no network) : " << std::fixed << std::setprecision(2)
                  << unconstrained.objective << " $/h\n";
        std::cout << "  cost of the network limits     : "
                  << (money.objective - unconstrained.objective) << " $/h\n";
        return money.converged ? 0 : 2;
    }

    int exit_code = 0;

    if (options.method == "gs" || options.method == "both") {
        solver_options.max_iterations = options.max_iterations > 0 ? options.max_iterations : 2000;
        const pf::PowerFlowResult gauss = pf::solve_gauss_seidel(network, solver_options);
        std::cout << pf::format_summary(network, gauss) << '\n';
        std::cout << pf::format_convergence(gauss) << '\n';
        if (options.show_tables) {
            std::cout << "bus voltages and injections\n"
                      << pf::format_bus_table(network, gauss) << '\n';
            std::cout << "generation by bus\n"
                      << pf::format_generator_table(network, gauss) << '\n';
            std::cout << "branch flows\n"
                      << pf::format_branch_table(network, gauss) << '\n';
        }
        if (!gauss.status.converged) {
            exit_code = 2;
        }

        if (options.method == "both") {
            solver_options.max_iterations = options.max_iterations > 0 ? options.max_iterations : 30;
            const pf::PowerFlowResult newton = pf::solve_newton_raphson(network, solver_options);
            std::cout << format_comparison(newton, gauss);
        }
        return exit_code;
    }

    solver_options.max_iterations = options.max_iterations > 0 ? options.max_iterations : 30;
    const pf::PowerFlowResult result = pf::solve_newton_raphson(network, solver_options);

    std::cout << pf::format_summary(network, result) << '\n';
    std::cout << pf::format_convergence(result) << '\n';

    if (options.show_tables) {
        std::cout << "bus voltages and injections\n" << pf::format_bus_table(network, result) << '\n';
        std::cout << "generation by bus\n"
                  << pf::format_generator_table(network, result) << '\n';
        std::cout << "branch flows\n" << pf::format_branch_table(network, result) << '\n';
    }

    if (!options.csv_directory.empty()) {
        pf::write_csv_reports(options.csv_directory, network, result);
        std::cout << "wrote bus.csv, branch.csv, gen.csv and summary.csv to "
                  << options.csv_directory << '\n';
    }
    if (!options.json_path.empty()) {
        pf::write_json_report(options.json_path, network, result);
        std::cout << "wrote " << options.json_path << '\n';
    }
    if (!options.svg_path.empty()) {
        write_text_file(options.svg_path, pf::render_svg(network, result));
        std::cout << "wrote " << options.svg_path << '\n';
    }

    if (!result.status.converged) {
        exit_code = 2;
    }
    return exit_code;
}

}  // namespace

int main(int argc, char** argv) {
    CliOptions options;
    try {
        options = parse_arguments(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n\n";
        print_usage(std::cerr);
        return 1;
    }

    if (options.show_help) {
        print_usage(std::cout);
        return 0;
    }
    if (options.show_version) {
        std::cout << "powerflow " << kVersion << '\n';
        return 0;
    }
    if (options.benchmark_buses > 0) {
        try {
            return run_benchmark(options.benchmark_buses);
        } catch (const std::exception& error) {
            std::cerr << "error: " << error.what() << '\n';
            return 1;
        }
    }
    if (options.case_path.empty()) {
        print_usage(std::cerr);
        return 1;
    }

    try {
        return run(options);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
