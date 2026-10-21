// N-1 security analysis.
//
// A network operator has to know what happens when a single element is lost: the
// grid must survive any one branch or generator going away. That is exactly what
// this does -- take every element out of service in turn, re-solve the power
// flow, and rank the outages by how badly they break the operating limits.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "powerflow/network.hpp"
#include "powerflow/power_flow.hpp"

namespace pf {

struct ContingencyOptions {
    SolverOptions solver;
    // Also trip every online generator in turn (the classical N-1 "generator
    // contingency"). The slack bus is never tripped, it is what absorbs the loss.
    bool include_generators = false;
    // How many entries the formatted report lists.
    std::size_t worst_count = 10;
};

struct ContingencyOutcome {
    std::string label;  // e.g. "branch 4-7" or "generator at bus 2"
    bool solved = false;
    std::string note;  // why it could not be solved, when it could not

    int overloaded_branches = 0;
    int voltage_violations = 0;
    int isolated_buses = 0;
    double worst_loading = 0.0;  // pu of rate A
    double min_voltage = 0.0;
    double max_voltage = 0.0;
    double p_loss_mw = 0.0;
};

struct ContingencyReport {
    std::string case_name;
    std::size_t elements = 0;
    std::size_t branch_elements = 0;
    std::size_t generator_elements = 0;
    std::size_t solved = 0;
    std::size_t secure = 0;  // solved and inside every limit
    // Violations the case already has before any outage. They are excluded from
    // the per-contingency counts, otherwise a network whose base case already
    // breaks a limit would report every single outage as a problem.
    std::size_t base_overloaded_branches = 0;
    std::size_t base_voltage_violations = 0;
    bool base_converged = false;
    std::vector<ContingencyOutcome> outcomes;  // worst first
};

// Runs the scan on `base`. The base case itself is not modified.
ContingencyReport run_contingency_scan(const Network& base,
                                       const ContingencyOptions& options = {});

std::string format_contingency_report(const ContingencyReport& report,
                                      const ContingencyOptions& options = {});

}  // namespace pf
