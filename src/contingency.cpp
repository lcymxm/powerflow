#include "powerflow/contingency.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace pf {
namespace {

// Which limits the un-contingenced case already breaks. A contingency only
// counts the *new* violations, otherwise a case whose data already sits outside
// its own limits would flag every outage.
struct Baseline {
    std::vector<char> overloaded_branch;
    std::vector<char> voltage_violation;

    std::size_t overloaded_count() const {
        return static_cast<std::size_t>(std::count(overloaded_branch.begin(),
                                                   overloaded_branch.end(), '\x01'));
    }
    std::size_t voltage_count() const {
        return static_cast<std::size_t>(std::count(voltage_violation.begin(),
                                                   voltage_violation.end(), '\x01'));
    }
};

Baseline make_baseline(const Network& net, const PowerFlowResult& result) {
    Baseline baseline;
    baseline.overloaded_branch.assign(net.branch_count(), 0);
    baseline.voltage_violation.assign(net.bus_count(), 0);

    for (Index k = 0; k < net.branch_count(); ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        if (net.branch(k).rate_a > 0.0 && result.loading[k] > 1.0) {
            baseline.overloaded_branch[k] = 1;
        }
    }
    for (Index i = 0; i < net.bus_count(); ++i) {
        if (!result.connected[i]) {
            continue;
        }
        const Bus& bus = net.bus(i);
        if ((bus.vmax > 0.0 && result.vm[i] > bus.vmax + 1e-9) ||
            (bus.vmin > 0.0 && result.vm[i] < bus.vmin - 1e-9)) {
            baseline.voltage_violation[i] = 1;
        }
    }
    return baseline;
}

// Copies the case, then demotes any voltage-controlled bus that no longer has an
// online generator to a PQ bus. Without this a generator outage would leave a PV
// bus with an unconstrained reactive injection, which is not what an operator
// would see.
void demote_buses_without_generation(CaseData& data) {
    for (Bus& bus : data.buses) {
        if (bus.type != BusType::PV) {
            continue;
        }
        const bool still_has_generator =
            std::any_of(data.generators.begin(), data.generators.end(),
                        [&bus](const Generator& generator) {
                            return generator.online && generator.bus == bus.id;
                        });
        if (!still_has_generator) {
            bus.type = BusType::PQ;
        }
    }
}

ContingencyOutcome summarise(const Network& net, const PowerFlowResult& result,
                             const Baseline& baseline, const std::string& label) {
    ContingencyOutcome outcome;
    outcome.label = label;
    outcome.solved = result.status.converged;
    if (!outcome.solved) {
        outcome.note = result.status.message;
        return outcome;
    }

    outcome.isolated_buses = static_cast<int>(net.isolated_buses().size());

    for (Index k = 0; k < net.branch_count(); ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        const double loading = result.loading[k];
        if (loading > outcome.worst_loading) {
            outcome.worst_loading = loading;
        }
        const bool overloaded = net.branch(k).rate_a > 0.0 && loading > 1.0;
        const bool already = baseline.overloaded_branch[k] != 0;
        if (overloaded && !already) {
            ++outcome.overloaded_branches;
        }
    }

    for (Index i = 0; i < net.bus_count(); ++i) {
        if (!result.connected[i]) {
            continue;
        }
        const Bus& bus = net.bus(i);
        const bool out_of_limits = (bus.vmax > 0.0 && result.vm[i] > bus.vmax + 1e-9) ||
                                   (bus.vmin > 0.0 && result.vm[i] < bus.vmin - 1e-9);
        const bool already = baseline.voltage_violation[i] != 0;
        if (out_of_limits && !already) {
            ++outcome.voltage_violations;
        }
    }

    outcome.min_voltage = result.min_voltage();
    outcome.max_voltage = result.max_voltage();
    outcome.p_loss_mw = result.p_loss_mw;
    return outcome;
}

int severity(const ContingencyOutcome& outcome) {
    return outcome.overloaded_branches + outcome.voltage_violations;
}

}  // namespace

ContingencyReport run_contingency_scan(const Network& base, const ContingencyOptions& options) {
    ContingencyReport report;
    report.case_name = base.data().name;
    const CaseData& data = base.data();

    const PowerFlowResult base_result = solve_newton_raphson(base, options.solver);
    report.base_converged = base_result.status.converged;
    const Baseline baseline = make_baseline(base, base_result);
    report.base_overloaded_branches = baseline.overloaded_count();
    report.base_voltage_violations = baseline.voltage_count();

    for (Index k = 0; k < data.branches.size(); ++k) {
        if (!data.branches[k].in_service) {
            continue;
        }
        CaseData modified = data;
        modified.branches[k].in_service = false;

        const std::string label = "branch " + std::to_string(data.branches[k].from) + "-" +
                                  std::to_string(data.branches[k].to);
        ++report.elements;
        ++report.branch_elements;

        try {
            const Network tripped(modified);
            const PowerFlowResult result = solve_newton_raphson(tripped, options.solver);
            report.outcomes.push_back(summarise(tripped, result, baseline, label));
        } catch (const std::exception& error) {
            ContingencyOutcome outcome;
            outcome.label = label;
            outcome.note = error.what();
            report.outcomes.push_back(outcome);
        }
    }

    if (options.include_generators) {
        const int slack_bus_id = data.buses[base.slack_bus()].id;
        for (Index g = 0; g < data.generators.size(); ++g) {
            if (!data.generators[g].online || data.generators[g].bus == slack_bus_id) {
                continue;
            }
            CaseData modified = data;
            modified.generators[g].online = false;
            demote_buses_without_generation(modified);

            const std::string label = "generator at bus " + std::to_string(data.generators[g].bus);
            ++report.elements;
            ++report.generator_elements;

            try {
                const Network tripped(modified);
                const PowerFlowResult result = solve_newton_raphson(tripped, options.solver);
                report.outcomes.push_back(summarise(tripped, result, baseline, label));
            } catch (const std::exception& error) {
                ContingencyOutcome outcome;
                outcome.label = label;
                outcome.note = error.what();
                report.outcomes.push_back(outcome);
            }
        }
    }

    std::stable_sort(report.outcomes.begin(), report.outcomes.end(),
                     [](const ContingencyOutcome& a, const ContingencyOutcome& b) {
                         if (a.solved != b.solved) {
                             return !a.solved;  // the ones that did not even solve first
                         }
                         if (severity(a) != severity(b)) {
                             return severity(a) > severity(b);
                         }
                         return a.worst_loading > b.worst_loading;
                     });

    for (const ContingencyOutcome& outcome : report.outcomes) {
        if (outcome.solved) {
            ++report.solved;
            if (severity(outcome) == 0 && outcome.isolated_buses == 0) {
                ++report.secure;
            }
        }
    }
    return report;
}

std::string format_contingency_report(const ContingencyReport& report,
                                      const ContingencyOptions& options) {
    std::ostringstream out;

    out << "N-1 contingency scan : " << report.case_name << ", " << report.elements
        << " elements\n";
    out << "  branch outages    : " << report.branch_elements << '\n';
    if (report.generator_elements > 0) {
        out << "  generator outages : " << report.generator_elements << '\n';
    }
    out << "  solved            : " << report.solved << " / " << report.elements << '\n';
    out << "  secure            : " << report.secure << '\n';
    out << "  with violations   : " << (report.solved - report.secure) << '\n';
    if (report.solved < report.elements) {
        out << "  failed to solve   : " << (report.elements - report.solved)
            << "  (islanded or lost the slack)\n";
    }
    if (report.base_overloaded_branches > 0 || report.base_voltage_violations > 0) {
        out << "  base case already : " << report.base_overloaded_branches
            << " branch overload(s), " << report.base_voltage_violations
            << " voltage excursion(s); these are excluded from the counts below\n";
    }
    out << '\n';

    const std::size_t shown = std::min(options.worst_count, report.outcomes.size());
    if (shown == 0) {
        out << "  no contingencies to report\n";
        return out.str();
    }

    out << "worst " << shown << '\n';
    out << std::left << std::setw(24) << "element" << std::right << std::setw(8) << "solved"
        << std::setw(12) << "overloaded" << std::setw(9) << "voltage" << std::setw(10)
        << "max load" << std::setw(10) << "min Vm" << '\n';
    out << std::string(73, '-') << '\n';

    const auto previous_flags = out.flags();
    for (std::size_t i = 0; i < shown; ++i) {
        const ContingencyOutcome& outcome = report.outcomes[i];
        std::string label = outcome.label;
        if (label.size() > 23) {
            label.resize(23);
        }
        out << std::left << std::setw(24) << label << std::right << std::setw(8)
            << (outcome.solved ? "yes" : "no");
        if (outcome.solved) {
            out << std::setw(12) << outcome.overloaded_branches << std::setw(9)
                << outcome.voltage_violations << std::setw(10) << std::fixed
                << std::setprecision(2) << outcome.worst_loading << std::setw(10)
                << std::setprecision(4) << outcome.min_voltage;
        } else {
            out << std::setw(12) << "-" << std::setw(9) << "-" << std::setw(10) << "-"
                << std::setw(10) << "-";
        }
        out << '\n';
    }
    out.flags(previous_flags);

    if (report.solved < report.elements) {
        out << "\nfailures\n";
        for (const ContingencyOutcome& outcome : report.outcomes) {
            if (outcome.solved) {
                continue;
            }
            out << "  " << std::left << std::setw(24) << outcome.label << outcome.note << '\n';
        }
    }

    const bool all_secure = report.secure == report.elements;
    out << '\n'
        << (all_secure ? "verdict: the system survives every single outage examined\n"
                       : "verdict: at least one single outage breaks an operating limit\n");
    return out.str();
}

}  // namespace pf
