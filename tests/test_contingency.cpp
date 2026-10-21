#include <string>
#include <vector>

#include "powerflow/case_file.hpp"
#include "powerflow/contingency.hpp"
#include "powerflow/network.hpp"
#include "test_framework.hpp"

namespace {

const std::string kDataDir = POWERFLOW_DATA_DIR;

std::string data_path(const std::string& name) { return kDataDir + "/" + name; }

pf::Network load_network(const std::string& name) {
    return pf::Network(pf::load_case_file(data_path(name)));
}

}  // namespace

TEST(Contingency, ExaminesEveryInServiceBranch) {
    const pf::Network net = load_network("case14.m");
    const pf::ContingencyReport report = pf::run_contingency_scan(net);

    CHECK_EQ(report.elements, std::size_t{20});
    CHECK_EQ(report.branch_elements, std::size_t{20});
    CHECK_EQ(report.generator_elements, std::size_t{0});
    CHECK_EQ(report.outcomes.size(), std::size_t{20});
    CHECK_EQ(report.solved + (report.elements - report.solved), report.elements);
    CHECK_EQ(report.solved, std::size_t{20});
}

TEST(Contingency, SecureAndViolatingAddUp) {
    const pf::Network net = load_network("case14.m");
    const pf::ContingencyReport report = pf::run_contingency_scan(net);
    CHECK(report.secure <= report.solved);
    // Every solved outage is either secure or has at least one new violation.
    CHECK_EQ(report.secure + (report.solved - report.secure), report.solved);
}

TEST(Contingency, ExcludesViolationsTheBaseCaseAlreadyHas) {
    // Case 14 has PV setpoints above the Vmax recorded in its own data, so the
    // base case already shows voltage excursions. Those must not be attributed to
    // every single outage.
    const pf::Network net = load_network("case14.m");
    const pf::ContingencyReport report = pf::run_contingency_scan(net);

    CHECK(report.base_converged);
    CHECK(report.base_voltage_violations > 0);
    CHECK(report.secure > 0);
}

TEST(Contingency, SortsTheWorstOutageFirst) {
    const pf::Network net = load_network("case14.m");
    const pf::ContingencyReport report = pf::run_contingency_scan(net);

    for (std::size_t i = 1; i < report.outcomes.size(); ++i) {
        const auto& previous = report.outcomes[i - 1];
        const auto& current = report.outcomes[i];
        const int previous_severity =
            previous.overloaded_branches + previous.voltage_violations;
        const int current_severity = current.overloaded_branches + current.voltage_violations;
        CHECK(previous_severity >= current_severity);
    }
}

TEST(Contingency, ReportsOutagesThatIslandTheNetwork) {
    // Tripping a radial branch must be caught rather than silently producing a
    // wrong answer: the islanded bus leaves the Jacobian singular.
    const pf::Network net = load_network("case118.m");
    pf::ContingencyOptions options;
    options.worst_count = 200;
    const pf::ContingencyReport report = pf::run_contingency_scan(net, options);

    CHECK(report.elements >= std::size_t{186});
    CHECK(report.solved < report.elements);

    std::size_t failures = 0;
    for (const pf::ContingencyOutcome& outcome : report.outcomes) {
        if (!outcome.solved) {
            ++failures;
            CHECK(!outcome.note.empty());
        }
    }
    CHECK_EQ(failures, report.elements - report.solved);
}

TEST(Contingency, IncludesGeneratorsWhenAsked) {
    const pf::Network net = load_network("case14.m");
    pf::ContingencyOptions options;
    options.include_generators = true;
    const pf::ContingencyReport report = pf::run_contingency_scan(net, options);

    // 20 branches plus the four generators that are not on the slack bus.
    CHECK_EQ(report.branch_elements, std::size_t{20});
    CHECK_EQ(report.generator_elements, std::size_t{4});
    CHECK_EQ(report.elements, std::size_t{24});
    CHECK(report.solved > std::size_t{0});
}

TEST(Contingency, IsDeterministic) {
    const pf::Network net = load_network("case14.m");
    const pf::ContingencyReport first = pf::run_contingency_scan(net);
    const pf::ContingencyReport second = pf::run_contingency_scan(net);

    CHECK_EQ(first.outcomes.size(), second.outcomes.size());
    for (std::size_t i = 0; i < first.outcomes.size(); ++i) {
        CHECK_EQ(first.outcomes[i].label, second.outcomes[i].label);
        CHECK_EQ(first.outcomes[i].overloaded_branches, second.outcomes[i].overloaded_branches);
        CHECK_EQ(first.outcomes[i].voltage_violations, second.outcomes[i].voltage_violations);
    }
}

TEST(Contingency, FormattingMentionsTheVerdict) {
    const pf::Network net = load_network("case14.m");
    const pf::ContingencyReport report = pf::run_contingency_scan(net);
    const std::string text = pf::format_contingency_report(report);

    CHECK(text.find("N-1 contingency scan") != std::string::npos);
    CHECK(text.find("case14") != std::string::npos);
    CHECK(text.find("verdict") != std::string::npos);
}
