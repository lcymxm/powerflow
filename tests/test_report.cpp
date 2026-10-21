#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "powerflow/case_file.hpp"
#include "powerflow/network.hpp"
#include "powerflow/power_flow.hpp"
#include "powerflow/report.hpp"
#include "powerflow/svg_renderer.hpp"
#include "test_framework.hpp"

namespace {

const std::string kDataDir = POWERFLOW_DATA_DIR;

std::string data_path(const std::string& name) { return kDataDir + "/" + name; }

std::size_t count_occurrences(const std::string& haystack, const std::string& needle) {
    std::size_t count = 0;
    std::size_t position = 0;
    while ((position = haystack.find(needle, position)) != std::string::npos) {
        ++count;
        position += needle.size();
    }
    return count;
}

std::size_t count_lines(const std::string& text) {
    std::size_t lines = 0;
    for (const char c : text) {
        if (c == '\n') {
            ++lines;
        }
    }
    return lines;
}

std::string read_file(const std::filesystem::path& path) {
    std::ifstream stream(path, std::ios::binary);
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}

// Kept next to the test binary rather than in the system temp directory: on
// Windows the temp directory can be locked down, and CTest sets the working
// directory to the build tree, which is always writable.
std::filesystem::path scratch_directory() {
    const std::filesystem::path directory =
        std::filesystem::current_path() / "powerflow_test_artifacts";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    return directory;
}

struct SolvedCase {
    explicit SolvedCase(const std::string& name)
        : network(pf::load_case_file(data_path(name))), result(pf::solve_newton_raphson(network)) {}

    pf::Network network;
    pf::PowerFlowResult result;
};

}  // namespace

TEST(Report, SummaryMentionsTheCaseAndTheLosses) {
    const SolvedCase solved("case14.m");
    const std::string text = pf::format_summary(solved.network, solved.result);

    CHECK(text.find("case14") != std::string::npos);
    CHECK(text.find("Newton-Raphson") != std::string::npos);
    CHECK(text.find("converged") != std::string::npos);
    CHECK(text.find("losses") != std::string::npos);
    CHECK(text.find("13.393") != std::string::npos);
}

TEST(Report, BusTableHasOneRowPerBus) {
    const SolvedCase solved("case14.m");
    const std::string table = pf::format_bus_table(solved.network, solved.result);

    // one header row, one rule, then one row per bus
    CHECK_EQ(count_lines(table), static_cast<std::size_t>(14 + 2));
    CHECK(table.find("Vm (pu)") != std::string::npos);
}

TEST(Report, BranchTableHasOneRowPerBranch) {
    const SolvedCase solved("case14.m");
    const std::string table = pf::format_branch_table(solved.network, solved.result);
    CHECK_EQ(count_lines(table), static_cast<std::size_t>(20 + 2));
    CHECK(table.find("Loading (%)") != std::string::npos);
}

TEST(Report, GeneratorTableListsGeneratingBusesOnly) {
    const SolvedCase solved("case14.m");
    const std::string table = pf::format_generator_table(solved.network, solved.result);
    CHECK_EQ(count_lines(table), static_cast<std::size_t>(5 + 2));
}

TEST(Report, ConvergenceTableShowsEveryIteration) {
    const SolvedCase solved("case14.m");
    const std::string text = pf::format_convergence(solved.result);
    CHECK_EQ(count_lines(text),
             static_cast<std::size_t>(solved.result.mismatch_history.size() + 1 + 1 + 1));
}

TEST(Report, CsvExportWritesFourFiles) {
    const SolvedCase solved("case14.m");
    const std::filesystem::path directory = scratch_directory() / "csv";
    pf::write_csv_reports(directory.string(), solved.network, solved.result);

    for (const char* name : {"bus.csv", "branch.csv", "gen.csv", "summary.csv"}) {
        CHECK(std::filesystem::exists(directory / name));
    }

    const std::string bus_csv = read_file(directory / "bus.csv");
    CHECK_EQ(count_lines(bus_csv), static_cast<std::size_t>(14 + 1));
    CHECK(bus_csv.find("vm_pu") != std::string::npos);

    const std::string branch_csv = read_file(directory / "branch.csv");
    CHECK_EQ(count_lines(branch_csv), static_cast<std::size_t>(20 + 1));
}

TEST(Report, JsonExportLooksLikeJson) {
    const SolvedCase solved("case30.m");
    const std::filesystem::path path = scratch_directory() / "case30.json";
    pf::write_json_report(path.string(), solved.network, solved.result);

    const std::string text = read_file(path);
    CHECK(!text.empty());
    CHECK_EQ(text.front(), '{');
    CHECK(text.find("\"buses\"") != std::string::npos);
    CHECK(text.find("\"branches\"") != std::string::npos);
    CHECK(text.find("\"iterations\"") != std::string::npos);
    CHECK(text.find("nan") == std::string::npos);
    CHECK(text.find("inf") == std::string::npos);
}

TEST(Report, SvgContainsOneCirclePerBus) {
    const SolvedCase solved("case14.m");
    const std::string svg = pf::render_svg(solved.network, solved.result);

    CHECK(svg.find("<?xml") == 0);
    CHECK(svg.find("<svg") != std::string::npos);
    CHECK(svg.find("</svg>") != std::string::npos);
    CHECK_EQ(count_occurrences(svg, "<circle "), static_cast<std::size_t>(14));
    CHECK_EQ(count_occurrences(svg, "<line "),
             static_cast<std::size_t>(20 + 1 + 3));  // branches, rule, loading samples
    CHECK(svg.find("nan") == std::string::npos);
    CHECK(svg.find("-nan") == std::string::npos);
}

TEST(Report, SvgIsDeterministic) {
    const SolvedCase solved("case14.m");
    const std::string first = pf::render_svg(solved.network, solved.result);
    const std::string second = pf::render_svg(solved.network, solved.result);
    CHECK_EQ(first, second);
}

TEST(Report, SvgRendersLargeCasesWithoutLabels) {
    const SolvedCase solved("case118.m");
    const std::string svg = pf::render_svg(solved.network, solved.result);
    CHECK_EQ(count_occurrences(svg, "<circle "), static_cast<std::size_t>(118));
    // Above 70 buses the labels are suppressed to keep the diagram readable.
    CHECK(svg.find(">1</text>") == std::string::npos);
}
