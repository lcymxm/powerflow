// Text, CSV and JSON output for a solved case.
#pragma once

#include <string>
#include <vector>

#include "powerflow/network.hpp"
#include "powerflow/power_flow.hpp"

namespace pf {

// Human-readable blocks, meant for a terminal.
std::string format_summary(const Network& net, const PowerFlowResult& result);
std::string format_convergence(const PowerFlowResult& result);
std::string format_bus_table(const Network& net, const PowerFlowResult& result);
std::string format_generator_table(const Network& net, const PowerFlowResult& result);
std::string format_branch_table(const Network& net, const PowerFlowResult& result);

// Writes bus.csv, gen.csv, branch.csv and summary.csv into `directory`,
// creating the directory if needed.
void write_csv_reports(const std::string& directory, const Network& net,
                       const PowerFlowResult& result);

// Writes a single JSON document holding the summary and all per-bus and
// per-branch values.
void write_json_report(const std::string& path, const Network& net,
                       const PowerFlowResult& result);

}  // namespace pf
