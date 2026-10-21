#include "powerflow/case_file.hpp"

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pf {
namespace {

// Minimum number of columns we insist on for each matrix. Real MATPOWER files
// carry more (cost data, ramp rates, ...); we simply ignore the extras.
constexpr std::size_t kBusMinCols = 13;
constexpr std::size_t kGenMinCols = 10;
constexpr std::size_t kBranchMinCols = 13;

std::string strip_comments(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    bool in_comment = false;
    for (const char c : text) {
        if (c == '\n') {
            in_comment = false;
            out.push_back(c);
            continue;
        }
        if (in_comment) {
            continue;
        }
        if (c == '%') {
            in_comment = true;
            continue;
        }
        out.push_back(c);
    }
    return out;
}

bool is_identifier_char(char c) {
    const auto uc = static_cast<unsigned char>(c);
    return std::isalnum(uc) != 0 || c == '_';
}

// Finds `key` as a whole token so that `mpc.bus` does not match `mpc.bus_name`
// and `mpc.gen` does not match `mpc.gencost`.
std::size_t find_key(const std::string& text, const std::string& key) {
    std::size_t pos = 0;
    while ((pos = text.find(key, pos)) != std::string::npos) {
        const bool left_ok = pos == 0 || !is_identifier_char(text[pos - 1]);
        const std::size_t after = pos + key.size();
        const bool right_ok = after >= text.size() || !is_identifier_char(text[after]);
        if (left_ok && right_ok) {
            return pos;
        }
        pos += key.size();
    }
    return std::string::npos;
}

void skip_spaces(const std::string& text, std::size_t& i) {
    while (i < text.size() && std::isspace(static_cast<unsigned char>(text[i])) != 0) {
        ++i;
    }
}

double parse_scalar(const std::string& text, std::size_t& i, const std::string& key,
                    const std::string& name) {
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(text.c_str() + i, &end);
    if (end == text.c_str() + i || errno == ERANGE || !std::isfinite(value)) {
        throw std::runtime_error(name + ": could not read a number for " + key);
    }
    i = static_cast<std::size_t>(end - text.c_str());
    return value;
}

std::string matrix_body(const std::string& text, const std::string& key, const std::string& name) {
    const std::size_t pos = find_key(text, key);
    if (pos == std::string::npos) {
        throw std::runtime_error(name + ": missing required section " + key);
    }

    std::size_t i = pos + key.size();
    skip_spaces(text, i);
    if (i >= text.size() || text[i] != '=') {
        throw std::runtime_error(name + ": expected '=' after " + key);
    }
    ++i;
    skip_spaces(text, i);
    if (i >= text.size() || text[i] != '[') {
        throw std::runtime_error(name + ": expected '[' after " + key + " =");
    }
    ++i;

    const std::size_t start = i;
    int depth = 1;
    while (i < text.size()) {
        if (text[i] == '[') {
            ++depth;
        } else if (text[i] == ']') {
            --depth;
            if (depth == 0) {
                return text.substr(start, i - start);
            }
        }
        ++i;
    }
    throw std::runtime_error(name + ": unterminated matrix for " + key);
}

std::vector<double> parse_row(const std::string& row, const std::string& key,
                              const std::string& name, std::size_t row_number) {
    std::vector<double> values;
    std::string token;

    const auto flush = [&]() {
        if (token.empty()) {
            return;
        }
        char* end = nullptr;
        errno = 0;
        const double value = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(value)) {
            throw std::runtime_error(name + ": " + key + " row " + std::to_string(row_number) +
                                     " contains a non-numeric value '" + token + "'");
        }
        values.push_back(value);
        token.clear();
    };

    for (const char c : row) {
        if (std::isspace(static_cast<unsigned char>(c)) != 0 || c == ',') {
            flush();
        } else {
            token.push_back(c);
        }
    }
    flush();
    return values;
}

std::vector<std::vector<double>> parse_matrix(const std::string& text, const std::string& key,
                                              const std::string& name, std::size_t min_cols) {
    const std::string body = matrix_body(text, key, name);
    std::vector<std::vector<double>> rows;
    std::string chunk;

    const auto flush = [&]() {
        const std::size_t row_number = rows.size() + 1;
        std::vector<double> values = parse_row(chunk, key, name, row_number);
        chunk.clear();
        if (values.empty()) {
            return;
        }
        if (values.size() < min_cols) {
            throw std::runtime_error(name + ": " + key + " row " + std::to_string(row_number) +
                                     " has " + std::to_string(values.size()) +
                                     " columns, at least " + std::to_string(min_cols) + " required");
        }
        rows.push_back(std::move(values));
    };

    for (const char c : body) {
        if (c == ';') {
            flush();
        } else {
            chunk.push_back(c);
        }
    }
    flush();
    return rows;
}

double parse_base_mva(const std::string& text, const std::string& name) {
    const std::size_t pos = find_key(text, "mpc.baseMVA");
    if (pos == std::string::npos) {
        throw std::runtime_error(name + ": missing required section mpc.baseMVA");
    }
    std::size_t i = pos + std::string("mpc.baseMVA").size();
    skip_spaces(text, i);
    if (i >= text.size() || text[i] != '=') {
        throw std::runtime_error(name + ": expected '=' after mpc.baseMVA");
    }
    ++i;
    skip_spaces(text, i);
    const double base = parse_scalar(text, i, "mpc.baseMVA", name);
    if (base <= 0.0) {
        throw std::runtime_error(name + ": mpc.baseMVA must be positive, got " +
                                 std::to_string(base));
    }
    return base;
}

BusType decode_bus_type(double raw, std::size_t row_number, const std::string& name) {
    const int code = static_cast<int>(std::lround(raw));
    switch (code) {
        case 1:
            return BusType::PQ;
        case 2:
            return BusType::PV;
        case 3:
            return BusType::Slack;
        case 4:
            return BusType::Isolated;
        default:
            throw std::runtime_error(name + ": mpc.bus row " + std::to_string(row_number) +
                                     " has unsupported bus type " + std::to_string(code) +
                                     " (expected 1..4)");
    }
}

}  // namespace

CaseData load_case_file(const std::string& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        throw std::runtime_error("cannot open case file: " + path);
    }
    std::ostringstream buffer;
    buffer << stream.rdbuf();

    std::string name = path;
    const std::size_t slash = name.find_last_of("/\\");
    if (slash != std::string::npos) {
        name = name.substr(slash + 1);
    }
    const std::size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) {
        name = name.substr(0, dot);
    }

    CaseData data = parse_case_text(buffer.str(), name);
    data.source_path = path;
    return data;
}

CaseData parse_case_text(const std::string& text, const std::string& name) {
    const std::string clean = strip_comments(text);

    CaseData data;
    data.name = name;
    data.base_mva = parse_base_mva(clean, name);

    const auto bus_rows = parse_matrix(clean, "mpc.bus", name, kBusMinCols);
    const auto gen_rows = parse_matrix(clean, "mpc.gen", name, kGenMinCols);
    const auto branch_rows = parse_matrix(clean, "mpc.branch", name, kBranchMinCols);

    if (bus_rows.empty()) {
        throw std::runtime_error(name + ": mpc.bus is empty");
    }

    data.buses.reserve(bus_rows.size());
    for (std::size_t r = 0; r < bus_rows.size(); ++r) {
        const auto& row = bus_rows[r];
        Bus bus;
        bus.id = static_cast<int>(std::lround(row[0]));
        bus.type = decode_bus_type(row[1], r + 1, name);
        bus.pd = row[2];
        bus.qd = row[3];
        bus.gs = row[4];
        bus.bs = row[5];
        bus.vm = row[7];
        bus.va = row[8];
        bus.vmax = row[11];
        bus.vmin = row[12];
        data.buses.push_back(bus);
    }

    data.generators.reserve(gen_rows.size());
    for (const auto& row : gen_rows) {
        Generator gen;
        gen.bus = static_cast<int>(std::lround(row[0]));
        gen.pg = row[1];
        gen.qg = row[2];
        gen.qmax = row[3];
        gen.qmin = row[4];
        gen.vg = row[5];
        gen.online = row[7] > 0.5;
        gen.pmax = row[8];
        gen.pmin = row[9];
        data.generators.push_back(gen);
    }

    data.branches.reserve(branch_rows.size());
    for (const auto& row : branch_rows) {
        Branch branch;
        branch.from = static_cast<int>(std::lround(row[0]));
        branch.to = static_cast<int>(std::lround(row[1]));
        branch.r = row[2];
        branch.x = row[3];
        branch.b = row[4];
        branch.rate_a = row[5];
        branch.ratio = row[8];
        branch.shift = row[9];
        branch.in_service = row[10] > 0.5;
        data.branches.push_back(branch);
    }

    // `mpc.gencost` is optional: only an optimal power flow needs it. Polynomial
    // costs (model 2) with up to three coefficients are understood; piecewise
    // linear costs are not, and simply leave the case without costs.
    if (find_key(clean, "mpc.gencost") != std::string::npos) {
        const auto cost_rows = parse_matrix(clean, "mpc.gencost", name, 4);
        bool complete = cost_rows.size() >= data.generators.size();
        for (std::size_t g = 0; g < data.generators.size() && g < cost_rows.size(); ++g) {
            const auto& row = cost_rows[g];
            const int model = static_cast<int>(std::lround(row[0]));
            const int coefficients = static_cast<int>(std::lround(row[3]));
            if (model != 2 || coefficients < 1 || coefficients > 3 ||
                row.size() < static_cast<std::size_t>(4 + coefficients)) {
                complete = false;
                continue;
            }
            Generator& gen = data.generators[g];
            // Highest order first: c2, c1, c0.
            if (coefficients == 3) {
                gen.cost_c2 = row[4];
                gen.cost_c1 = row[5];
                gen.cost_c0 = row[6];
            } else if (coefficients == 2) {
                gen.cost_c1 = row[4];
                gen.cost_c0 = row[5];
            } else {
                gen.cost_c0 = row[4];
            }
        }
        data.has_generation_costs = complete;
    }

    return data;
}

}  // namespace pf
