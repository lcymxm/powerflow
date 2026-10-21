#include "powerflow/report.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pf {
namespace {

std::string fixed(double value, int precision) {
    // Snap anything that would print as a signed zero: "-0.000" reads like a real
    // (if tiny) result, while the value is pure numerical noise around zero.
    const double quantum = 0.5 * std::pow(10.0, -precision);
    if (std::abs(value) < quantum) {
        value = 0.0;
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

std::string scientific(double value, int precision = 3) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(precision) << value;
    return out.str();
}

std::string json_number(double value) {
    if (!std::isfinite(value)) {
        return "null";
    }
    std::ostringstream out;
    out << std::setprecision(12) << value;
    return out.str();
}

std::string json_string(const std::string& text) {
    std::string out;
    out.reserve(text.size() + 2);
    out.push_back('"');
    for (const char c : text) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\t':
                out += "\\t";
                break;
            default:
                out.push_back(c);
                break;
        }
    }
    out.push_back('"');
    return out;
}

// Fixed-width columns, right aligned, with a dashed rule under the header.
class Table {
public:
    void header(std::vector<std::string> cells) { rows_.push_back(std::move(cells)); }
    void row(std::vector<std::string> cells) { rows_.push_back(std::move(cells)); }

    std::string render() const {
        if (rows_.empty()) {
            return {};
        }
        std::size_t columns = 0;
        for (const auto& cells : rows_) {
            columns = std::max(columns, cells.size());
        }
        std::vector<std::size_t> width(columns, 0);
        for (const auto& cells : rows_) {
            for (std::size_t c = 0; c < cells.size(); ++c) {
                width[c] = std::max(width[c], cells[c].size());
            }
        }

        std::ostringstream out;
        for (std::size_t r = 0; r < rows_.size(); ++r) {
            for (std::size_t c = 0; c < columns; ++c) {
                const std::string cell = c < rows_[r].size() ? rows_[r][c] : std::string{};
                if (c > 0) {
                    out << "  ";
                }
                out << std::right << std::setw(static_cast<int>(width[c])) << cell;
            }
            out << '\n';
            if (r == 0) {
                for (std::size_t c = 0; c < columns; ++c) {
                    if (c > 0) {
                        out << "  ";
                    }
                    out << std::string(width[c], '-');
                }
                out << '\n';
            }
        }
        return out.str();
    }

private:
    std::vector<std::vector<std::string>> rows_;
};

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

}  // namespace

std::string format_summary(const Network& net, const PowerFlowResult& result) {
    std::ostringstream out;
    const CaseData& data = net.data();

    out << "case          : " << data.name << "  (" << net.bus_count() << " buses, "
        << net.branch_count() << " branches, " << data.generators.size() << " generators)\n";

    out << "solver        : " << to_string(result.method) << '\n';

    out << "status        : ";
    if (result.status.converged) {
        out << "converged in " << result.status.iterations << " iteration"
            << (result.status.iterations == 1 ? "" : "s") << ", max mismatch "
            << scientific(result.status.final_mismatch) << " pu";
    } else {
        out << "NOT CONVERGED (" << result.status.message << ")";
    }
    out << '\n';

    if (result.q_limit_switches > 0) {
        out << "q limits      : " << result.q_limit_switches
            << " PV bus(es) demoted to PQ\n";
    }

    out << "load          : " << fixed(net.total_load_mw(), 3) << " MW  /  "
        << fixed(net.total_load_mvar(), 3) << " MVAr\n";
    out << "generation    : " << fixed(result.total_generation_mw(), 3) << " MW  /  "
        << fixed(result.total_generation_mvar(), 3) << " MVAr\n";
    out << "losses        : " << fixed(result.p_loss_mw, 3) << " MW  /  "
        << fixed(result.q_loss_mvar, 3) << " MVAr\n";
    out << "voltage range : " << fixed(result.min_voltage(), 4) << " .. "
        << fixed(result.max_voltage(), 4) << " pu\n";

    const Index worst = result.worst_loaded_branch();
    if (worst != kInvalidIndex) {
        const Branch& branch = net.branch(worst);
        out << "worst branch  : " << branch.from << "-" << branch.to << " at "
            << fixed(result.loading[worst] * 100.0, 1) << " % of rate A\n";
    }
    return out.str();
}

std::string format_convergence(const PowerFlowResult& result) {
    const bool newton = result.method == SolverKind::NewtonRaphson;
    std::ostringstream out;
    out << (newton ? "convergence: max |dP|, |dQ| in pu" : "convergence: max |dV| in pu") << '\n';

    Table table;
    table.header({"iteration", newton ? "|dP|,|dQ|" : "|dV|", "ratio"});
    double previous = 0.0;
    for (std::size_t i = 0; i < result.mismatch_history.size(); ++i) {
        const double value = result.mismatch_history[i];
        std::string ratio = "-";
        if (i > 0 && previous > 0.0) {
            ratio = scientific(value / previous, 2);
        }
        table.row({std::to_string(i + 1), scientific(value), ratio});
        previous = value;
    }
    out << table.render();
    return out.str();
}

std::string format_bus_table(const Network& net, const PowerFlowResult& result) {
    const double base = net.base_mva();
    Table table;
    table.header({"Bus", "Type", "Vm (pu)", "Va (deg)", "Pgen (MW)", "Qgen (MVAr)",
                  "Pload (MW)", "Qload (MVAr)", "Flag"});

    for (Index i = 0; i < net.bus_count(); ++i) {
        const Bus& bus = net.bus(i);
        const BusType type = net.bus_type(i);
        std::string flag;
        if (!result.connected[i]) {
            flag = "isolated";
        } else if (type == BusType::PQ) {
            // Only PQ buses have a solved voltage; on a PV or slack bus the
            // magnitude is a setpoint, so "outside Vmax" there is a data problem
            // rather than a result worth flagging.
            if (bus.vmax > 0.0 && result.vm[i] > bus.vmax + 1e-9) {
                flag = "over-Vmax";
            } else if (bus.vmin > 0.0 && result.vm[i] < bus.vmin - 1e-9) {
                flag = "under-Vmin";
            }
        }

        table.row({std::to_string(bus.id), to_string(type), fixed(result.vm[i], 4),
                   fixed(result.va[i], 3), fixed(result.s_gen[i].real() * base, 3),
                   fixed(result.s_gen[i].imag() * base, 3), fixed(bus.pd, 3),
                   fixed(bus.qd, 3), flag});
    }
    return table.render();
}

std::string format_generator_table(const Network& net, const PowerFlowResult& result) {
    const double base = net.base_mva();

    Table table;
    table.header({"Bus", "Pg (MW)", "Qg (MVAr)", "Qmax", "Qmin", "Flag"});

    for (Index i = 0; i < net.bus_count(); ++i) {
        double qmax = 0.0;
        double qmin = 0.0;
        std::size_t count = 0;
        for (const Generator& gen : net.data().generators) {
            if (gen.online && gen.bus == net.bus(i).id) {
                qmax += gen.qmax;
                qmin += gen.qmin;
                ++count;
            }
        }
        if (count == 0) {
            continue;
        }

        const double pg = result.s_gen[i].real() * base;
        const double qg = result.s_gen[i].imag() * base;
        std::string flag;
        if (qg >= qmax - 1e-3) {
            flag = "at Qmax";
        } else if (qg <= qmin + 1e-3) {
            flag = "at Qmin";
        }

        table.row({std::to_string(net.bus(i).id), fixed(pg, 3), fixed(qg, 3), fixed(qmax, 3),
                   fixed(qmin, 3), flag});
    }
    return table.render();
}

std::string format_branch_table(const Network& net, const PowerFlowResult& result) {
    const double base = net.base_mva();

    Table table;
    table.header({"From", "To", "P (MW)", "Q (MVAr)", "Ploss (MW)", "Qloss (MVAr)", "Loading (%)",
                  "Status"});

    for (Index k = 0; k < net.branch_count(); ++k) {
        const Branch& branch = net.branch(k);
        const bool active = net.branch_in_service(k);
        const Complex loss = result.sf[k] + result.st[k];

        std::string loading = "-";
        if (active && branch.rate_a > 0.0) {
            loading = fixed(result.loading[k] * 100.0, 1);
        }

        table.row({std::to_string(branch.from), std::to_string(branch.to),
                   fixed(result.sf[k].real() * base, 3), fixed(result.sf[k].imag() * base, 3),
                   fixed(loss.real() * base, 3), fixed(loss.imag() * base, 3), loading,
                   active ? "on" : "out"});
    }
    return table.render();
}

void write_csv_reports(const std::string& directory, const Network& net,
                       const PowerFlowResult& result) {
    std::error_code code;
    std::filesystem::create_directories(directory, code);
    if (code) {
        throw std::runtime_error("could not create output directory '" + directory +
                                 "': " + code.message());
    }

    const double base = net.base_mva();

    {
        std::ostringstream out;
        out << "bus,type,vm_pu,va_deg,pg_mw,qg_mvar,pd_mw,qd_mvar,vmax_pu,vmin_pu,connected\n";
        for (Index i = 0; i < net.bus_count(); ++i) {
            const Bus& bus = net.bus(i);
            out << bus.id << ',' << to_string(net.bus_type(i)) << ',' << fixed(result.vm[i], 6)
                << ',' << fixed(result.va[i], 6) << ',' << fixed(result.s_gen[i].real() * base, 6)
                << ',' << fixed(result.s_gen[i].imag() * base, 6) << ',' << fixed(bus.pd, 6)
                << ',' << fixed(bus.qd, 6) << ',' << fixed(bus.vmax, 6) << ','
                << fixed(bus.vmin, 6) << ',' << (result.connected[i] ? 1 : 0) << '\n';
        }
        write_text_file((std::filesystem::path(directory) / "bus.csv").string(), out.str());
    }

    {
        std::ostringstream out;
        out << "from,to,p_from_mw,q_from_mvar,p_to_mw,q_to_mvar,p_loss_mw,q_loss_mvar,"
               "loading_pct,in_service\n";
        for (Index k = 0; k < net.branch_count(); ++k) {
            const Branch& branch = net.branch(k);
            const Complex loss = result.sf[k] + result.st[k];
            const double loading = branch.rate_a > 0.0 ? result.loading[k] * 100.0 : 0.0;
            out << branch.from << ',' << branch.to << ',' << fixed(result.sf[k].real() * base, 6)
                << ',' << fixed(result.sf[k].imag() * base, 6) << ','
                << fixed(result.st[k].real() * base, 6) << ','
                << fixed(result.st[k].imag() * base, 6) << ',' << fixed(loss.real() * base, 6)
                << ',' << fixed(loss.imag() * base, 6) << ',' << fixed(loading, 4) << ','
                << (net.branch_in_service(k) ? 1 : 0) << '\n';
        }
        write_text_file((std::filesystem::path(directory) / "branch.csv").string(), out.str());
    }

    {
        std::ostringstream out;
        out << "bus,pg_mw,qg_mvar,vm_pu,va_deg\n";
        for (Index i = 0; i < net.bus_count(); ++i) {
            double online = 0.0;
            for (const Generator& gen : net.data().generators) {
                if (gen.online && gen.bus == net.bus(i).id) {
                    online += 1.0;
                }
            }
            if (online == 0.0) {
                continue;
            }
            out << net.bus(i).id << ',' << fixed(result.s_gen[i].real() * base, 6) << ','
                << fixed(result.s_gen[i].imag() * base, 6) << ',' << fixed(result.vm[i], 6) << ','
                << fixed(result.va[i], 6) << '\n';
        }
        write_text_file((std::filesystem::path(directory) / "gen.csv").string(), out.str());
    }

    {
        std::ostringstream out;
        out << "case,buses,branches,generators,solver,converged,iterations,max_mismatch_pu,"
               "load_mw,load_mvar,generation_mw,generation_mvar,p_loss_mw,q_loss_mvar,"
               "vm_min,vm_max,q_limit_switches\n";
        out << net.data().name << ',' << net.bus_count() << ',' << net.branch_count() << ','
            << net.data().generators.size() << ',' << to_string(result.method) << ','
            << (result.status.converged ? 1 : 0) << ',' << result.status.iterations << ','
            << json_number(result.status.final_mismatch) << ','
            << fixed(net.total_load_mw(), 6) << ',' << fixed(net.total_load_mvar(), 6) << ','
            << fixed(result.total_generation_mw(), 6) << ','
            << fixed(result.total_generation_mvar(), 6) << ',' << fixed(result.p_loss_mw, 6)
            << ',' << fixed(result.q_loss_mvar, 6) << ',' << fixed(result.min_voltage(), 6) << ','
            << fixed(result.max_voltage(), 6) << ',' << result.q_limit_switches << '\n';
        write_text_file((std::filesystem::path(directory) / "summary.csv").string(), out.str());
    }
}

void write_json_report(const std::string& path, const Network& net,
                       const PowerFlowResult& result) {
    const double base = net.base_mva();

    std::ostringstream out;
    out << "{\n";
    out << "  \"case\": " << json_string(net.data().name) << ",\n";
    out << "  \"source_path\": " << json_string(net.data().source_path) << ",\n";
    out << "  \"base_mva\": " << json_number(base) << ",\n";
    out << "  \"solver\": " << json_string(to_string(result.method)) << ",\n";
    out << "  \"converged\": " << (result.status.converged ? "true" : "false") << ",\n";
    out << "  \"iterations\": " << result.status.iterations << ",\n";
    out << "  \"max_mismatch_pu\": " << json_number(result.status.final_mismatch) << ",\n";
    out << "  \"message\": " << json_string(result.status.message) << ",\n";
    out << "  \"q_limit_switches\": " << result.q_limit_switches << ",\n";

    out << "  \"summary\": {\n";
    out << "    \"bus_count\": " << net.bus_count() << ",\n";
    out << "    \"branch_count\": " << net.branch_count() << ",\n";
    out << "    \"load_mw\": " << json_number(net.total_load_mw()) << ",\n";
    out << "    \"load_mvar\": " << json_number(net.total_load_mvar()) << ",\n";
    out << "    \"generation_mw\": " << json_number(result.total_generation_mw()) << ",\n";
    out << "    \"generation_mvar\": " << json_number(result.total_generation_mvar()) << ",\n";
    out << "    \"p_loss_mw\": " << json_number(result.p_loss_mw) << ",\n";
    out << "    \"q_loss_mvar\": " << json_number(result.q_loss_mvar) << ",\n";
    out << "    \"vm_min\": " << json_number(result.min_voltage()) << ",\n";
    out << "    \"vm_max\": " << json_number(result.max_voltage()) << "\n";
    out << "  },\n";

    out << "  \"buses\": [\n";
    for (Index i = 0; i < net.bus_count(); ++i) {
        const Bus& bus = net.bus(i);
        out << "    {\"bus\": " << bus.id << ", \"type\": " << json_string(to_string(net.bus_type(i)))
            << ", \"vm\": " << json_number(result.vm[i])
            << ", \"va\": " << json_number(result.va[i])
            << ", \"pg_mw\": " << json_number(result.s_gen[i].real() * base)
            << ", \"qg_mvar\": " << json_number(result.s_gen[i].imag() * base)
            << ", \"pd_mw\": " << json_number(bus.pd)
            << ", \"qd_mvar\": " << json_number(bus.qd) << "}"
            << (i + 1 < net.bus_count() ? "," : "") << "\n";
    }
    out << "  ],\n";

    out << "  \"branches\": [\n";
    for (Index k = 0; k < net.branch_count(); ++k) {
        const Branch& branch = net.branch(k);
        const Complex loss = result.sf[k] + result.st[k];
        out << "    {\"from\": " << branch.from << ", \"to\": " << branch.to
            << ", \"p_from_mw\": " << json_number(result.sf[k].real() * base)
            << ", \"q_from_mvar\": " << json_number(result.sf[k].imag() * base)
            << ", \"p_to_mw\": " << json_number(result.st[k].real() * base)
            << ", \"q_to_mvar\": " << json_number(result.st[k].imag() * base)
            << ", \"p_loss_mw\": " << json_number(loss.real() * base)
            << ", \"q_loss_mvar\": " << json_number(loss.imag() * base)
            << ", \"loading_pct\": " << json_number(result.loading[k] * 100.0)
            << ", \"in_service\": " << (net.branch_in_service(k) ? "true" : "false") << "}"
            << (k + 1 < net.branch_count() ? "," : "") << "\n";
    }
    out << "  ]\n";
    out << "}\n";

    write_text_file(path, out.str());
}

}  // namespace pf
