#include "powerflow/power_flow.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "powerflow/jacobian.hpp"

namespace pf {

const char* to_string(SolverKind kind) {
    switch (kind) {
        case SolverKind::NewtonRaphson:
            return "Newton-Raphson";
        case SolverKind::GaussSeidel:
            return "Gauss-Seidel";
    }
    return "?";
}

namespace {

// Everything the iteration needs that is not part of the immutable Network:
// the current voltage estimate, the bus classification (which changes when a
// PV bus hits a reactive limit) and the specified injections.
struct WorkingState {
    std::vector<double> vm;       // pu
    std::vector<double> va;       // radians
    std::vector<BusType> type;
    std::vector<Complex> s_spec;  // pu
};

WorkingState make_working_state(const Network& net) {
    const std::size_t nb = net.bus_count();
    const double base = net.base_mva();

    WorkingState state;
    state.vm.assign(nb, 1.0);
    state.va.assign(nb, 0.0);
    state.type.assign(nb, BusType::PQ);
    state.s_spec.assign(nb, Complex{});

    for (Index i = 0; i < nb; ++i) {
        const Bus& bus = net.bus(i);
        state.vm[i] = bus.vm > 0.0 ? bus.vm : 1.0;
        state.va[i] = bus.va * kDegToRad;
        state.type[i] = net.bus_type(i);
        state.s_spec[i] = Complex{-bus.pd, -bus.qd} / base;
    }

    std::vector<bool> magnitude_from_generator(nb, false);
    for (const Generator& gen : net.data().generators) {
        if (!gen.online || !net.has_bus(gen.bus)) {
            continue;
        }
        const Index i = net.index_of(gen.bus);
        state.s_spec[i] += Complex{gen.pg, gen.qg} / base;
        const bool holds_voltage =
            state.type[i] == BusType::PV || state.type[i] == BusType::Slack;
        if (holds_voltage && !magnitude_from_generator[i] && gen.vg > 0.0) {
            state.vm[i] = gen.vg;
            magnitude_from_generator[i] = true;
        }
    }
    return state;
}

std::vector<Complex> voltage_vector(const std::vector<double>& vm,
                                    const std::vector<double>& va_rad) {
    std::vector<Complex> v(vm.size());
    for (std::size_t i = 0; i < vm.size(); ++i) {
        v[i] = std::polar(vm[i], va_rad[i]);
    }
    return v;
}

// S = V .* conj(Y V)
std::vector<Complex> calc_injection(const ComplexMatrix& ybus, const std::vector<Complex>& v) {
    const std::vector<Complex> current = matvec(ybus, v);
    std::vector<Complex> s(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        s[i] = v[i] * std::conj(current[i]);
    }
    return s;
}

std::vector<Index> collect_non_slack(const std::vector<BusType>& types) {
    std::vector<Index> out;
    for (Index i = 0; i < types.size(); ++i) {
        if (types[i] == BusType::PV || types[i] == BusType::PQ) {
            out.push_back(i);
        }
    }
    return out;
}

std::vector<Index> collect_pq(const std::vector<BusType>& types) {
    std::vector<Index> out;
    for (Index i = 0; i < types.size(); ++i) {
        if (types[i] == BusType::PQ) {
            out.push_back(i);
        }
    }
    return out;
}

std::vector<double> compute_mismatch(const std::vector<Complex>& s_spec,
                                     const std::vector<Complex>& s_calc,
                                     const std::vector<Index>& non_slack,
                                     const std::vector<Index>& pq) {
    std::vector<double> mismatch(non_slack.size() + pq.size(), 0.0);
    for (std::size_t k = 0; k < non_slack.size(); ++k) {
        mismatch[k] = s_spec[non_slack[k]].real() - s_calc[non_slack[k]].real();
    }
    for (std::size_t k = 0; k < pq.size(); ++k) {
        mismatch[non_slack.size() + k] = s_spec[pq[k]].imag() - s_calc[pq[k]].imag();
    }
    return mismatch;
}

double largest_abs(const std::vector<double>& values) {
    double worst = 0.0;
    for (const double value : values) {
        worst = std::max(worst, std::abs(value));
    }
    return worst;
}

// Runs Newton iterations until convergence or the iteration budget runs out.
SolverStatus run_newton(const Network& net, WorkingState& state, const SolverOptions& options,
                        std::vector<double>& history) {
    SolverStatus status;
    const std::vector<Index> non_slack = collect_non_slack(state.type);
    const std::vector<Index> pq = collect_pq(state.type);
    const std::size_t n_theta = non_slack.size();

    if (non_slack.empty() && pq.empty()) {
        status.converged = true;
        status.message = "nothing to solve (single slack bus)";
        return status;
    }

    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        const std::vector<Complex> v = voltage_vector(state.vm, state.va);
        const std::vector<Complex> s_calc = calc_injection(net.ybus(), v);
        const std::vector<double> mismatch =
            compute_mismatch(state.s_spec, s_calc, non_slack, pq);
        const double norm = largest_abs(mismatch);

        history.push_back(norm);
        status.iterations = iteration;
        status.final_mismatch = norm;

        if (norm < options.tolerance) {
            status.converged = true;
            status.message = "converged";
            return status;
        }

        std::vector<double> correction;
        try {
            const RealMatrix jac =
                build_power_flow_jacobian(net.ybus(), state.vm, state.va, non_slack, pq);
            LuSolver<double> solver(jac);
            correction = solver.solve(mismatch);
        } catch (const std::exception& error) {
            status.converged = false;
            status.message = std::string("linear solve failed: ") + error.what();
            return status;
        }

        for (std::size_t k = 0; k < n_theta; ++k) {
            state.va[non_slack[k]] += correction[k];
        }
        for (std::size_t k = 0; k < pq.size(); ++k) {
            state.vm[pq[k]] += correction[n_theta + k];
        }
    }

    status.converged = false;
    status.message = "reached the iteration limit (" + std::to_string(options.max_iterations) +
                     ") without meeting the tolerance";
    return status;
}

// Sums the reactive limits of the online generators sitting on one bus.
struct ReactiveLimit {
    double qmax = 0.0;
    double qmin = 0.0;
    bool present = false;
};

ReactiveLimit reactive_limit_of(const Network& net, Index bus_index) {
    ReactiveLimit limit;
    for (const Generator& gen : net.data().generators) {
        if (!gen.online || gen.bus != net.bus(bus_index).id) {
            continue;
        }
        limit.qmax += gen.qmax;
        limit.qmin += gen.qmin;
        limit.present = true;
    }
    return limit;
}

// Demotes every PV bus whose solved reactive output leaves its limits down to a
// PQ bus with Q pinned at the violated bound. Returns true when anything moved.
bool apply_reactive_limits(const Network& net, WorkingState& state, int& switches) {
    const double base = net.base_mva();
    const std::vector<Complex> v = voltage_vector(state.vm, state.va);
    const std::vector<Complex> s_calc = calc_injection(net.ybus(), v);

    bool changed = false;
    for (Index i = 0; i < net.bus_count(); ++i) {
        if (state.type[i] != BusType::PV) {
            continue;
        }
        const ReactiveLimit limit = reactive_limit_of(net, i);
        if (!limit.present) {
            continue;
        }

        const double q_gen_mvar = s_calc[i].imag() * base + net.bus(i).qd;
        constexpr double kLimitTolerance = 1e-5;
        double pinned = 0.0;
        bool violated = false;
        if (q_gen_mvar > limit.qmax + kLimitTolerance) {
            pinned = limit.qmax;
            violated = true;
        } else if (q_gen_mvar < limit.qmin - kLimitTolerance) {
            pinned = limit.qmin;
            violated = true;
        }
        if (!violated) {
            continue;
        }

        state.s_spec[i] = Complex{state.s_spec[i].real(), (pinned - net.bus(i).qd) / base};
        state.type[i] = BusType::PQ;
        ++switches;
        changed = true;
    }
    return changed;
}

void finalise_result(const Network& net, const WorkingState& state, PowerFlowResult& result) {
    const double base = net.base_mva();
    const std::size_t nb = net.bus_count();
    const std::size_t nl = net.branch_count();

    result.base_mva = base;
    result.connected.assign(nb, true);
    for (const Index i : net.isolated_buses()) {
        result.connected[i] = false;
    }

    result.vm = state.vm;
    result.va.assign(nb, 0.0);
    for (Index i = 0; i < nb; ++i) {
        result.va[i] = state.va[i] * kRadToDeg;
    }

    const std::vector<Complex> v = voltage_vector(state.vm, state.va);
    const std::vector<Complex> s_calc = calc_injection(net.ybus(), v);

    result.s_inj = s_calc;
    result.s_gen.assign(nb, Complex{});
    result.s_load.assign(nb, Complex{});
    for (Index i = 0; i < nb; ++i) {
        result.s_load[i] = Complex{net.bus(i).pd, net.bus(i).qd} / base;
        result.s_gen[i] = s_calc[i] + result.s_load[i];
    }
    for (const Index i : net.isolated_buses()) {
        result.s_inj[i] = Complex{};
        result.s_gen[i] = Complex{};
    }

    result.sf.assign(nl, Complex{});
    result.st.assign(nl, Complex{});
    result.loading.assign(nl, 0.0);

    Complex total_loss{};
    for (Index k = 0; k < nl; ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        const BranchAdmittance& adm = net.branch_admittances()[k];
        const Index f = net.index_of(net.branch(k).from);
        const Index t = net.index_of(net.branch(k).to);

        result.sf[k] = v[f] * std::conj(adm.yff * v[f] + adm.yft * v[t]);
        result.st[k] = v[t] * std::conj(adm.ytf * v[f] + adm.ytt * v[t]);
        total_loss += result.sf[k] + result.st[k];

        const double rating = net.branch(k).rate_a;
        if (rating > 0.0) {
            result.loading[k] =
                std::max(std::abs(result.sf[k]), std::abs(result.st[k])) * base / rating;
        }
    }
    result.p_loss_mw = total_loss.real() * base;
    result.q_loss_mvar = total_loss.imag() * base;
}

}  // namespace

PowerFlowResult solve_newton_raphson(const Network& net, const SolverOptions& options) {
    PowerFlowResult result;
    result.method = SolverKind::NewtonRaphson;

    WorkingState state = make_working_state(net);

    SolverStatus status;
    const int rounds = options.enforce_q_limits ? options.max_q_limit_switches + 1 : 1;
    for (int round = 0; round < rounds; ++round) {
        status = run_newton(net, state, options, result.mismatch_history);
        if (!status.converged || !options.enforce_q_limits) {
            break;
        }
        if (!apply_reactive_limits(net, state, result.q_limit_switches)) {
            break;
        }
    }

    result.status = status;
    finalise_result(net, state, result);
    return result;
}

PowerFlowResult solve_gauss_seidel(const Network& net, const SolverOptions& options) {
    const std::size_t nb = net.bus_count();
    const ComplexMatrix& ybus = net.ybus();

    PowerFlowResult result;
    result.method = SolverKind::GaussSeidel;

    WorkingState state = make_working_state(net);
    std::vector<Complex> v = voltage_vector(state.vm, state.va);

    // Buses whose voltage magnitude is held: PV and slack.
    std::vector<double> held_magnitude(nb, 0.0);
    for (Index i = 0; i < nb; ++i) {
        if (state.type[i] == BusType::PV || state.type[i] == BusType::Slack) {
            held_magnitude[i] = state.vm[i];
        }
    }

    const std::vector<Index> non_slack = collect_non_slack(state.type);

    SolverStatus status;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        double max_delta = 0.0;
        bool aborted = false;

        for (const Index i : non_slack) {
            Complex off_diagonal{};
            for (Index j = 0; j < nb; ++j) {
                if (j != i) {
                    off_diagonal += ybus(i, j) * v[j];
                }
            }

            const Complex diagonal = ybus(i, i);
            if (std::abs(diagonal) < 1e-14) {
                status.converged = false;
                status.message = "bus " + std::to_string(net.bus(i).id) +
                                 " has zero self-admittance, Gauss-Seidel cannot proceed";
                aborted = true;
                break;
            }

            Complex specified = state.s_spec[i];
            const bool is_pv = state.type[i] == BusType::PV;
            if (is_pv) {
                // The reactive output of a PV bus is a result, not an input:
                // take it from the current voltage estimate.
                const Complex current = off_diagonal + diagonal * v[i];
                specified = Complex{state.s_spec[i].real(), (v[i] * std::conj(current)).imag()};
            }

            Complex updated = (std::conj(specified / v[i]) - off_diagonal) / diagonal;
            if (is_pv) {
                const double magnitude = std::abs(updated);
                if (magnitude > 0.0) {
                    updated *= held_magnitude[i] / magnitude;
                }
            }

            max_delta = std::max(max_delta, std::abs(updated - v[i]));
            v[i] = updated;
        }

        result.mismatch_history.push_back(max_delta);
        status.iterations = iteration;
        status.final_mismatch = max_delta;

        if (aborted) {
            break;
        }
        if (max_delta < options.tolerance) {
            status.converged = true;
            status.message = "converged";
            break;
        }
    }

    if (!status.converged && status.message.empty()) {
        status.message = "reached the iteration limit (" +
                         std::to_string(options.max_iterations) +
                         ") without meeting the tolerance";
    }

    for (Index i = 0; i < nb; ++i) {
        state.vm[i] = std::abs(v[i]);
        state.va[i] = std::arg(v[i]);
    }

    result.status = status;
    finalise_result(net, state, result);
    return result;
}

double PowerFlowResult::total_generation_mw() const {
    double sum = 0.0;
    for (const Complex& s : s_gen) {
        sum += s.real();
    }
    return sum * base_mva;
}

double PowerFlowResult::total_generation_mvar() const {
    double sum = 0.0;
    for (const Complex& s : s_gen) {
        sum += s.imag();
    }
    return sum * base_mva;
}

double PowerFlowResult::min_voltage() const {
    double best = 0.0;
    bool first = true;
    for (std::size_t i = 0; i < vm.size(); ++i) {
        if (i < connected.size() && !connected[i]) {
            continue;
        }
        if (first || vm[i] < best) {
            best = vm[i];
            first = false;
        }
    }
    return best;
}

double PowerFlowResult::max_voltage() const {
    double best = 0.0;
    bool first = true;
    for (std::size_t i = 0; i < vm.size(); ++i) {
        if (i < connected.size() && !connected[i]) {
            continue;
        }
        if (first || vm[i] > best) {
            best = vm[i];
            first = false;
        }
    }
    return best;
}

Index PowerFlowResult::worst_loaded_branch() const {
    Index worst = kInvalidIndex;
    double peak = 0.0;
    for (std::size_t k = 0; k < loading.size(); ++k) {
        if (loading[k] > peak) {
            peak = loading[k];
            worst = k;
        }
    }
    return worst;
}

}  // namespace pf
