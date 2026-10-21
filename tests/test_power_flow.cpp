#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <utility>
#include <vector>

#include "powerflow/case_file.hpp"
#include "powerflow/jacobian.hpp"
#include "powerflow/network.hpp"
#include "powerflow/power_flow.hpp"
#include "test_framework.hpp"

namespace {

const std::string kDataDir = POWERFLOW_DATA_DIR;

std::string data_path(const std::string& name) { return kDataDir + "/" + name; }

pf::Network load_network(const std::string& name) {
    return pf::Network(pf::load_case_file(data_path(name)));
}

pf::SolverOptions tight_options() {
    pf::SolverOptions options;
    options.tolerance = 1e-12;
    options.max_iterations = 50;
    return options;
}

// Reference active/reactive losses in MW / MVAr, produced by an independent
// numpy implementation of the same equations (see tests/README.md).
struct LossReference {
    const char* name;
    double p_loss_mw;
    double q_loss_mvar;
};

const LossReference kReferences[] = {
    {"case9.m", 4.64102147, -92.16012522},
    {"case14.m", 13.39327236, 30.12238805},
    {"case30.m", 2.44380314, -6.56273064},
    {"case57.m", 27.86375151, 6.32797220},
    {"case118.m", 132.86287189, -557.94742312},
};

}  // namespace

TEST(PowerFlow, Ieee14ConvergesQuickly) {
    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());
    CHECK(result.status.converged);
    CHECK(result.status.iterations <= 6);
    CHECK(result.status.final_mismatch < 1e-12);
    CHECK_EQ(result.mismatch_history.size(), static_cast<std::size_t>(result.status.iterations));
}

TEST(PowerFlow, LossesMatchReferenceForEveryCase) {
    for (const LossReference& reference : kReferences) {
        const pf::Network net = load_network(reference.name);
        const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());
        CHECK(result.status.converged);
        CHECK_NEAR(result.p_loss_mw, reference.p_loss_mw, 1e-5);
        CHECK_NEAR(result.q_loss_mvar, reference.q_loss_mvar, 1e-4);
    }
}

TEST(PowerFlow, Ieee14GenerationCoversLoadPlusLosses) {
    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());
    const double base = net.base_mva();

    // Generation covers the load, the branch losses and whatever the bus shunts
    // absorb. Case 14 has a 19 MVAr capacitor at bus 9, so it covers *less* than
    // load plus losses: the capacitor supplies 21.18 MVAr locally.
    pf::Complex shunt{0.0, 0.0};
    for (std::size_t i = 0; i < net.bus_count(); ++i) {
        shunt += (result.vm[i] * result.vm[i]) * std::conj(net.bus_shunts()[i]);
    }

    CHECK_NEAR(result.total_generation_mw(),
               net.total_load_mw() + result.p_loss_mw + shunt.real() * base, 1e-6);
    CHECK_NEAR(result.total_generation_mvar(),
               net.total_load_mvar() + result.q_loss_mvar + shunt.imag() * base, 1e-6);

    CHECK_NEAR(result.total_generation_mw(), 272.39327236, 1e-5);
    CHECK_NEAR(result.total_generation_mvar(), 82.43754387, 1e-5);
    CHECK_NEAR(shunt.imag() * base, -21.18484418, 1e-5);
}

TEST(PowerFlow, Ieee14BusVoltagesMatchReference) {
    const double expected_vm[14] = {1.0600000000, 1.0450000000, 1.0100000000, 1.0176708537,
                                    1.0195138598, 1.0700000000, 1.0615195325, 1.0900000000,
                                    1.0559317206, 1.0509846250, 1.0569065185, 1.0551885632,
                                    1.0503817136, 1.0355299459};
    const double expected_va[14] = {0.0000000000,  -4.9825891420, -12.7250999383, -10.3129010923,
                                    -8.7738538983, -14.2209464637, -13.3596273653, -13.3596273653,
                                    -14.9385212952, -15.0972884631, -14.7906220313, -15.0755845204,
                                    -15.1562763362, -16.0336445292};

    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());

    for (std::size_t i = 0; i < 14; ++i) {
        CHECK_NEAR(result.vm[i], expected_vm[i], 1e-8);
        CHECK_NEAR(result.va[i], expected_va[i], 1e-8);
    }
}

TEST(PowerFlow, Ieee14BranchFlowsMatchReference) {
    // Flow leaving the "from" terminal, MW / MVAr.
    const double expected_p[5] = {156.88289053, 75.51038183, 73.23757924, 56.13149594,
                                  41.51621502};
    const double expected_q[5] = {-20.40429168, 3.85499114, 3.56020295, -1.55035040, 1.17099786};

    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());

    for (std::size_t k = 0; k < 5; ++k) {
        CHECK_NEAR(result.sf[k].real() * 100.0, expected_p[k], 1e-5);
        CHECK_NEAR(result.sf[k].imag() * 100.0, expected_q[k], 1e-5);
    }
}

TEST(PowerFlow, BusInjectionEqualsBranchFlowsPlusShunts) {
    // Conservation check that ties the Y-bus, the branch model and the branch
    // flows together: the sum of the bus injections must equal the sum of the
    // branch losses plus the power absorbed by the bus shunts.
    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());

    pf::Complex bus_total{0.0, 0.0};
    for (const pf::Complex& s : result.s_inj) {
        bus_total += s;
    }

    pf::Complex branch_total{0.0, 0.0};
    for (std::size_t k = 0; k < net.branch_count(); ++k) {
        branch_total += result.sf[k] + result.st[k];
    }

    pf::Complex shunt_total{0.0, 0.0};
    for (std::size_t i = 0; i < net.bus_count(); ++i) {
        const double magnitude = result.vm[i] * result.vm[i];
        shunt_total += magnitude * std::conj(net.bus_shunts()[i]);
    }

    CHECK_NEAR(std::abs(bus_total - branch_total - shunt_total), 0.0, 1e-12);
}

TEST(PowerFlow, JacobianMatchesFiniteDifferences) {
    // The analytic Jacobian is the single most error-prone piece of algebra in
    // the solver, so compare every entry against a central difference of the
    // mismatch function at a realistic operating point.
    const pf::Network net = load_network("case14.m");
    const std::size_t nb = net.bus_count();

    std::vector<double> vm(nb, 0.0);
    std::vector<double> va(nb, 0.0);
    for (std::size_t i = 0; i < nb; ++i) {
        vm[i] = net.bus(i).vm;
        va[i] = net.bus(i).va * pf::kDegToRad;
    }

    const std::vector<pf::Index>& non_slack = net.non_slack_buses();
    const std::vector<pf::Index>& pq = net.pq_buses();
    const std::size_t n_theta = non_slack.size();
    const std::size_t n = n_theta + pq.size();

    const pf::RealMatrix jac = pf::build_power_flow_jacobian(net.ybus(), vm, va, non_slack, pq);
    CHECK_EQ(jac.rows(), n);
    CHECK_EQ(jac.cols(), n);

    // The analytic Jacobian is d(P_calc, Q_calc)/d(theta, V); its entries are the
    // derivatives of the *calculated* injections, so the residual used for the
    // finite differences must be the calculated injection itself, not a mismatch.
    const auto residual = [&](const std::vector<double>& vm_x, const std::vector<double>& va_x) {
        std::vector<pf::Complex> v(nb);
        for (std::size_t i = 0; i < nb; ++i) {
            v[i] = std::polar(vm_x[i], va_x[i]);
        }
        const std::vector<pf::Complex> current = pf::matvec(net.ybus(), v);
        std::vector<double> values(n, 0.0);
        for (std::size_t k = 0; k < n_theta; ++k) {
            const pf::Complex s = v[non_slack[k]] * std::conj(current[non_slack[k]]);
            values[k] = s.real();
        }
        for (std::size_t k = 0; k < pq.size(); ++k) {
            const pf::Complex s = v[pq[k]] * std::conj(current[pq[k]]);
            values[n_theta + k] = s.imag();
        }
        return values;
    };

    std::vector<double> x(n, 0.0);
    for (std::size_t k = 0; k < n_theta; ++k) {
        x[k] = va[non_slack[k]];
    }
    for (std::size_t k = 0; k < pq.size(); ++k) {
        x[n_theta + k] = vm[pq[k]];
    }

    const auto state_of = [&](const std::vector<double>& xv) {
        std::vector<double> vm_out = vm;
        std::vector<double> va_out = va;
        for (std::size_t k = 0; k < n_theta; ++k) {
            va_out[non_slack[k]] = xv[k];
        }
        for (std::size_t k = 0; k < pq.size(); ++k) {
            vm_out[pq[k]] = xv[n_theta + k];
        }
        return std::make_pair(vm_out, va_out);
    };

    const double step = 1e-6;
    for (std::size_t c = 0; c < n; ++c) {
        std::vector<double> plus = x;
        std::vector<double> minus = x;
        plus[c] += step;
        minus[c] -= step;

        const auto forward = state_of(plus);
        const auto backward = state_of(minus);
        const std::vector<double> f_plus = residual(forward.first, forward.second);
        const std::vector<double> f_minus = residual(backward.first, backward.second);

        for (std::size_t r = 0; r < n; ++r) {
            const double numeric = (f_plus[r] - f_minus[r]) / (2.0 * step);
            const double tolerance = 2e-5 * std::max(1.0, std::abs(numeric));
            CHECK_NEAR(jac(r, c), numeric, tolerance);
        }
    }
}

TEST(PowerFlow, SparseAndDenseLinearSolversAgree) {
    // The two factorisations are completely different code paths, so matching
    // results to the last digit is a strong check on both.
    for (const char* name : {"case14.m", "case57.m", "case118.m"}) {
        const pf::Network net = load_network(name);

        pf::SolverOptions sparse_options = tight_options();
        sparse_options.linear_solver = pf::LinearSolver::Sparse;
        const pf::PowerFlowResult sparse = pf::solve_newton_raphson(net, sparse_options);

        pf::SolverOptions dense_options = tight_options();
        dense_options.linear_solver = pf::LinearSolver::Dense;
        const pf::PowerFlowResult dense = pf::solve_newton_raphson(net, dense_options);

        CHECK(sparse.status.converged);
        CHECK(dense.status.converged);
        CHECK_EQ(sparse.status.iterations, dense.status.iterations);
        CHECK_NEAR(sparse.p_loss_mw, dense.p_loss_mw, 1e-9);
        CHECK_NEAR(sparse.q_loss_mvar, dense.q_loss_mvar, 1e-6);
        for (std::size_t i = 0; i < net.bus_count(); ++i) {
            CHECK_NEAR(sparse.vm[i], dense.vm[i], 1e-10);
            CHECK_NEAR(sparse.va[i], dense.va[i], 1e-8);
        }
    }
}

TEST(PowerFlow, SparsePathReportsACompactFactorisation) {
    const pf::Network net = load_network("case118.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());

    CHECK(result.linear_solver == pf::LinearSolver::Sparse);
    CHECK(result.jacobian_nonzeros > 0);
    // The factors stay far below the dense n^2 the old path would have stored.
    const std::size_t n = net.non_slack_buses().size() + net.pq_buses().size();
    CHECK(result.factor_nonzeros < n * n / 2);
}

TEST(PowerFlow, GaussSeidelAgreesWithNewton) {
    const pf::Network net = load_network("case14.m");

    pf::SolverOptions options = tight_options();
    options.max_iterations = 20000;
    const pf::PowerFlowResult gauss = pf::solve_gauss_seidel(net, options);
    const pf::PowerFlowResult newton = pf::solve_newton_raphson(net, tight_options());

    CHECK(gauss.status.converged);
    CHECK(gauss.status.iterations > newton.status.iterations);
    for (std::size_t i = 0; i < net.bus_count(); ++i) {
        CHECK_NEAR(gauss.vm[i], newton.vm[i], 1e-6);
        CHECK_NEAR(gauss.va[i], newton.va[i], 1e-5);
    }
    CHECK_NEAR(gauss.p_loss_mw, newton.p_loss_mw, 1e-4);
    CHECK_NEAR(gauss.q_loss_mvar, newton.q_loss_mvar, 1e-3);
}

TEST(PowerFlow, ReportsFailureWhenTheBudgetRunsOut) {
    const pf::Network net = load_network("case14.m");
    pf::SolverOptions options;
    options.max_iterations = 1;
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, options);
    CHECK_FALSE(result.status.converged);
    CHECK(!result.status.message.empty());
}

TEST(PowerFlow, Ieee14HasNoRatedBranches) {
    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());
    CHECK_EQ(result.worst_loaded_branch(), pf::kInvalidIndex);
    CHECK_NEAR(result.min_voltage(), 1.0100000000, 1e-8);
    CHECK_NEAR(result.max_voltage(), 1.0900000000, 1e-8);
}

TEST(PowerFlow, Ieee118HandlesReactiveLimits) {
    const pf::Network net = load_network("case118.m");
    pf::SolverOptions options = tight_options();
    options.enforce_q_limits = true;
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, options);

    CHECK(result.status.converged);
    CHECK(result.q_limit_switches > 0);
    // Losses are not monotone in the reactive dispatch, but the constraint must
    // actually hold afterwards: every generator bus except the slack has to sit
    // inside its reactive limits.
    const double base = net.base_mva();
    for (std::size_t i = 0; i < net.bus_count(); ++i) {
        if (i == net.slack_bus()) {
            continue;
        }
        double qmax = 0.0;
        double qmin = 0.0;
        bool present = false;
        for (const pf::Generator& gen : net.data().generators) {
            if (gen.online && gen.bus == net.bus(i).id) {
                qmax += gen.qmax;
                qmin += gen.qmin;
                present = true;
            }
        }
        if (!present) {
            continue;
        }
        const double qg = result.s_gen[i].imag() * base;
        CHECK(qg <= qmax + 1e-3);
        CHECK(qg >= qmin - 1e-3);
    }
}

TEST(PowerFlow, Ieee14HasSlackOnlyGenerationOnBusOne) {
    const pf::Network net = load_network("case14.m");
    const pf::PowerFlowResult result = pf::solve_newton_raphson(net, tight_options());

    // Bus 1 picks up whatever the rest of the system does not supply.
    const double expected_slack =
        result.total_generation_mw() - 40.0 - 0.0 - 0.0 - 0.0;
    CHECK_NEAR(result.s_gen[0].real() * 100.0, expected_slack, 1e-6);

    // Every non-generator bus must end up with zero generation.
    for (std::size_t i = 0; i < net.bus_count(); ++i) {
        const bool has_generator = net.bus(i).id == 1 || net.bus(i).id == 2 ||
                                   net.bus(i).id == 3 || net.bus(i).id == 6 ||
                                   net.bus(i).id == 8;
        if (!has_generator) {
            CHECK_NEAR(std::abs(result.s_gen[i]), 0.0, 1e-10);
        }
    }
}
