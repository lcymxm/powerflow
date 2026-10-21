#include "powerflow/dc_opf.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "powerflow/dense.hpp"

namespace pf {
namespace {

struct DcModel {
    std::size_t buses = 0;
    std::size_t generators = 0;
    std::vector<std::size_t> generator_index;  // into net.data().generators
    DenseMatrix<double> bbus;                  // MW per radian
    DenseMatrix<double> flow_row;              // MW per radian, one row per branch
    std::vector<double> flow_offset;           // MW, from phase shifters
    std::vector<char> rated;
    std::vector<double> rating;     // MW
    std::vector<double> demand_mw;  // per bus
    std::vector<double> offset_mw;  // injection-like term from phase shifters
    std::size_t slack = 0;
};

DcModel build_model(const Network& net) {
    const CaseData& data = net.data();
    const double base = net.base_mva();

    DcModel model;
    model.buses = net.bus_count();
    model.slack = net.slack_bus();
    model.bbus = DenseMatrix<double>(model.buses, model.buses, 0.0);
    model.flow_row = DenseMatrix<double>(data.branches.size(), model.buses, 0.0);
    model.flow_offset.assign(data.branches.size(), 0.0);
    model.rated.assign(data.branches.size(), 0);
    model.rating.assign(data.branches.size(), 0.0);
    model.demand_mw.assign(model.buses, 0.0);
    model.offset_mw.assign(model.buses, 0.0);

    for (Index i = 0; i < model.buses; ++i) {
        model.demand_mw[i] = net.bus(i).pd;
    }

    for (Index k = 0; k < data.branches.size(); ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        const Branch& br = data.branches[k];
        if (br.x == 0.0) {
            throw std::runtime_error("dc opf: branch " + std::to_string(k + 1) +
                                     " has zero reactance");
        }
        const double tap = br.ratio != 0.0 ? br.ratio : 1.0;
        const double susceptance_mw = base / (br.x * tap);  // MW per radian
        const Index f = net.index_of(br.from);
        const Index t = net.index_of(br.to);

        model.flow_row(k, f) += susceptance_mw;
        model.flow_row(k, t) -= susceptance_mw;
        model.flow_offset[k] = -susceptance_mw * br.shift * kDegToRad;

        model.bbus(f, f) += susceptance_mw;
        model.bbus(f, t) -= susceptance_mw;
        model.bbus(t, f) -= susceptance_mw;
        model.bbus(t, t) += susceptance_mw;

        model.offset_mw[f] += model.flow_offset[k];
        model.offset_mw[t] -= model.flow_offset[k];

        if (br.rate_a > 0.0) {
            model.rated[k] = 1;
            model.rating[k] = br.rate_a;
        }
    }

    for (Index g = 0; g < data.generators.size(); ++g) {
        if (data.generators[g].online) {
            model.generator_index.push_back(g);
        }
    }
    model.generators = model.generator_index.size();
    return model;
}

std::vector<double> multiply(const DenseMatrix<double>& matrix, const std::vector<double>& x) {
    return matvec(matrix, x);
}

std::vector<double> transpose_multiply(const DenseMatrix<double>& matrix,
                                       const std::vector<double>& y) {
    std::vector<double> out(matrix.cols(), 0.0);
    for (std::size_t r = 0; r < matrix.rows(); ++r) {
        const double value = y[r];
        if (value == 0.0) {
            continue;
        }
        for (std::size_t c = 0; c < matrix.cols(); ++c) {
            out[c] += matrix(r, c) * value;
        }
    }
    return out;
}

double norm_inf(const std::vector<double>& values) {
    double worst = 0.0;
    for (const double value : values) {
        worst = std::max(worst, std::abs(value));
    }
    return worst;
}

// Small proximal term added to the objective, 1/2 eps ||x||^2, included in both
// the KKT matrix and the stationarity residual. The formulation below has no
// angle variables, so the Hessian is already positive definite whenever a
// generator has a quadratic cost; the term only guards the purely linear ones.
constexpr double kProximal = 1e-6;

// Solves the KKT system after equilibrating rows and columns.
//
// Near convergence the diagonal weights D = z / s span many orders of magnitude
// (an active limit has s -> 0 while its multiplier stays finite), so the raw KKT
// matrix is far too badly scaled for a relative singularity test: partial
// pivoting reports it singular when it is merely stiff. Scaling every row and
// column to unit inf-norm fixes both the false singularity and the accuracy of
// the step.
class EquilibratedSolver {
public:
    explicit EquilibratedSolver(const DenseMatrix<double>& matrix)
        : size_(matrix.rows()), row_scale_(matrix.rows(), 1.0), column_scale_(matrix.cols(), 1.0) {
        DenseMatrix<double> scaled = matrix;
        for (std::size_t r = 0; r < size_; ++r) {
            double largest = 0.0;
            for (std::size_t c = 0; c < size_; ++c) {
                largest = std::max(largest, std::abs(matrix(r, c)));
            }
            row_scale_[r] = largest > 0.0 ? 1.0 / largest : 1.0;
        }
        for (std::size_t c = 0; c < size_; ++c) {
            double largest = 0.0;
            for (std::size_t r = 0; r < size_; ++r) {
                largest = std::max(largest, std::abs(matrix(r, c)) * row_scale_[r]);
            }
            column_scale_[c] = largest > 0.0 ? 1.0 / largest : 1.0;
        }
        for (std::size_t r = 0; r < size_; ++r) {
            for (std::size_t c = 0; c < size_; ++c) {
                scaled(r, c) = matrix(r, c) * row_scale_[r] * column_scale_[c];
            }
        }
        solver_ = std::make_unique<LuSolver<double>>(scaled);
    }

    bool is_singular() const { return solver_->is_singular(); }

    std::vector<double> solve(const std::vector<double>& rhs) const {
        std::vector<double> scaled_rhs(size_, 0.0);
        for (std::size_t i = 0; i < size_; ++i) {
            scaled_rhs[i] = rhs[i] * row_scale_[i];
        }
        const std::vector<double> scaled = solver_->solve(scaled_rhs);
        std::vector<double> out(size_, 0.0);
        for (std::size_t i = 0; i < size_; ++i) {
            out[i] = scaled[i] * column_scale_[i];
        }
        return out;
    }

private:
    std::size_t size_ = 0;
    std::vector<double> row_scale_;
    std::vector<double> column_scale_;
    std::unique_ptr<LuSolver<double>> solver_;
};

// Primal-dual interior point for
//     minimise  1/2 x'Qx + c'x    subject to  Ax = b,  Gx <= h
// with Q diagonal and positive semidefinite. Mehrotra predictor-corrector.
//
// Each iteration solves the full KKT system
//     [ M  -A' ] [dx]   [rhs ]
//     [ A   0  ] [dy] = [-r_p],      M = Q + G' D G
// rather than eliminating dx first. Forming A M^-1 A' would square the condition
// number, and M is only positive *semi* definite when no inequality touches the
// angle variables, so the KKT form is both more accurate and more robust.
struct QpProblem {
    std::vector<double> q;  // diagonal of Q
    std::vector<double> c;
    DenseMatrix<double> a;
    std::vector<double> b;
    DenseMatrix<double> g;
    std::vector<double> h;
    // Optional starting point. An interior one matters a lot: starting at zero
    // can violate the lower generator bounds and stall the method immediately.
    std::vector<double> x_start;
};

struct QpSolution {
    bool converged = false;
    std::string message;
    int iterations = 0;
    std::vector<double> x;
    std::vector<double> y;  // equality duals
    std::vector<double> z;  // inequality duals
    double complementarity = 0.0;
    double equality_error = 0.0;
    double inequality_violation = 0.0;
};

QpSolution solve_qp(const QpProblem& problem, double tolerance, int max_iterations) {
    const std::size_t n = problem.q.size();
    const std::size_t m_eq = problem.a.rows();
    const std::size_t m_ineq = problem.g.rows();
    const std::size_t total = n + m_eq;

    QpSolution solution;
    solution.x.assign(n, 0.0);
    if (problem.x_start.size() == n) {
        solution.x = problem.x_start;
    }
    solution.y.assign(m_eq, 0.0);

    // A start strictly inside the inequalities: flat angles give zero branch
    // flows and the generator start points sit inside their limits.
    std::vector<double> s(m_ineq, 0.0);
    std::vector<double> z(m_ineq, 1.0);
    for (std::size_t i = 0; i < m_ineq; ++i) {
        double slack = problem.h[i];
        for (std::size_t j = 0; j < n; ++j) {
            slack -= problem.g(i, j) * solution.x[j];
        }
        // A bound that is active at the start would make D = z / s infinite.
        s[i] = std::max(slack, 1.0);
    }

    std::vector<double> r_d(n, 0.0);
    std::vector<double> r_p(m_eq, 0.0);
    std::vector<double> r_h(m_ineq, 0.0);

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        solution.iterations = iteration;

        const std::vector<double> ax = multiply(problem.a, solution.x);
        for (std::size_t i = 0; i < m_eq; ++i) {
            r_p[i] = ax[i] - problem.b[i];
        }
        const std::vector<double> gx = multiply(problem.g, solution.x);
        for (std::size_t i = 0; i < m_ineq; ++i) {
            r_h[i] = gx[i] + s[i] - problem.h[i];
        }
        const std::vector<double> aty = transpose_multiply(problem.a, solution.y);
        const std::vector<double> gtz = transpose_multiply(problem.g, z);
        for (std::size_t j = 0; j < n; ++j) {
            r_d[j] = (problem.q[j] + kProximal) * solution.x[j] + problem.c[j] - aty[j] - gtz[j];
        }

        double mu = 0.0;
        for (std::size_t i = 0; i < m_ineq; ++i) {
            mu += s[i] * z[i];
        }
        mu /= static_cast<double>(std::max<std::size_t>(m_ineq, 1));
        solution.complementarity = mu;
        solution.equality_error = norm_inf(r_p);
        solution.inequality_violation = std::max(0.0, norm_inf(r_h));

        if (mu < tolerance && norm_inf(r_p) < 1e-7 && norm_inf(r_d) < 1e-7) {
            solution.converged = true;
            solution.message = "converged";
            return solution;
        }

        // M = Q + G' D G changes every iteration through D = z / s.
        DenseMatrix<double> kkt(total, total, 0.0);
        for (std::size_t j = 0; j < n; ++j) {
            kkt(j, j) = std::max(problem.q[j], 0.0) + kProximal;
        }        for (std::size_t i = 0; i < m_ineq; ++i) {
            const double weight = std::min(z[i] / s[i], 1e10);
            if (weight == 0.0) {
                continue;
            }
            for (std::size_t r = 0; r < n; ++r) {
                const double gr = problem.g(i, r);
                if (gr == 0.0) {
                    continue;
                }
                for (std::size_t c = 0; c < n; ++c) {
                    const double gc = problem.g(i, c);
                    if (gc != 0.0) {
                        kkt(r, c) += weight * gr * gc;
                    }
                }
            }
        }
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t i = 0; i < m_eq; ++i) {
                kkt(r, n + i) = -problem.a(i, r);
                kkt(n + i, r) = problem.a(i, r);
            }
        }

        const EquilibratedSolver kkt_solver(kkt);
        if (kkt_solver.is_singular()) {
            solution.message = "the KKT system became singular";
            return solution;
        }
        const auto complementarity_rhs = [&](std::size_t i, double sigma, double mu_current,
                                             const std::vector<double>& ds_affine,
                                             const std::vector<double>& dz_affine) {
            double rc = s[i] * z[i] - sigma * mu_current;
            if (!ds_affine.empty()) {
                rc += ds_affine[i] * dz_affine[i];
            }
            return rc;
        };

        const auto direction = [&](double sigma, double mu_current,
                                   const std::vector<double>& ds_affine,
                                   const std::vector<double>& dz_affine, std::vector<double>& dx,
                                   std::vector<double>& dy, std::vector<double>& ds,
                                   std::vector<double>& dz) {
            std::vector<double> rhs(total, 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                rhs[j] = -r_d[j];
            }
            for (std::size_t i = 0; i < m_eq; ++i) {
                rhs[n + i] = -r_p[i];
            }
            std::vector<double> g_term(m_ineq, 0.0);
            for (std::size_t i = 0; i < m_ineq; ++i) {
                const double rc = complementarity_rhs(i, sigma, mu_current, ds_affine, dz_affine);
                g_term[i] = (z[i] / s[i]) * r_h[i] - rc / s[i];
            }
            const std::vector<double> gt = transpose_multiply(problem.g, g_term);
            for (std::size_t j = 0; j < n; ++j) {
                rhs[j] += gt[j];
            }

            const std::vector<double> step = kkt_solver.solve(rhs);
            dx.assign(step.begin(), step.begin() + static_cast<std::ptrdiff_t>(n));
            dy.assign(step.begin() + static_cast<std::ptrdiff_t>(n), step.end());

            ds.assign(m_ineq, 0.0);
            for (std::size_t i = 0; i < m_ineq; ++i) {
                double value = -r_h[i];
                for (std::size_t j = 0; j < n; ++j) {
                    value -= problem.g(i, j) * dx[j];
                }
                ds[i] = value;
            }
            dz.assign(m_ineq, 0.0);
            for (std::size_t i = 0; i < m_ineq; ++i) {
                const double rc = complementarity_rhs(i, sigma, mu_current, ds_affine, dz_affine);
                dz[i] = -rc / s[i] - (z[i] / s[i]) * ds[i];
            }
        };

        const auto step_length = [](const std::vector<double>& value,
                                    const std::vector<double>& direction_vector) {
            double alpha = 1.0;
            for (std::size_t i = 0; i < value.size(); ++i) {
                if (direction_vector[i] < 0.0) {
                    alpha = std::min(alpha, -value[i] / direction_vector[i]);
                }
            }
            return std::min(1.0, 0.995 * alpha);
        };

        std::vector<double> dx;
        std::vector<double> dy;
        std::vector<double> ds;
        std::vector<double> dz;
        const std::vector<double> empty;

        direction(0.0, mu, empty, empty, dx, dy, ds, dz);
        const double alpha_affine = std::min(step_length(s, ds), step_length(z, dz));

        // Copy before the corrector runs: the direction routine clears its output
        // vectors, and these double as the affine inputs for the second-order term.
        const std::vector<double> ds_affine = ds;
        const std::vector<double> dz_affine = dz;

        double mu_affine = 0.0;
        for (std::size_t i = 0; i < m_ineq; ++i) {
            mu_affine += (s[i] + alpha_affine * ds[i]) * (z[i] + alpha_affine * dz[i]);
        }
        mu_affine /= static_cast<double>(std::max<std::size_t>(m_ineq, 1));
        double sigma = 0.0;
        if (mu > 0.0) {
            sigma = std::clamp(std::pow(mu_affine / mu, 3.0), 0.0, 1.0);
        }

        direction(sigma, mu, ds_affine, dz_affine, dx, dy, ds, dz);
        const double alpha_primal = step_length(s, ds);
        const double alpha_dual = step_length(z, dz);

        for (std::size_t j = 0; j < n; ++j) {
            solution.x[j] += alpha_primal * dx[j];
        }
        for (std::size_t i = 0; i < m_eq; ++i) {
            solution.y[i] += alpha_dual * dy[i];
        }
        for (std::size_t i = 0; i < m_ineq; ++i) {
            s[i] = std::max(s[i] + alpha_primal * ds[i], 1e-300);
            z[i] = std::max(z[i] + alpha_dual * dz[i], 1e-300);
        }

    }

    solution.message = "reached the iteration limit";
    solution.z = z;
    return solution;
}

}  // namespace

DcOpfResult solve_dc_opf(const Network& net, const DcOpfOptions& options) {
    const CaseData& data = net.data();
    if (!data.has_generation_costs) {
        throw std::runtime_error(
            "dc opf: the case has no polynomial generator costs (mpc.gencost)");
    }

    const DcModel model = build_model(net);
    const std::size_t nb = model.buses;
    const std::size_t ng = model.generators;
    const std::size_t nl = data.branches.size();
    if (ng == 0) {
        throw std::runtime_error("dc opf: the case has no online generators");
    }

    // ---- power transfer distribution factors ---------------------------------
    // The angles are eliminated: a unit injection at bus i produces exactly the
    // branch flows stored in column i of `ptdf`. Every network constraint then
    // becomes a linear constraint on the generation alone, which keeps the
    // quadratic program small and well conditioned even when no branch has a
    // rating (the degenerate case that made an angle-based formulation singular).
    std::vector<Index> reduced_bus;
    std::vector<Index> reduced_index(nb, kInvalidIndex);
    reduced_bus.reserve(nb - 1);
    for (Index i = 0; i < nb; ++i) {
        if (i == model.slack) {
            continue;
        }
        reduced_index[i] = reduced_bus.size();
        reduced_bus.push_back(i);
    }
    const std::size_t nr = reduced_bus.size();
    if (nr == 0) {
        throw std::runtime_error("dc opf: the case has a single bus");
    }

    DenseMatrix<double> laplacian(nr, nr, 0.0);
    for (std::size_t r = 0; r < nr; ++r) {
        for (std::size_t c = 0; c < nr; ++c) {
            laplacian(r, c) = model.bbus(reduced_bus[r], reduced_bus[c]);
        }
    }
    const LuSolver<double> laplacian_solver(laplacian);
    if (laplacian_solver.is_singular()) {
        throw std::runtime_error("dc opf: the network is not connected");
    }

    std::vector<double> ptdf(nl * nb, 0.0);
    std::vector<double> rhs(nr, 0.0);
    std::vector<double> theta(nb, 0.0);
    for (const Index bus : reduced_bus) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        rhs[reduced_index[bus]] = 1.0;
        const std::vector<double> solved = laplacian_solver.solve(rhs);
        std::fill(theta.begin(), theta.end(), 0.0);
        for (std::size_t r = 0; r < nr; ++r) {
            theta[reduced_bus[r]] = solved[r];
        }
        for (std::size_t k = 0; k < nl; ++k) {
            if (!net.branch_in_service(k)) {
                continue;
            }
            double flow = 0.0;
            for (Index j = 0; j < nb; ++j) {
                flow += model.flow_row(k, j) * theta[j];
            }
            ptdf[k * nb + bus] = flow;
        }
    }

    // flow_k(Pg) = base_flow[k] + sum_g coefficient(k, g) * Pg_g
    std::vector<double> base_flow(nl, 0.0);
    std::vector<double> coefficient(nl * ng, 0.0);
    for (std::size_t k = 0; k < nl; ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        double flow = model.flow_offset[k];
        for (Index i = 0; i < nb; ++i) {
            flow -= ptdf[k * nb + i] * (model.demand_mw[i] + model.offset_mw[i]);
        }
        base_flow[k] = flow;
        for (std::size_t g = 0; g < ng; ++g) {
            const std::size_t bus = net.index_of(data.generators[model.generator_index[g]].bus);
            coefficient[k * ng + g] = ptdf[k * nb + bus];
        }
    }

    // ---- if the network-free dispatch is feasible, it is already optimal -----
    // The economic dispatch is the unconstrained lower bound on the cost. When it
    // also respects every branch rating, it *is* the OPF solution and every nodal
    // price equals the system marginal cost. Checking that first is exact, and it
    // avoids running a barrier method on the degenerate case where several units
    // share identical cost coefficients: there the optimal dispatch is a whole
    // face, the duals are not unique, and a KKT solve has nothing to converge to.
    {
        const DcOpfResult dispatch = solve_economic_dispatch(net, options);
        if (dispatch.converged && dispatch.pg_mw.size() == ng) {
            bool feasible = true;
            for (std::size_t g = 0; g < ng && feasible; ++g) {
                const Generator& generator = data.generators[model.generator_index[g]];
                const double upper = generator.pmax > 0.0 ? generator.pmax : 0.0;
                const double pg = dispatch.pg_mw[g];
                feasible = pg >= generator.pmin - 1e-6 && pg <= upper + 1e-6;
            }
            for (std::size_t k = 0; k < nl && feasible; ++k) {
                if (!net.branch_in_service(k) || model.rated[k] == 0) {
                    continue;
                }
                double flow = base_flow[k];
                for (std::size_t g = 0; g < ng; ++g) {
                    flow += coefficient[k * ng + g] * dispatch.pg_mw[g];
                }
                feasible = std::abs(flow) <= model.rating[k] + 1e-6;
            }

            if (feasible) {
                DcOpfResult result;
                result.converged = true;
                result.message = "uncongested: the economic dispatch respects every rating";
                result.iterations = 0;
                result.objective = dispatch.objective;
                result.total_load_mw = dispatch.total_load_mw;
                result.total_generation_mw = dispatch.total_generation_mw;
                result.max_balance_error_mw =
                    std::abs(dispatch.total_generation_mw - dispatch.total_load_mw);
                result.generators = model.generator_index;
                result.pg_mw = dispatch.pg_mw;
                const double lambda = dispatch.lmp.empty() ? 0.0 : dispatch.lmp[0];
                result.lmp.assign(nb, lambda);
                result.branch_flow_mw.assign(nl, 0.0);
                result.branch_congested.assign(nl, 0);

                std::vector<double> angle_rhs(nr, 0.0);
                std::vector<double> injection(nb, 0.0);
                for (Index i = 0; i < nb; ++i) {
                    injection[i] = -model.demand_mw[i] - model.offset_mw[i];
                }
                for (std::size_t g = 0; g < ng; ++g) {
                    const std::size_t bus =
                        net.index_of(data.generators[model.generator_index[g]].bus);
                    injection[bus] += dispatch.pg_mw[g];
                }
                for (std::size_t r = 0; r < nr; ++r) {
                    angle_rhs[r] = injection[reduced_bus[r]];
                }
                const std::vector<double> angles = laplacian_solver.solve(angle_rhs);
                result.theta_deg.assign(nb, 0.0);
                for (std::size_t r = 0; r < nr; ++r) {
                    result.theta_deg[reduced_bus[r]] = angles[r] * kRadToDeg;
                }
                for (std::size_t k = 0; k < nl; ++k) {
                    if (!net.branch_in_service(k)) {
                        continue;
                    }
                    double flow = base_flow[k];
                    for (std::size_t g = 0; g < ng; ++g) {
                        flow += coefficient[k * ng + g] * dispatch.pg_mw[g];
                    }
                    result.branch_flow_mw[k] = flow;
                }
                return result;
            }
        }
    }

    // ---- quadratic program in the generator outputs only ---------------------
    QpProblem problem;
    problem.q.assign(ng, 0.0);
    problem.c.assign(ng, 0.0);
    problem.x_start.assign(ng, 0.0);
    for (std::size_t g = 0; g < ng; ++g) {
        const Generator& generator = data.generators[model.generator_index[g]];
        problem.q[g] = 2.0 * generator.cost_c2;
        problem.c[g] = generator.cost_c1;
        const double pmax = generator.pmax > 0.0 ? generator.pmax : 0.0;
        problem.x_start[g] = 0.5 * (generator.pmin + pmax);
    }

    // Equality: the units must cover the load.
    problem.a = DenseMatrix<double>(1, ng, 0.0);
    problem.b.assign(1, 0.0);
    double total_demand = 0.0;
    for (const double demand : model.demand_mw) {
        total_demand += demand;
    }
    for (std::size_t g = 0; g < ng; ++g) {
        problem.a(0, g) = 1.0;
    }
    problem.b[0] = total_demand;

    std::size_t rated_count = 0;
    for (std::size_t k = 0; k < nl; ++k) {
        if (net.branch_in_service(k) && model.rated[k] != 0) {
            ++rated_count;
        }
    }
    const std::size_t m_ineq = 2 * ng + 2 * rated_count;
    problem.g = DenseMatrix<double>(m_ineq, ng, 0.0);
    problem.h.assign(m_ineq, 0.0);

    std::size_t row = 0;
    for (std::size_t g = 0; g < ng; ++g) {
        const Generator& generator = data.generators[model.generator_index[g]];
        const double pmax = generator.pmax > 0.0 ? generator.pmax : 0.0;
        problem.g(row, g) = 1.0;
        problem.h[row] = pmax;
        ++row;
        problem.g(row, g) = -1.0;
        problem.h[row] = -generator.pmin;
        ++row;
    }

    // Which inequality row holds which flow limit, so the duals can be turned
    // back into nodal prices afterwards.
    std::vector<std::size_t> upper_row(nl, std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> lower_row(nl, std::numeric_limits<std::size_t>::max());
    for (std::size_t k = 0; k < nl; ++k) {
        if (!net.branch_in_service(k) || model.rated[k] == 0) {
            continue;
        }
        for (std::size_t g = 0; g < ng; ++g) {
            problem.g(row, g) = coefficient[k * ng + g];
            problem.g(row + 1, g) = -coefficient[k * ng + g];
        }
        problem.h[row] = model.rating[k] - base_flow[k];
        problem.h[row + 1] = model.rating[k] + base_flow[k];
        upper_row[k] = row;
        lower_row[k] = row + 1;
        row += 2;
    }

    const QpSolution qp = solve_qp(problem, options.tolerance, options.max_iterations);

    DcOpfResult result;
    result.converged = qp.converged;
    result.message = qp.message;
    result.iterations = qp.iterations;
    result.complementarity = qp.complementarity;
    result.max_balance_error_mw = qp.equality_error;
    result.max_limit_violation_mw = qp.inequality_violation;
    result.generators = model.generator_index;
    result.total_load_mw = total_demand;

    result.pg_mw.assign(ng, 0.0);
    for (std::size_t g = 0; g < ng; ++g) {
        const Generator& generator = data.generators[model.generator_index[g]];
        result.pg_mw[g] = qp.x[g];
        result.total_generation_mw += result.pg_mw[g];
        result.objective += generator.cost_c2 * result.pg_mw[g] * result.pg_mw[g] +
                            generator.cost_c1 * result.pg_mw[g] + generator.cost_c0;
    }

    // Nodal prices: the balance dual plus the congestion rent each branch pushes
    // onto every bus through its PTDF column.
    const double lambda = qp.y.empty() ? 0.0 : qp.y[0];
    result.lmp.assign(nb, lambda);
    for (std::size_t k = 0; k < nl; ++k) {
        if (upper_row[k] == std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        const double mu_upper = upper_row[k] < qp.z.size() ? qp.z[upper_row[k]] : 0.0;
        const double mu_lower = lower_row[k] < qp.z.size() ? qp.z[lower_row[k]] : 0.0;
        const double rent = mu_upper - mu_lower;
        if (rent == 0.0) {
            continue;
        }
        for (Index i = 0; i < nb; ++i) {
            result.lmp[i] -= ptdf[k * nb + i] * rent;
        }
    }

    // Recover the angles and the branch flows for reporting.
    std::vector<double> injection(nb, 0.0);
    for (Index i = 0; i < nb; ++i) {
        injection[i] = -model.demand_mw[i] - model.offset_mw[i];
    }
    for (std::size_t g = 0; g < ng; ++g) {
        const std::size_t bus = net.index_of(data.generators[model.generator_index[g]].bus);
        injection[bus] += qp.x[g];
    }
    std::fill(rhs.begin(), rhs.end(), 0.0);
    for (std::size_t r = 0; r < nr; ++r) {
        rhs[r] = injection[reduced_bus[r]];
    }
    const std::vector<double> angles = laplacian_solver.solve(rhs);
    result.theta_deg.assign(nb, 0.0);
    for (std::size_t r = 0; r < nr; ++r) {
        result.theta_deg[reduced_bus[r]] = angles[r] * kRadToDeg;
    }

    result.branch_flow_mw.assign(nl, 0.0);
    result.branch_congested.assign(nl, 0);
    for (std::size_t k = 0; k < nl; ++k) {
        if (!net.branch_in_service(k)) {
            continue;
        }
        double flow = base_flow[k];
        for (std::size_t g = 0; g < ng; ++g) {
            flow += coefficient[k * ng + g] * qp.x[g];
        }
        result.branch_flow_mw[k] = flow;
        if (model.rated[k] != 0 && std::abs(flow) >= 0.999 * model.rating[k]) {
            result.branch_congested[k] = 1;
        }
    }
    return result;
}

DcOpfResult solve_economic_dispatch(const Network& net, const DcOpfOptions& options) {
    const CaseData& data = net.data();
    if (!data.has_generation_costs) {
        throw std::runtime_error(
            "economic dispatch: the case has no polynomial generator costs (mpc.gencost)");
    }
    (void)options;

    DcOpfResult result;

    std::vector<std::size_t> index;
    for (std::size_t g = 0; g < data.generators.size(); ++g) {
        if (data.generators[g].online) {
            index.push_back(g);
        }
    }
    result.generators = index;

    double demand = 0.0;
    for (const Bus& bus : data.buses) {
        demand += bus.pd;
    }
    result.total_load_mw = demand;

    // Lambda iteration: Pg = clamp((lambda - c1) / (2 c2)), bisecting on lambda
    // until the units meet the load.
    double low = -1e6;
    double high = 1e6;
    std::vector<double> pg(index.size(), 0.0);
    const auto dispatch = [&](double lambda) {
        double total = 0.0;
        for (std::size_t g = 0; g < index.size(); ++g) {
            const Generator& generator = data.generators[index[g]];
            const double pmax = generator.pmax > 0.0 ? generator.pmax : 0.0;
            double value = 0.0;
            if (generator.cost_c2 > 0.0) {
                value = (lambda - generator.cost_c1) / (2.0 * generator.cost_c2);
            } else if (lambda > generator.cost_c1) {
                value = pmax;
            }
            value = std::clamp(value, generator.pmin, pmax);
            pg[g] = value;
            total += value;
        }
        return total;
    };

    for (int iteration = 0; iteration < 200; ++iteration) {
        const double middle = 0.5 * (low + high);
        if (dispatch(middle) < demand) {
            low = middle;
        } else {
            high = middle;
        }
    }
    dispatch(0.5 * (low + high));

    result.pg_mw = pg;
    for (std::size_t g = 0; g < index.size(); ++g) {
        const Generator& generator = data.generators[index[g]];
        result.total_generation_mw += pg[g];
        result.objective += generator.cost_c2 * pg[g] * pg[g] + generator.cost_c1 * pg[g] +
                            generator.cost_c0;
    }
    result.lmp.assign(data.buses.size(), 0.5 * (low + high));
    result.theta_deg.assign(data.buses.size(), 0.0);
    result.branch_flow_mw.assign(data.branches.size(), 0.0);
    result.branch_congested.assign(data.branches.size(), 0);
    result.converged = std::abs(result.total_generation_mw - demand) < 1e-6;
    result.message = result.converged ? "converged" : "load not met";
    result.iterations = 200;
    return result;
}

std::string format_dc_opf_report(const Network& net, const DcOpfResult& result) {
    const CaseData& data = net.data();
    std::ostringstream out;

    out << "DC optimal power flow : " << data.name << '\n';
    out << "  status        : ";
    if (!result.converged) {
        out << "NOT CONVERGED (" << result.message << ")\n";
        return out.str();
    }
    if (result.iterations == 0) {
        out << "optimal, uncongested (the network-free dispatch already respects every "
               "branch rating)\n";
    } else {
        out << "converged in " << result.iterations << " interior-point iterations\n";
    }

    out << "  objective     : " << std::fixed << std::setprecision(2) << result.objective
        << " $/h\n";
    out << "  load          : " << std::setprecision(3) << result.total_load_mw << " MW\n";
    out << "  generation    : " << result.total_generation_mw << " MW\n";
    out << "  balance error : " << std::scientific << std::setprecision(2)
        << result.max_balance_error_mw << " MW\n";
    out << "  limit margin  : " << result.max_limit_violation_mw << " MW\n";
    out << "  barrier       : " << result.complementarity << '\n';

    double lowest = result.lmp.empty() ? 0.0 : result.lmp[0];
    double highest = lowest;
    for (const double price : result.lmp) {
        lowest = std::min(lowest, price);
        highest = std::max(highest, price);
    }
    out << "  LMP range     : " << std::fixed << std::setprecision(2) << lowest << " .. "
        << highest << " $/MWh\n";

    std::size_t congested = 0;
    for (const char flag : result.branch_congested) {
        congested += flag != 0 ? 1 : 0;
    }
    out << "  congested     : " << congested << " branch(es)\n\n";

    out << "dispatch\n";
    out << std::left << std::setw(6) << "bus" << std::right << std::setw(12) << "Pg (MW)"
        << std::setw(12) << "Pmin" << std::setw(12) << "Pmax" << std::setw(14) << "marginal cost"
        << '\n';
    out << std::string(56, '-') << '\n';
    const auto previous_flags = out.flags();
    for (std::size_t g = 0; g < result.generators.size(); ++g) {
        const Generator& generator = data.generators[result.generators[g]];
        const double marginal = generator.cost_c1 + 2.0 * generator.cost_c2 * result.pg_mw[g];
        out << std::left << std::setw(6) << generator.bus << std::right << std::setw(12)
            << std::fixed << std::setprecision(3) << result.pg_mw[g] << std::setw(12)
            << generator.pmin << std::setw(12) << (generator.pmax > 0.0 ? generator.pmax : 0.0)
            << std::setw(14) << std::setprecision(2) << marginal << '\n';
    }
    out.flags(previous_flags);
    return out.str();
}

}  // namespace pf
