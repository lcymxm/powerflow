#include "powerflow/jacobian.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace pf {
namespace {

// Where each bus lands in the Jacobian's row and column numbering. kInvalidIndex
// means "not part of that block".
struct Layout {
    std::vector<Index> theta_row;  // dP row of a bus in non_slack
    std::vector<Index> q_row;      // dQ row of a bus in pq
    std::vector<Index> theta_col;  // dVa column of a bus in non_slack
    std::vector<Index> vm_col;     // dVm column of a bus in pq
};

Layout make_layout(std::size_t bus_count, const std::vector<Index>& non_slack,
                   const std::vector<Index>& pq) {
    Layout layout;
    layout.theta_row.assign(bus_count, kInvalidIndex);
    layout.q_row.assign(bus_count, kInvalidIndex);
    layout.theta_col.assign(bus_count, kInvalidIndex);
    layout.vm_col.assign(bus_count, kInvalidIndex);

    const std::size_t n_theta = non_slack.size();
    for (std::size_t k = 0; k < n_theta; ++k) {
        layout.theta_row[non_slack[k]] = k;
        layout.theta_col[non_slack[k]] = k;
    }
    for (std::size_t k = 0; k < pq.size(); ++k) {
        layout.q_row[pq[k]] = n_theta + k;
        layout.vm_col[pq[k]] = n_theta + k;
    }
    return layout;
}

}  // namespace

SparseMatrix<double> assemble_sparse_jacobian(const SparseMatrix<Complex>& ybus,
                                              const std::vector<double>& vm,
                                              const std::vector<double>& va_rad,
                                              const std::vector<Index>& non_slack,
                                              const std::vector<Index>& pq) {
    const std::size_t nb = ybus.rows();
    if (ybus.cols() != nb) {
        throw std::invalid_argument("assemble_sparse_jacobian: Y-bus must be square");
    }
    if (vm.size() != nb || va_rad.size() != nb) {
        throw std::invalid_argument("assemble_sparse_jacobian: voltage vectors must have one "
                                    "entry per bus");
    }
    if (non_slack.size() < pq.size()) {
        throw std::invalid_argument("assemble_sparse_jacobian: pq must be a subset of non_slack");
    }

    std::vector<Complex> voltage(nb, Complex{});
    for (Index i = 0; i < nb; ++i) {
        if (vm[i] == 0.0) {
            throw std::invalid_argument("assemble_sparse_jacobian: bus " + std::to_string(i + 1) +
                                        " has zero voltage magnitude");
        }
        voltage[i] = std::polar(vm[i], va_rad[i]);
    }

    const std::vector<Complex> current = matvec(ybus, voltage);
    std::vector<Complex> s_calc(nb, Complex{});
    for (Index i = 0; i < nb; ++i) {
        s_calc[i] = voltage[i] * std::conj(current[i]);
    }

    const Layout layout = make_layout(nb, non_slack, pq);
    const std::size_t n = non_slack.size() + pq.size();

    SparseMatrix<double>::Builder builder(n, n);
    builder.reserve(ybus.nonzeros() * 2 + n * 2);

    for (const Index i : non_slack) {
        const double vi = vm[i];
        const double p_i = s_calc[i].real();
        const double q_i = s_calc[i].imag();
        const Index row_theta = layout.theta_row[i];
        const Index row_q = layout.q_row[i];

        // Walking the structural nonzeros of row i of the Y-bus covers every
        // Jacobian entry that can be nonzero: an entry is proportional to G_ij
        // and B_ij off the diagonal, and the diagonal is always stored.
        for (std::size_t k = ybus.row_begin(i); k < ybus.row_end(i); ++k) {
            const Index j = ybus.column(k);
            const double g_ij = ybus.value(k).real();
            const double b_ij = ybus.value(k).imag();

            double h = 0.0;        // dP_i / dVa_j
            double n_entry = 0.0;  // dP_i / dVm_j
            double m_entry = 0.0;  // dQ_i / dVa_j
            double l = 0.0;        // dQ_i / dVm_j

            if (i == j) {
                h = -q_i - b_ij * vi * vi;
                n_entry = p_i / vi + g_ij * vi;
                m_entry = p_i - g_ij * vi * vi;
                l = q_i / vi - b_ij * vi;
            } else {
                const double theta = va_rad[i] - va_rad[j];
                const double sin_theta = std::sin(theta);
                const double cos_theta = std::cos(theta);
                h = vi * vm[j] * (g_ij * sin_theta - b_ij * cos_theta);
                n_entry = vi * (g_ij * cos_theta + b_ij * sin_theta);
                m_entry = -vi * vm[j] * (g_ij * cos_theta + b_ij * sin_theta);
                l = vi * (g_ij * sin_theta - b_ij * cos_theta);
            }

            const Index theta_col = layout.theta_col[j];
            if (theta_col != kInvalidIndex) {
                builder.add(row_theta, theta_col, h);
                if (row_q != kInvalidIndex) {
                    builder.add(row_q, theta_col, m_entry);
                }
            }
            const Index vm_col = layout.vm_col[j];
            if (vm_col != kInvalidIndex) {
                builder.add(row_theta, vm_col, n_entry);
                if (row_q != kInvalidIndex) {
                    builder.add(row_q, vm_col, l);
                }
            }
        }
    }

    return builder.build();
}

DenseMatrix<double> build_power_flow_jacobian(const SparseMatrix<Complex>& ybus,
                                              const std::vector<double>& vm,
                                              const std::vector<double>& va_rad,
                                              const std::vector<Index>& non_slack,
                                              const std::vector<Index>& pq) {
    return to_dense(assemble_sparse_jacobian(ybus, vm, va_rad, non_slack, pq));
}

}  // namespace pf
