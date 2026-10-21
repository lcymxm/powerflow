#include "powerflow/jacobian.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace pf {

RealMatrix build_power_flow_jacobian(const ComplexMatrix& ybus, const std::vector<double>& vm,
                                     const std::vector<double>& va_rad,
                                     const std::vector<Index>& non_slack,
                                     const std::vector<Index>& pq) {
    const std::size_t nb = ybus.rows();
    if (ybus.cols() != nb) {
        throw std::invalid_argument("build_power_flow_jacobian: Y-bus must be square");
    }
    if (vm.size() != nb || va_rad.size() != nb) {
        throw std::invalid_argument("build_power_flow_jacobian: voltage vectors must have one "
                                    "entry per bus");
    }
    if (non_slack.size() < pq.size()) {
        throw std::invalid_argument("build_power_flow_jacobian: pq must be a subset of non_slack");
    }

    // Calculated injections at the operating point: S = V .* conj(Y V).
    std::vector<Complex> voltage(nb, Complex{});
    for (Index i = 0; i < nb; ++i) {
        if (vm[i] == 0.0) {
            throw std::invalid_argument("build_power_flow_jacobian: bus " + std::to_string(i + 1) +
                                        " has zero voltage magnitude");
        }
        voltage[i] = std::polar(vm[i], va_rad[i]);
    }

    std::vector<Complex> current(nb, Complex{});
    for (Index r = 0; r < nb; ++r) {
        Complex acc{};
        for (Index c = 0; c < nb; ++c) {
            acc += ybus(r, c) * voltage[c];
        }
        current[r] = acc;
    }

    std::vector<Complex> s_calc(nb, Complex{});
    for (Index i = 0; i < nb; ++i) {
        s_calc[i] = voltage[i] * std::conj(current[i]);
    }

    const std::size_t n_theta = non_slack.size();
    const std::size_t n_vm = pq.size();
    RealMatrix jac(n_theta + n_vm, n_theta + n_vm, 0.0);

    for (std::size_t a = 0; a < n_theta; ++a) {
        const Index i = non_slack[a];
        const double vi = vm[i];
        const double p_i = s_calc[i].real();
        const double q_i = s_calc[i].imag();

        for (std::size_t b = 0; b < n_theta; ++b) {
            const Index j = non_slack[b];
            if (i == j) {
                jac(a, b) = -q_i - ybus(i, i).imag() * vi * vi;
            } else {
                const double theta = va_rad[i] - va_rad[j];
                jac(a, b) = vi * vm[j] * (ybus(i, j).real() * std::sin(theta) -
                                          ybus(i, j).imag() * std::cos(theta));
            }
        }

        for (std::size_t b = 0; b < n_vm; ++b) {
            const Index j = pq[b];
            if (i == j) {
                jac(a, n_theta + b) = p_i / vi + ybus(i, i).real() * vi;
            } else {
                const double theta = va_rad[i] - va_rad[j];
                jac(a, n_theta + b) = vi * (ybus(i, j).real() * std::cos(theta) +
                                            ybus(i, j).imag() * std::sin(theta));
            }
        }
    }

    for (std::size_t a = 0; a < n_vm; ++a) {
        const Index i = pq[a];
        const double vi = vm[i];
        const double p_i = s_calc[i].real();
        const double q_i = s_calc[i].imag();

        for (std::size_t b = 0; b < n_theta; ++b) {
            const Index j = non_slack[b];
            if (i == j) {
                jac(n_theta + a, b) = p_i - ybus(i, i).real() * vi * vi;
            } else {
                const double theta = va_rad[i] - va_rad[j];
                jac(n_theta + a, b) = -vi * vm[j] * (ybus(i, j).real() * std::cos(theta) +
                                                     ybus(i, j).imag() * std::sin(theta));
            }
        }

        for (std::size_t b = 0; b < n_vm; ++b) {
            const Index j = pq[b];
            if (i == j) {
                jac(n_theta + a, n_theta + b) = q_i / vi - ybus(i, i).imag() * vi;
            } else {
                const double theta = va_rad[i] - va_rad[j];
                jac(n_theta + a, n_theta + b) =
                    vi * (ybus(i, j).real() * std::sin(theta) - ybus(i, j).imag() * std::cos(theta));
            }
        }
    }

    return jac;
}

}  // namespace pf
