// Small dense linear-algebra helpers: a matrix container and a LU solver with
// partial pivoting. Everything is header-only and templated on the scalar type
// so that the same code serves the real-valued Newton Jacobian and, if ever
// needed, complex-valued systems.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "powerflow/types.hpp"

namespace pf {

template <typename Scalar>
class DenseMatrix {
public:
    DenseMatrix() = default;

    DenseMatrix(std::size_t rows, std::size_t cols, Scalar value = Scalar{})
        : rows_(rows), cols_(cols), data_(rows * cols, value) {}

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    bool empty() const noexcept { return data_.empty(); }

    Scalar& operator()(std::size_t r, std::size_t c) noexcept { return data_[r * cols_ + c]; }
    const Scalar& operator()(std::size_t r, std::size_t c) const noexcept {
        return data_[r * cols_ + c];
    }

    const std::vector<Scalar>& storage() const noexcept { return data_; }
    std::vector<Scalar>& storage() noexcept { return data_; }

    static DenseMatrix identity(std::size_t n) {
        DenseMatrix m(n, n, Scalar{});
        for (std::size_t i = 0; i < n; ++i) {
            m(i, i) = Scalar{1};
        }
        return m;
    }

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<Scalar> data_;
};

using RealMatrix = DenseMatrix<double>;
using ComplexMatrix = DenseMatrix<Complex>;

template <typename Scalar>
std::vector<Scalar> matvec(const DenseMatrix<Scalar>& m, const std::vector<Scalar>& v) {
    if (v.size() != m.cols()) {
        throw std::invalid_argument("matvec: vector length does not match matrix columns");
    }
    std::vector<Scalar> out(m.rows(), Scalar{});
    for (std::size_t r = 0; r < m.rows(); ++r) {
        Scalar acc{};
        for (std::size_t c = 0; c < m.cols(); ++c) {
            acc += m(r, c) * v[c];
        }
        out[r] = acc;
    }
    return out;
}

// Dense LU factorisation with partial (row) pivoting.
//
// The factorisation is computed once and can then be reused for several
// right-hand sides, which is what the Newton loop wants: one factorisation per
// iteration, one solve per iteration.
template <typename Scalar>
class LuSolver {
public:
    explicit LuSolver(DenseMatrix<Scalar> a)
        : n_(a.rows()), lu_(std::move(a)), pivot_(n_), singular_(false), scale_(0.0) {
        if (lu_.rows() != lu_.cols()) {
            throw std::invalid_argument("LuSolver: matrix must be square");
        }
        for (const Scalar& value : lu_.storage()) {
            scale_ = std::max(scale_, std::abs(value));
        }
        factorise();
    }

    std::size_t size() const noexcept { return n_; }
    bool is_singular() const noexcept { return singular_; }

    // Solves `A x = b`. Throws when the factorisation hit a (numerically) zero
    // pivot, which for a power-flow Jacobian almost always means the network
    // splits into islands or a bus has no connection at all.
    std::vector<Scalar> solve(const std::vector<Scalar>& b) const {
        if (b.size() != n_) {
            throw std::invalid_argument("LuSolver::solve: right-hand side has wrong length");
        }
        if (singular_) {
            throw std::runtime_error("LuSolver::solve: matrix is singular");
        }

        std::vector<Scalar> x = b;

        // Forward substitution with the row interchanges applied on the fly.
        //
        // The order matters: row `i` must be swapped into place *before* it is
        // reduced, because the swap moves a right-hand-side entry that earlier
        // rows have not touched yet. Applying all the swaps up front and then
        // reducing column by column looks equivalent but silently uses the wrong
        // entry whenever two or more interchanges occur.
        for (std::size_t i = 0; i < n_; ++i) {
            if (pivot_[i] != i) {
                std::swap(x[i], x[pivot_[i]]);
            }
            Scalar acc = x[i];
            for (std::size_t j = 0; j < i; ++j) {
                acc -= lu_(i, j) * x[j];
            }
            x[i] = acc;
        }

        for (std::size_t k = n_; k > 0; --k) {
            const std::size_t row = k - 1;
            Scalar acc = x[row];
            for (std::size_t j = row + 1; j < n_; ++j) {
                acc -= lu_(row, j) * x[j];
            }
            x[row] = acc / lu_(row, row);
        }
        return x;
    }

private:
    void factorise() {
        if (n_ == 0) {
            return;
        }
        const double threshold = 1e-14 * (scale_ > 0.0 ? scale_ : 1.0);

        for (std::size_t k = 0; k < n_; ++k) {
            std::size_t best_row = k;
            double best_value = std::abs(lu_(k, k));
            for (std::size_t i = k + 1; i < n_; ++i) {
                const double candidate = std::abs(lu_(i, k));
                if (candidate > best_value) {
                    best_value = candidate;
                    best_row = i;
                }
            }
            pivot_[k] = best_row;

            if (best_value <= threshold) {
                singular_ = true;
                // Keep going so that `pivot_` stays valid; the solve path will
                // refuse to run. Avoid dividing by a zero pivot.
                if (std::abs(lu_(k, k)) == 0.0) {
                    continue;
                }
            }

            if (best_row != k) {
                for (std::size_t j = 0; j < n_; ++j) {
                    std::swap(lu_(k, j), lu_(best_row, j));
                }
            }

            const Scalar diagonal = lu_(k, k);
            for (std::size_t i = k + 1; i < n_; ++i) {
                lu_(i, k) /= diagonal;
                const Scalar factor = lu_(i, k);
                if (factor == Scalar{}) {
                    continue;
                }
                for (std::size_t j = k + 1; j < n_; ++j) {
                    lu_(i, j) -= factor * lu_(k, j);
                }
            }
        }
    }

    std::size_t n_ = 0;
    DenseMatrix<Scalar> lu_;
    std::vector<std::size_t> pivot_;
    bool singular_ = false;
    double scale_ = 0.0;
};

}  // namespace pf
