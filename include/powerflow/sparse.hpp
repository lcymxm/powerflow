// Compressed sparse row matrix and the helpers the sparse solver needs.
//
// The power-flow Jacobian and the bus admittance matrix are both extremely
// sparse (a few nonzeros per bus), so storing them densely throws away almost
// all of the memory and most of the arithmetic. This type keeps only the
// structural nonzeros.
#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "powerflow/dense.hpp"
#include "powerflow/types.hpp"

namespace pf {

template <typename Scalar>
class SparseMatrix {
public:
    SparseMatrix() = default;

    SparseMatrix(std::size_t rows, std::size_t cols)
        : rows_(rows), cols_(cols), row_start_(rows + 1, 0) {}

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    std::size_t nonzeros() const noexcept { return columns_.size(); }

    // Value at (r, c); zero when the entry was never assembled. The lookup is a
    // binary search inside the row, so it is fine for inspection but not for
    // inner loops -- those should walk row_begin()/row_end().
    Scalar at(std::size_t r, std::size_t c) const {
        const std::size_t begin = row_start_[r];
        const std::size_t end = row_start_[r + 1];
        const auto* first = columns_.data() + begin;
        const auto* last = columns_.data() + end;
        const auto* found = std::lower_bound(first, last, c);
        if (found == last || *found != c) {
            return Scalar{};
        }
        return values_[static_cast<std::size_t>(found - columns_.data())];
    }

    Scalar operator()(std::size_t r, std::size_t c) const { return at(r, c); }

    std::size_t row_begin(std::size_t r) const noexcept { return row_start_[r]; }
    std::size_t row_end(std::size_t r) const noexcept { return row_start_[r + 1]; }
    std::size_t column(std::size_t k) const noexcept { return columns_[k]; }
    const Scalar& value(std::size_t k) const noexcept { return values_[k]; }

    std::vector<Scalar> multiply(const std::vector<Scalar>& x) const {
        if (x.size() != cols_) {
            throw std::invalid_argument("SparseMatrix::multiply: dimension mismatch");
        }
        std::vector<Scalar> out(rows_, Scalar{});
        for (std::size_t r = 0; r < rows_; ++r) {
            Scalar acc{};
            for (std::size_t k = row_start_[r]; k < row_start_[r + 1]; ++k) {
                acc += values_[k] * x[columns_[k]];
            }
            out[r] = acc;
        }
        return out;
    }

    // Triplet assembly: duplicates are summed and columns end up sorted.
    class Builder {
    public:
        Builder(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols) {}

        void add(std::size_t r, std::size_t c, Scalar value) {
            if (r >= rows_ || c >= cols_) {
                throw std::out_of_range("SparseMatrix::Builder::add: index out of range");
            }
            rows_index_.push_back(r);
            columns_index_.push_back(c);
            values_.push_back(value);
        }

        void reserve(std::size_t count) {
            rows_index_.reserve(count);
            columns_index_.reserve(count);
            values_.reserve(count);
        }

        std::size_t entry_count() const noexcept { return values_.size(); }

        SparseMatrix build() const {
            SparseMatrix matrix(rows_, cols_);
            std::vector<std::size_t> order(values_.size());
            std::iota(order.begin(), order.end(), std::size_t{0});
            std::sort(order.begin(), order.end(), [this](std::size_t a, std::size_t b) {
                if (rows_index_[a] != rows_index_[b]) {
                    return rows_index_[a] < rows_index_[b];
                }
                return columns_index_[a] < columns_index_[b];
            });

            std::vector<std::size_t> merged_rows;
            std::vector<std::size_t> merged_cols;
            std::vector<Scalar> merged_values;
            merged_rows.reserve(order.size());
            merged_cols.reserve(order.size());
            merged_values.reserve(order.size());

            for (const std::size_t k : order) {
                const std::size_t r = rows_index_[k];
                const std::size_t c = columns_index_[k];
                if (!merged_rows.empty() && merged_rows.back() == r && merged_cols.back() == c) {
                    merged_values.back() += values_[k];
                } else {
                    merged_rows.push_back(r);
                    merged_cols.push_back(c);
                    merged_values.push_back(values_[k]);
                }
            }

            matrix.row_start_.assign(rows_ + 1, 0);
            for (const std::size_t r : merged_rows) {
                matrix.row_start_[r + 1] += 1;
            }
            for (std::size_t r = 0; r < rows_; ++r) {
                matrix.row_start_[r + 1] += matrix.row_start_[r];
            }
            matrix.columns_ = std::move(merged_cols);
            matrix.values_ = std::move(merged_values);
            return matrix;
        }

    private:
        std::size_t rows_ = 0;
        std::size_t cols_ = 0;
        std::vector<std::size_t> rows_index_;
        std::vector<std::size_t> columns_index_;
        std::vector<Scalar> values_;
    };

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<std::size_t> row_start_;
    std::vector<std::size_t> columns_;
    std::vector<Scalar> values_;
};

template <typename Scalar>
std::vector<Scalar> matvec(const SparseMatrix<Scalar>& matrix, const std::vector<Scalar>& x) {
    return matrix.multiply(x);
}

// Sparsity pattern of A + A^T without the diagonal, i.e. the undirected graph a
// symmetric ordering should be computed on.
template <typename Scalar>
std::vector<std::vector<std::size_t>> symmetric_pattern(const SparseMatrix<Scalar>& matrix) {
    std::vector<std::vector<std::size_t>> adjacency(matrix.rows());
    for (std::size_t r = 0; r < matrix.rows(); ++r) {
        for (std::size_t k = matrix.row_begin(r); k < matrix.row_end(r); ++k) {
            const std::size_t c = matrix.column(k);
            if (c == r) {
                continue;
            }
            adjacency[r].push_back(c);
            adjacency[c].push_back(r);
        }
    }
    for (auto& neighbours : adjacency) {
        std::sort(neighbours.begin(), neighbours.end());
        neighbours.erase(std::unique(neighbours.begin(), neighbours.end()), neighbours.end());
    }
    return adjacency;
}

template <typename Scalar>
DenseMatrix<Scalar> to_dense(const SparseMatrix<Scalar>& matrix) {
    DenseMatrix<Scalar> out(matrix.rows(), matrix.cols(), Scalar{});
    for (std::size_t r = 0; r < matrix.rows(); ++r) {
        for (std::size_t k = matrix.row_begin(r); k < matrix.row_end(r); ++k) {
            out(r, matrix.column(k)) = matrix.value(k);
        }
    }
    return out;
}

}  // namespace pf
