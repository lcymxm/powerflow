// Sparse LU factorisation with a fill-reducing ordering.
//
// The dense factorisation in dense.hpp is O(n^3) time and O(n^2) memory. Both
// are wasted on a power-flow Jacobian, which has only a handful of nonzeros per
// row: what dominates the cost is *fill-in*, the new nonzeros that elimination
// creates. Two ingredients keep it small:
//
//   * a minimum-degree ordering, which eliminates the least-connected buses
//     first (see minimum_degree_order), and
//   * threshold partial pivoting, which keeps the diagonal pivot whenever it is
//     at least `threshold` times the largest candidate in its column. That
//     preserves the ordering's fill reduction while still refusing tiny pivots.
#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "powerflow/sparse.hpp"

namespace pf {

// Greedy minimum-degree ordering of an undirected graph given as adjacency
// lists. Returns the elimination order: order[i] is the node eliminated i-th.
//
// This is the classical algorithm: repeatedly pick the surviving node with the
// fewest surviving neighbours, append it to the order, and connect all of its
// surviving neighbours to each other (the clique that elimination creates).
std::vector<std::size_t> minimum_degree_order(
    const std::vector<std::vector<std::size_t>>& adjacency);

struct SparseLuStats {
    std::size_t nnz_ordered = 0;  // nonzeros after the ordering, before elimination
    std::size_t nnz_factors = 0;  // nonzeros of L (unit diagonal excluded) plus U
    std::size_t row_swaps = 0;    // how often threshold pivoting moved a row
};

class SparseLu {
public:
    // `threshold` is the pivoting rule from UMFPACK: a diagonal pivot is kept
    // when |a_kk| >= threshold * max|a_ik|. Zero means "never move a row".
    explicit SparseLu(const SparseMatrix<double>& matrix, bool use_ordering = true,
                      double threshold = 0.1);

    // Reuses an ordering computed for an earlier matrix with the same structure.
    // The Newton loop assembles a new Jacobian every iteration but its pattern
    // never changes, so the ordering only has to be computed once.
    SparseLu(const SparseMatrix<double>& matrix, const std::vector<std::size_t>& ordering,
             double threshold = 0.1);

    bool is_singular() const noexcept { return singular_; }
    const SparseLuStats& stats() const noexcept { return stats_; }
    // order()[i] is the original index of the row/column eliminated i-th.
    const std::vector<std::size_t>& order() const noexcept { return order_; }

    std::vector<double> solve(const std::vector<double>& rhs) const;

private:
    void build(const SparseMatrix<double>& matrix, const std::vector<std::size_t>& ordering);
    void rebuild_column_index();
    void factorise();

    std::size_t n_ = 0;
    // Each row is a column-ascending list without duplicates. Elimination keeps
    // them sorted by merging, which avoids the per-entry allocation overhead that
    // a node-based container would pay on every fill-in.
    using Row = std::vector<std::pair<std::size_t, double>>;
    std::vector<Row> rows_;
    // rows touching each column, ascending. Without this every pivot would have
    // to scan all rows, which is O(n^2 log n) and dominates everything else.
    std::vector<std::vector<std::size_t>> column_rows_;
    std::vector<std::size_t> pivot_;
    std::vector<std::size_t> order_;
    bool singular_ = false;
    double threshold_ = 0.1;
    SparseLuStats stats_;
};

}  // namespace pf
