#include "powerflow/sparse_lu.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pf {
namespace {

// Position of `column` in a sorted row, or row.size() when it is absent.
std::size_t find_column(const std::vector<std::pair<std::size_t, double>>& row,
                        std::size_t column) {
    const auto it = std::lower_bound(row.begin(), row.end(), column,
                                     [](const std::pair<std::size_t, double>& entry, std::size_t c) {
                                         return entry.first < c;
                                     });
    if (it == row.end() || it->first != column) {
        return row.size();
    }
    return static_cast<std::size_t>(it - row.begin());
}

}  // namespace

namespace {

// Inserts `value` into a sorted vector unless it is already there.
void insert_sorted(std::vector<std::size_t>& values, std::size_t value) {
    const auto position = std::lower_bound(values.begin(), values.end(), value);
    if (position == values.end() || *position != value) {
        values.insert(position, value);
    }
}

}  // namespace

std::vector<std::size_t> minimum_degree_order(
    const std::vector<std::vector<std::size_t>>& adjacency) {
    const std::size_t n = adjacency.size();
    if (n == 0) {
        return {};
    }

    // The elimination graph. Neighbour lists stay sorted so membership tests and
    // insertions are cheap and the memory does not fragment.
    std::vector<std::vector<std::size_t>> graph(n);
    std::size_t max_degree = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (const std::size_t j : adjacency[i]) {
            if (j != i) {
                insert_sorted(graph[i], j);
            }
        }
        max_degree = std::max(max_degree, graph[i].size());
    }

    std::vector<char> alive(n, 1);
    std::vector<std::size_t> degree(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        degree[i] = graph[i].size();
    }

    // Bucket queue keyed by degree. Entries are validated on pop, which lets a
    // degree that later shrinks simply be re-inserted in the lower bucket.
    std::vector<std::vector<std::size_t>> buckets(max_degree + 1);
    for (std::size_t i = 0; i < n; ++i) {
        buckets[degree[i]].push_back(i);
    }
    std::size_t smallest = 0;

    std::vector<std::size_t> order;
    order.reserve(n);
    std::vector<std::size_t> surviving;

    for (std::size_t step = 0; step < n; ++step) {
        std::size_t pick = n;
        while (smallest < buckets.size()) {
            while (!buckets[smallest].empty()) {
                const std::size_t candidate = buckets[smallest].back();
                buckets[smallest].pop_back();
                if (alive[candidate] != 0 && degree[candidate] == smallest) {
                    pick = candidate;
                    break;
                }
            }
            if (pick != n) {
                break;
            }
            ++smallest;
        }
        if (pick == n) {
            for (std::size_t i = 0; i < n; ++i) {
                if (alive[i] != 0) {
                    pick = i;
                    break;
                }
            }
            if (pick == n) {
                break;  // should not happen
            }
        }

        order.push_back(pick);
        alive[pick] = 0;

        surviving.clear();
        for (const std::size_t neighbour : graph[pick]) {
            if (alive[neighbour] != 0) {
                surviving.push_back(neighbour);
            }
        }

        // Eliminating `pick` turns its surviving neighbourhood into a clique.
        for (const std::size_t a : surviving) {
            for (const std::size_t b : surviving) {
                if (a != b) {
                    insert_sorted(graph[a], b);
                }
            }
        }

        degree[pick] = 0;
        for (const std::size_t node : surviving) {
            std::size_t alive_neighbours = 0;
            for (const std::size_t neighbour : graph[node]) {
                if (alive[neighbour] != 0) {
                    ++alive_neighbours;
                }
            }
            if (alive_neighbours >= buckets.size()) {
                buckets.resize(alive_neighbours + 1);
            }
            degree[node] = alive_neighbours;
            buckets[alive_neighbours].push_back(node);
            if (alive_neighbours < smallest) {
                smallest = alive_neighbours;
            }
        }
    }
    return order;
}

SparseLu::SparseLu(const SparseMatrix<double>& matrix, bool use_ordering, double threshold)
    : threshold_(threshold) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("SparseLu: matrix must be square");
    }
    n_ = matrix.rows();
    if (n_ == 0) {
        return;
    }

    std::vector<std::size_t> order;
    if (use_ordering) {
        order = minimum_degree_order(symmetric_pattern(matrix));
    } else {
        order.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            order[i] = i;
        }
    }
    build(matrix, order);
}

SparseLu::SparseLu(const SparseMatrix<double>& matrix, const std::vector<std::size_t>& ordering,
                   double threshold)
    : threshold_(threshold) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument("SparseLu: matrix must be square");
    }
    n_ = matrix.rows();
    if (n_ == 0) {
        return;
    }
    if (ordering.size() != n_) {
        throw std::invalid_argument("SparseLu: ordering length does not match the matrix");
    }
    build(matrix, ordering);
}

void SparseLu::build(const SparseMatrix<double>& matrix, const std::vector<std::size_t>& ordering) {
    order_ = ordering;
    order_.resize(n_);

    // position[original] = index in the permuted system
    std::vector<std::size_t> position(n_, 0);
    for (std::size_t i = 0; i < n_; ++i) {
        position[order_[i]] = i;
    }

    rows_.assign(n_, {});
    for (std::size_t r = 0; r < n_; ++r) {
        for (std::size_t k = matrix.row_begin(r); k < matrix.row_end(r); ++k) {
            const double value = matrix.value(k);
            if (value == 0.0) {
                continue;
            }
            rows_[position[r]].emplace_back(position[matrix.column(k)], value);
        }
    }
    for (Row& row : rows_) {
        std::sort(row.begin(), row.end(),
                  [](const std::pair<std::size_t, double>& a,
                     const std::pair<std::size_t, double>& b) { return a.first < b.first; });
    }

    stats_.nnz_ordered = 0;
    for (const Row& row : rows_) {
        stats_.nnz_ordered += row.size();
    }

    rebuild_column_index();
    factorise();
}

void SparseLu::rebuild_column_index() {
    column_rows_.assign(n_, {});
    for (std::size_t r = 0; r < n_; ++r) {
        for (const auto& [column, value] : rows_[r]) {
            (void)value;
            column_rows_[column].push_back(r);
        }
    }
}

void SparseLu::factorise() {
    pivot_.assign(n_, 0);

    for (std::size_t k = 0; k < n_; ++k) {
        // Largest candidate in column k, considering only the rows that actually
        // hold a nonzero there.
        std::size_t best_row = k;
        double best_value = 0.0;
        bool have_candidate = false;
        for (const std::size_t i : column_rows_[k]) {
            if (i < k) {
                continue;
            }
            const std::size_t position = find_column(rows_[i], k);
            if (position == rows_[i].size()) {
                continue;  // stale entry left behind by an earlier update
            }
            const double magnitude = std::abs(rows_[i][position].second);
            if (magnitude <= 0.0) {
                continue;
            }
            if (!have_candidate || magnitude > best_value) {
                best_value = magnitude;
                best_row = i;
                have_candidate = true;
            }
        }
        if (!have_candidate) {
            singular_ = true;
            return;
        }

        // Threshold pivoting: keep the diagonal when it is good enough, so the
        // fill-reducing ordering survives.
        if (best_row != k) {
            const std::size_t position = find_column(rows_[k], k);
            if (position != rows_[k].size() &&
                std::abs(rows_[k][position].second) >= threshold_ * best_value) {
                best_row = k;
            }
        }

        pivot_[k] = best_row;
        if (best_row != k) {
            std::swap(rows_[k], rows_[best_row]);
            ++stats_.row_swaps;
            // A swap moves entries between rows, so the column index has to be
            // rebuilt. Threshold pivoting makes this rare.
            rebuild_column_index();
        }

        const std::size_t pivot_position = find_column(rows_[k], k);
        if (pivot_position == rows_[k].size() || rows_[k][pivot_position].second == 0.0) {
            singular_ = true;
            return;
        }
        const double pivot_value = rows_[k][pivot_position].second;

        // Everything right of the pivot in the pivot row, for the update below.
        const Row& pivot_row = rows_[k];
        const std::size_t pivot_tail = pivot_position + 1;

        // Safe to iterate directly: the merge only ever appends to columns
        // strictly greater than k.
        for (const std::size_t i : column_rows_[k]) {
            if (i <= k) {
                continue;
            }
            const std::size_t entry = find_column(rows_[i], k);
            if (entry == rows_[i].size()) {
                continue;
            }
            const double multiplier = rows_[i][entry].second / pivot_value;
            rows_[i][entry].second = multiplier;  // L stays in place

            // Merge row i with -(multiplier) * (tail of row k).
            Row& target = rows_[i];
            Row merged;
            merged.reserve(target.size() + (pivot_row.size() - pivot_tail) + 1);

            std::size_t a = 0;
            std::size_t b = pivot_tail;
            while (a < target.size() || b < pivot_row.size()) {
                if (b >= pivot_row.size() ||
                    (a < target.size() && target[a].first < pivot_row[b].first)) {
                    merged.push_back(target[a]);
                    ++a;
                } else if (a >= target.size() || pivot_row[b].first < target[a].first) {
                    const std::size_t column = pivot_row[b].first;
                    merged.emplace_back(column, -multiplier * pivot_row[b].second);
                    // Fill-in: record that row i now touches this column.
                    column_rows_[column].push_back(i);
                    ++b;
                } else {
                    merged.emplace_back(target[a].first,
                                        target[a].second - multiplier * pivot_row[b].second);
                    ++a;
                    ++b;
                }
            }
            target = std::move(merged);
        }
    }

    stats_.nnz_factors = 0;
    for (const Row& row : rows_) {
        stats_.nnz_factors += row.size();
    }
}

std::vector<double> SparseLu::solve(const std::vector<double>& rhs) const {
    if (rhs.size() != n_) {
        throw std::invalid_argument("SparseLu::solve: right-hand side has wrong length");
    }
    if (singular_) {
        throw std::runtime_error("SparseLu::solve: matrix is singular");
    }
    if (n_ == 0) {
        return {};
    }

    std::vector<double> x(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        x[i] = rhs[order_[i]];
    }

    // Forward substitution: row i is swapped into place, then reduced with the
    // multipliers stored to the left of its diagonal.
    for (std::size_t i = 0; i < n_; ++i) {
        if (pivot_[i] != i) {
            std::swap(x[i], x[pivot_[i]]);
        }
        double acc = x[i];
        for (const auto& [column, value] : rows_[i]) {
            if (column >= i) {
                break;
            }
            acc -= value * x[column];
        }
        x[i] = acc;
    }

    // Back substitution.
    for (std::size_t k = n_; k > 0; --k) {
        const std::size_t row = k - 1;
        double acc = x[row];
        for (std::size_t p = rows_[row].size(); p > 0; --p) {
            const auto& [column, value] = rows_[row][p - 1];
            if (column <= row) {
                break;
            }
            acc -= value * x[column];
        }
        const std::size_t diagonal = find_column(rows_[row], row);
        if (diagonal == rows_[row].size()) {
            throw std::runtime_error("SparseLu::solve: missing diagonal entry");
        }
        x[row] = acc / rows_[row][diagonal].second;
    }

    std::vector<double> solution(n_);
    for (std::size_t i = 0; i < n_; ++i) {
        solution[order_[i]] = x[i];
    }
    return solution;
}

}  // namespace pf
