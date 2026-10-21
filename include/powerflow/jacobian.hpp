// Assembly of the polar power-flow Jacobian.
//
// Two entry points share one assembly routine:
//   * assemble_sparse_jacobian() is what the Newton solver consumes, and
//   * build_power_flow_jacobian() returns the same matrix densely, which is the
//     form the finite-difference cross-check in the tests wants.
#pragma once

#include <vector>

#include "powerflow/dense.hpp"
#include "powerflow/sparse.hpp"
#include "powerflow/types.hpp"

namespace pf {

// Builds
//
//     dP/dVa  dP/dVm        rows: non_slack (for dP), then pq (for dQ)
//     dQ/dVa  dQ/dVm        cols: non_slack (for dVa), then pq (for dVm)
//
// evaluated at the operating point given by `vm` (pu) and `va_rad` (radians).
// `non_slack` must list every PV and PQ bus, `pq` every PQ bus. Both lists are
// used in the order given, which fixes the row and column order of the result.
//
// Only entries the Y-bus structure can make nonzero are assembled, so the result
// needs O(nnz(Y-bus)) storage rather than O(n^2).
//
// Throws std::invalid_argument when the index lists are inconsistent, or when
// any voltage magnitude is zero (the polar form divides by it on the diagonal).
SparseMatrix<double> assemble_sparse_jacobian(const SparseMatrix<Complex>& ybus,
                                              const std::vector<double>& vm,
                                              const std::vector<double>& va_rad,
                                              const std::vector<Index>& non_slack,
                                              const std::vector<Index>& pq);

// The same matrix, dense. Used by the tests and by the ``--linear dense`` option.
DenseMatrix<double> build_power_flow_jacobian(const SparseMatrix<Complex>& ybus,
                                              const std::vector<double>& vm,
                                              const std::vector<double>& va_rad,
                                              const std::vector<Index>& non_slack,
                                              const std::vector<Index>& pq);

}  // namespace pf
