// The analytic Jacobian of the polar power-flow equations, exposed as a public
// utility so it can be verified independently against a finite-difference
// approximation (see tests/test_power_flow.cpp).
#pragma once

#include <vector>

#include "powerflow/dense.hpp"
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
// Throws std::invalid_argument when the index lists are inconsistent, or when
// any voltage magnitude is zero (the polar form divides by it on the diagonal).
RealMatrix build_power_flow_jacobian(const ComplexMatrix& ybus, const std::vector<double>& vm,
                                     const std::vector<double>& va_rad,
                                     const std::vector<Index>& non_slack,
                                     const std::vector<Index>& pq);

}  // namespace pf
