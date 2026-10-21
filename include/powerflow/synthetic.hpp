// Synthetic test networks.
//
// The shipped IEEE cases stop at 118 buses, which is far too small to show what
// the sparse factorisation buys: at that size a dense LU is perfectly fine. This
// generator builds a meshed grid of arbitrary size with a realistic sparsity
// pattern (each bus touches at most four neighbours), so the scaling difference
// can be measured and reported.
#pragma once

#include <string>

#include "powerflow/types.hpp"

namespace pf {

// A roughly square grid of `bus_count` buses. Bus 1 is the slack, every fortieth
// bus carries a voltage-controlled generator, and every bus has a small load.
// Branches have low impedance so the network stays well conditioned, which is
// what a transmission grid looks like in per unit.
//
// Throws std::invalid_argument when bus_count is smaller than 4.
CaseData make_synthetic_case(int bus_count, unsigned seed = 1);

}  // namespace pf
