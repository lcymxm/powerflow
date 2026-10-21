// DC optimal power flow.
//
// Minimises the total generating cost subject to a lossless DC network model:
//
//   minimise    sum_g  c2_g Pg^2 + c1_g Pg + c0_g
//   subject to   Bbus * theta - Cg * Pg = -Pd          (nodal power balance)
//                theta_ref = 0                          (angle reference)
//                Pmin <= Pg <= Pmax                     (generator limits)
//                |Bf * theta| <= rateA                  (branch flow limits)
//
// That is a convex quadratic program, which this file solves with a primal-dual
// interior-point method (Mehrotra predictor-corrector). The dual variables of the
// nodal balance rows are the locational marginal prices.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "powerflow/network.hpp"

namespace pf {

struct DcOpfOptions {
    double tolerance = 1e-9;  // on the complementarity and the residuals
    int max_iterations = 100;
};

struct DcOpfResult {
    bool converged = false;
    std::string message;
    int iterations = 0;

    double objective = 0.0;                // $/h
    double complementarity = 0.0;          // final barrier measure
    double max_balance_error_mw = 0.0;     // max |Ax - b|
    double max_limit_violation_mw = 0.0;   // max violation of an inequality

    std::vector<double> theta_deg;         // one per bus
    std::vector<std::size_t> generators;   // indices into the case's generators
    std::vector<double> pg_mw;             // one per entry of `generators`
    std::vector<double> lmp;               // one per bus, $/MWh
    std::vector<double> branch_flow_mw;    // one per branch, from -> to
    std::vector<char> branch_congested;    // flow sits on its rating

    double total_load_mw = 0.0;
    double total_generation_mw = 0.0;
};

// Solves the DC OPF for `net`. Requires the case to carry generation costs
// (see CaseData::has_generation_costs); throws otherwise.
DcOpfResult solve_dc_opf(const Network& net, const DcOpfOptions& options = {});

// Active power that an economic dispatch would assign, ignoring the network.
// Useful as a sanity bound for the OPF objective.
DcOpfResult solve_economic_dispatch(const Network& net, const DcOpfOptions& options = {});

std::string format_dc_opf_report(const Network& net, const DcOpfResult& result);

}  // namespace pf
