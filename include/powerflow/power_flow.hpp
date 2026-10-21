// AC power-flow solvers. Newton-Raphson (polar form, analytic Jacobian) is the
// main solver; Gauss-Seidel is kept alongside it as a reference implementation
// that converges through a completely different iteration scheme.
#pragma once

#include <string>
#include <vector>

#include "powerflow/dense.hpp"
#include "powerflow/network.hpp"
#include "powerflow/types.hpp"

namespace pf {

enum class SolverKind {
    NewtonRaphson,
    GaussSeidel,
};

const char* to_string(SolverKind kind);

struct SolverOptions {
    double tolerance = 1e-8;      // pu, applied to the largest power mismatch
    int max_iterations = 30;      // Newton needs 4-6 for the IEEE cases
    bool enforce_q_limits = false;
    int max_q_limit_switches = 12;
};

struct SolverStatus {
    bool converged = false;
    int iterations = 0;
    double final_mismatch = 0.0;
    std::string message;
};

struct PowerFlowResult {
    SolverKind method = SolverKind::NewtonRaphson;
    SolverStatus status;
    double base_mva = 100.0;

    std::vector<double> vm;       // pu, one entry per bus
    std::vector<double> va;       // degrees
    std::vector<bool> connected;  // false for electrically isolated buses
    std::vector<Complex> s_inj;   // pu net complex injection per bus
    std::vector<Complex> s_gen;   // pu generation per bus
    std::vector<Complex> s_load;  // pu load per bus

    std::vector<Complex> sf;      // pu flow leaving the "from" terminal
    std::vector<Complex> st;      // pu flow leaving the "to" terminal
    std::vector<double> loading;  // |S| / rateA, zero for unrated branches

    std::vector<double> mismatch_history;
    double p_loss_mw = 0.0;
    double q_loss_mvar = 0.0;
    int q_limit_switches = 0;

    double total_generation_mw() const;
    double total_generation_mvar() const;
    double min_voltage() const;
    double max_voltage() const;
    // Index of the most heavily loaded branch, or kInvalidIndex when no branch
    // in the case carries a rating.
    Index worst_loaded_branch() const;
};

PowerFlowResult solve_newton_raphson(const Network& net, const SolverOptions& options = {});
PowerFlowResult solve_gauss_seidel(const Network& net, const SolverOptions& options = {});

}  // namespace pf
