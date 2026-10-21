#include "powerflow/synthetic.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace pf {
namespace {

// Small deterministic generator so a given seed always produces the same grid.
class Lcg {
public:
    explicit Lcg(std::uint32_t seed) : state_(seed) {}

    double next() {
        state_ = state_ * 1664525u + 1013904223u;
        return static_cast<double>(state_ >> 8) / static_cast<double>(1u << 24);
    }

    double between(double low, double high) { return low + (high - low) * next(); }

private:
    std::uint32_t state_;
};

}  // namespace

CaseData make_synthetic_case(int bus_count, unsigned seed) {
    if (bus_count < 4) {
        throw std::invalid_argument("make_synthetic_case: need at least 4 buses");
    }

    CaseData data;
    data.name = "synthetic-" + std::to_string(bus_count);
    data.base_mva = 100.0;

    const int columns = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(bus_count))));
    const int rows = (bus_count + columns - 1) / columns;

    Lcg random(seed);

    data.buses.reserve(static_cast<std::size_t>(bus_count));
    for (int index = 0; index < bus_count; ++index) {
        Bus bus;
        bus.id = index + 1;
        bus.type = BusType::PQ;
        bus.vm = 1.0;
        bus.va = 0.0;
        bus.vmax = 1.1;
        bus.vmin = 0.9;
        // A small, roughly uniform load keeps the far corner of the grid from
        // sagging into the nose of the P-V curve.
        bus.pd = random.between(6.0, 14.0);
        bus.qd = random.between(1.5, 3.5);
        data.buses.push_back(bus);
    }

    data.buses[0].type = BusType::Slack;
    data.buses[0].pd = 0.0;
    data.buses[0].qd = 0.0;

    // Generation has to roughly match the load, otherwise the slack would have
    // to carry a few gigawatts across the grid and the case would sit past the
    // nose of the P-V curve. Every `kGeneratorEvery`-th bus is a PV bus.
    constexpr int kGeneratorEvery = 25;
    double total_load_mw = 0.0;
    for (const Bus& bus : data.buses) {
        total_load_mw += bus.pd;
    }
    int generator_buses = 0;
    for (int index = 1; index < bus_count; ++index) {
        if (index % kGeneratorEvery == 0) {
            ++generator_buses;
        }
    }
    if (generator_buses == 0) {
        generator_buses = 1;
    }
    const double per_generator_mw = total_load_mw / static_cast<double>(generator_buses);

    for (int index = 0; index < bus_count; ++index) {
        const bool is_slack = index == 0;
        const bool is_generator = !is_slack && (index % kGeneratorEvery == 0);
        if (is_generator) {
            data.buses[static_cast<std::size_t>(index)].type = BusType::PV;
        }
        if (!is_slack && !is_generator) {
            continue;  // no generator at a plain load bus
        }

        Generator generator;
        generator.bus = index + 1;
        // Everyone covers an equal share of the load; the slack then absorbs the
        // losses found by the solve.
        generator.pg = per_generator_mw;
        generator.qmax = 80.0;
        generator.qmin = -80.0;
        generator.vg = is_slack ? 1.02 : 1.01;
        generator.online = true;
        generator.pmax = per_generator_mw * 2.0;
        generator.pmin = 0.0;
        data.buses[static_cast<std::size_t>(index)].vm = generator.vg;
        data.generators.push_back(generator);
    }

    const auto impedance = [&random]() {
        return std::pair<double, double>{random.between(0.004, 0.010), random.between(0.02, 0.05)};
    };

    data.branches.reserve(static_cast<std::size_t>(bus_count) * 2);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int index = row * columns + column;
            if (index >= bus_count) {
                continue;
            }

            // Right neighbour.
            const int right = index + 1;
            if (column + 1 < columns && right < bus_count) {
                const auto [r, x] = impedance();
                Branch branch;
                branch.from = index + 1;
                branch.to = right + 1;
                branch.r = r;
                branch.x = x;
                branch.in_service = true;
                data.branches.push_back(branch);
            }

            // Down neighbour.
            const int down = index + columns;
            if (row + 1 < rows && down < bus_count) {
                const auto [r, x] = impedance();
                Branch branch;
                branch.from = index + 1;
                branch.to = down + 1;
                branch.r = r;
                branch.x = x;
                branch.in_service = true;
                data.branches.push_back(branch);
            }
        }
    }

    return data;
}

}  // namespace pf
