// Core value types shared by every module of the power-flow engine.
#pragma once

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace pf {

using Complex = std::complex<double>;
using Index = std::size_t;

inline constexpr double kDegToRad = 3.14159265358979323846 / 180.0;
inline constexpr double kRadToDeg = 180.0 / 3.14159265358979323846;
inline constexpr Index kInvalidIndex = static_cast<Index>(-1);

// MATPOWER bus type codes, kept numerically identical to the case format.
enum class BusType : int {
    PQ = 1,
    PV = 2,
    Slack = 3,
    Isolated = 4,
};

inline const char* to_string(BusType type) {
    switch (type) {
        case BusType::PQ:
            return "PQ";
        case BusType::PV:
            return "PV";
        case BusType::Slack:
            return "SL";
        case BusType::Isolated:
            return "ISO";
    }
    return "?";
}

// One row of `mpc.bus`. Powers are in MW / MVAr, voltages in pu and degrees.
struct Bus {
    int id = 0;
    BusType type = BusType::PQ;
    double pd = 0.0;
    double qd = 0.0;
    double gs = 0.0;
    double bs = 0.0;
    double vm = 1.0;
    double va = 0.0;
    double vmax = 1.1;
    double vmin = 0.9;
};

// One row of `mpc.gen`.
struct Generator {
    int bus = 0;
    double pg = 0.0;
    double qg = 0.0;
    double qmax = 0.0;
    double qmin = 0.0;
    double vg = 1.0;
    bool online = true;
    double pmax = 0.0;
    double pmin = 0.0;
};

// One row of `mpc.branch`. The tap is stored as magnitude plus phase shift in
// degrees; a magnitude of 0 means "no transformer, ratio 1".
struct Branch {
    int from = 0;
    int to = 0;
    double r = 0.0;
    double x = 0.0;
    double b = 0.0;
    double rate_a = 0.0;
    double ratio = 0.0;
    double shift = 0.0;
    bool in_service = true;
};

// Everything a power-flow study needs from a case file.
struct CaseData {
    std::string name;
    std::string source_path;
    double base_mva = 100.0;
    std::vector<Bus> buses;
    std::vector<Generator> generators;
    std::vector<Branch> branches;
};

}  // namespace pf
