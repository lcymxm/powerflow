#include "powerflow/network.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace pf {

Network::Network(CaseData data) : data_(std::move(data)), bus_count_(data_.buses.size()) {
    validate();
    classify();
    build_ybus();
}

void Network::validate() {
    if (bus_count_ == 0) {
        throw std::runtime_error(data_.name + ": the case contains no buses");
    }

    id_to_index_.reserve(bus_count_ * 2);
    for (Index i = 0; i < bus_count_; ++i) {
        const int id = data_.buses[i].id;
        const auto inserted = id_to_index_.emplace(id, i);
        if (!inserted.second) {
            throw std::runtime_error(data_.name + ": duplicate bus id " + std::to_string(id));
        }
    }

    for (Index k = 0; k < data_.branches.size(); ++k) {
        const Branch& br = data_.branches[k];
        const std::string where = data_.name + ": branch " + std::to_string(k + 1) + " (" +
                                  std::to_string(br.from) + "-" + std::to_string(br.to) + ") ";
        if (!has_bus(br.from)) {
            throw std::runtime_error(where + "references unknown bus " + std::to_string(br.from));
        }
        if (!has_bus(br.to)) {
            throw std::runtime_error(where + "references unknown bus " + std::to_string(br.to));
        }
        if (br.from == br.to) {
            throw std::runtime_error(where + "connects a bus to itself");
        }
        if (br.in_service && br.r == 0.0 && br.x == 0.0) {
            throw std::runtime_error(where + "has zero series impedance (r = x = 0), "
                                             "which this solver does not support");
        }
    }
}

void Network::classify() {
    std::vector<std::size_t> degree(bus_count_, 0);
    for (const Branch& br : data_.branches) {
        if (!br.in_service) {
            continue;
        }
        ++degree[index_of(br.from)];
        ++degree[index_of(br.to)];
    }

    bus_types_.assign(bus_count_, BusType::PQ);
    for (Index i = 0; i < bus_count_; ++i) {
        const BusType declared = data_.buses[i].type;
        // A bus that participates in no in-service branch is electrically
        // isolated: it has no solution to contribute and would make the
        // Jacobian singular.
        bus_types_[i] = (declared == BusType::Isolated || degree[i] == 0) ? BusType::Isolated
                                                                         : declared;
    }

    slack_ = kInvalidIndex;
    pv_.clear();
    pq_.clear();
    isolated_.clear();
    non_slack_.clear();

    std::size_t slack_count = 0;
    for (Index i = 0; i < bus_count_; ++i) {
        switch (bus_types_[i]) {
            case BusType::Slack:
                ++slack_count;
                slack_ = i;
                break;
            case BusType::PV:
                pv_.push_back(i);
                break;
            case BusType::PQ:
                pq_.push_back(i);
                break;
            case BusType::Isolated:
                isolated_.push_back(i);
                break;
        }
    }

    if (slack_count != 1) {
        throw std::runtime_error(data_.name + ": expected exactly one slack bus, found " +
                                 std::to_string(slack_count));
    }

    non_slack_.reserve(pv_.size() + pq_.size());
    non_slack_.insert(non_slack_.end(), pv_.begin(), pv_.end());
    non_slack_.insert(non_slack_.end(), pq_.begin(), pq_.end());
}

void Network::build_ybus() {
    ybus_ = ComplexMatrix(bus_count_, bus_count_, Complex{});
    bus_shunts_.assign(bus_count_, Complex{});
    branch_adm_.assign(data_.branches.size(), BranchAdmittance{});
    branch_on_.assign(data_.branches.size(), false);

    for (Index k = 0; k < data_.branches.size(); ++k) {
        const Branch& br = data_.branches[k];
        if (!br.in_service) {
            continue;
        }
        branch_on_[k] = true;

        const Complex series = Complex{1.0, 0.0} / Complex{br.r, br.x};
        const Complex half_charging{0.0, br.b / 2.0};
        const double tap_magnitude = br.ratio != 0.0 ? br.ratio : 1.0;
        const Complex tap = std::polar(tap_magnitude, br.shift * kDegToRad);

        BranchAdmittance& adm = branch_adm_[k];
        adm.shunt_from = half_charging;
        adm.shunt_to = half_charging;
        adm.ytt = series + half_charging;
        adm.yff = adm.ytt / (tap * std::conj(tap));
        adm.yft = -series / std::conj(tap);
        adm.ytf = -series / tap;

        const Index f = index_of(br.from);
        const Index t = index_of(br.to);
        ybus_(f, f) += adm.yff;
        ybus_(f, t) += adm.yft;
        ybus_(t, f) += adm.ytf;
        ybus_(t, t) += adm.ytt;
    }

    for (Index i = 0; i < bus_count_; ++i) {
        const Bus& b = data_.buses[i];
        bus_shunts_[i] = Complex{b.gs, b.bs} / data_.base_mva;
        ybus_(i, i) += bus_shunts_[i];
    }
}

Index Network::index_of(int bus_id) const {
    const auto it = id_to_index_.find(bus_id);
    if (it == id_to_index_.end()) {
        throw std::out_of_range("unknown bus id " + std::to_string(bus_id) + " in case " +
                                data_.name);
    }
    return it->second;
}

double Network::total_load_mw() const {
    double sum = 0.0;
    for (const Bus& b : data_.buses) {
        sum += b.pd;
    }
    return sum;
}

double Network::total_load_mvar() const {
    double sum = 0.0;
    for (const Bus& b : data_.buses) {
        sum += b.qd;
    }
    return sum;
}

double Network::total_generation_mw() const {
    double sum = 0.0;
    for (const Generator& g : data_.generators) {
        if (g.online) {
            sum += g.pg;
        }
    }
    return sum;
}

}  // namespace pf
