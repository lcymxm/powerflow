// Network model: resolves bus ids, classifies buses and assembles the complex
// bus admittance matrix (Y-bus) together with the per-branch 2x2 admittance
// blocks needed for the branch-flow post-processing.
#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "powerflow/dense.hpp"
#include "powerflow/types.hpp"

namespace pf {

// Terminal admittances of one branch, using the MATPOWER branch model:
// a series impedance, half of the total line-charging susceptance at each end,
// and an ideal complex tap ratio on the "from" side.
struct BranchAdmittance {
    Complex yff{};
    Complex yft{};
    Complex ytf{};
    Complex ytt{};
    Complex shunt_from{};
    Complex shunt_to{};
};

class Network {
public:
    // Takes ownership of a copy of the case, so the Network never dangles even
    // when it is constructed from a temporary.
    //
    // Throws std::runtime_error when the case is not solvable: duplicate bus
    // ids, a branch pointing at a missing bus, zero series impedance, or a
    // number of slack buses other than one.
    explicit Network(CaseData data);

    const CaseData& data() const noexcept { return data_; }
    double base_mva() const noexcept { return data_.base_mva; }
    std::size_t bus_count() const noexcept { return bus_count_; }
    std::size_t branch_count() const noexcept { return data_.branches.size(); }

    bool has_bus(int bus_id) const noexcept { return id_to_index_.count(bus_id) != 0; }
    // Zero-based index of a one-based case bus id; throws std::out_of_range.
    Index index_of(int bus_id) const;

    const Bus& bus(Index i) const { return data_.buses[i]; }
    const Branch& branch(Index k) const { return data_.branches[k]; }
    BusType bus_type(Index i) const { return bus_types_[i]; }

    Index slack_bus() const noexcept { return slack_; }
    const std::vector<Index>& pv_buses() const noexcept { return pv_; }
    const std::vector<Index>& pq_buses() const noexcept { return pq_; }
    const std::vector<Index>& isolated_buses() const noexcept { return isolated_; }
    // PV buses followed by PQ buses, i.e. every bus whose voltage angle is unknown.
    const std::vector<Index>& non_slack_buses() const noexcept { return non_slack_; }

    const ComplexMatrix& ybus() const noexcept { return ybus_; }
    // Bus shunt admittance in pu, i.e. (Gs + j*Bs) / baseMVA.
    const std::vector<Complex>& bus_shunts() const noexcept { return bus_shunts_; }
    const std::vector<BranchAdmittance>& branch_admittances() const noexcept { return branch_adm_; }
    bool branch_in_service(Index k) const noexcept { return branch_on_[k]; }

    double total_load_mw() const;
    double total_load_mvar() const;
    double total_generation_mw() const;

private:
    void validate();
    void classify();
    void build_ybus();

    CaseData data_;
    std::size_t bus_count_ = 0;
    std::unordered_map<int, Index> id_to_index_;
    std::vector<BusType> bus_types_;
    Index slack_ = kInvalidIndex;
    std::vector<Index> pv_;
    std::vector<Index> pq_;
    std::vector<Index> isolated_;
    std::vector<Index> non_slack_;
    ComplexMatrix ybus_;
    std::vector<Complex> bus_shunts_;
    std::vector<BranchAdmittance> branch_adm_;
    std::vector<bool> branch_on_;
};

}  // namespace pf
