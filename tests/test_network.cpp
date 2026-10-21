#include <cmath>
#include <complex>
#include <string>

#include "powerflow/case_file.hpp"
#include "powerflow/network.hpp"
#include "test_framework.hpp"

namespace {

const std::string kDataDir = POWERFLOW_DATA_DIR;

std::string data_path(const std::string& name) { return kDataDir + "/" + name; }

// Two buses, one line, one slack generator: the smallest solvable network.
pf::CaseData make_two_bus_case() {
    pf::CaseData data;
    data.name = "synthetic";
    data.base_mva = 100.0;

    pf::Bus slack;
    slack.id = 1;
    slack.type = pf::BusType::Slack;
    slack.vm = 1.0;

    pf::Bus load;
    load.id = 2;
    load.type = pf::BusType::PQ;
    load.pd = 50.0;
    load.qd = 20.0;
    load.vm = 1.0;

    data.buses = {slack, load};

    pf::Generator generator;
    generator.bus = 1;
    generator.pg = 50.0;
    generator.qmax = 100.0;
    generator.qmin = -100.0;
    generator.vg = 1.0;
    generator.online = true;
    data.generators = {generator};

    pf::Branch branch;
    branch.from = 1;
    branch.to = 2;
    branch.r = 0.01;
    branch.x = 0.05;
    branch.b = 0.01;
    branch.in_service = true;
    data.branches = {branch};

    return data;
}

// Wrapper so CHECK_THROWS_AS can call the constructor without hitting the most
// vexing parse.
pf::Network build(const pf::CaseData& data) { return pf::Network(data); }

}  // namespace

TEST(Network, ClassifiesIeee14) {
    const pf::Network net(pf::load_case_file(data_path("case14.m")));
    CHECK_EQ(net.bus_count(), std::size_t{14});
    CHECK_EQ(net.branch_count(), std::size_t{20});
    CHECK_EQ(net.slack_bus(), pf::Index{0});
    CHECK_EQ(net.pv_buses().size(), std::size_t{4});
    CHECK_EQ(net.pq_buses().size(), std::size_t{9});
    CHECK_EQ(net.non_slack_buses().size(), std::size_t{13});
    CHECK(net.isolated_buses().empty());
    CHECK(net.bus_type(pf::Index{0}) == pf::BusType::Slack);
    CHECK(net.bus_type(pf::Index{1}) == pf::BusType::PV);
    CHECK(net.bus_type(pf::Index{3}) == pf::BusType::PQ);
}

TEST(Network, YbusIsSymmetric) {
    const pf::Network net(pf::load_case_file(data_path("case14.m")));
    const pf::SparseMatrix<pf::Complex>& y = net.ybus();
    for (std::size_t i = 0; i < net.bus_count(); ++i) {
        for (std::size_t j = 0; j < net.bus_count(); ++j) {
            CHECK_NEAR(std::abs(y(i, j) - y(j, i)), 0.0, 1e-12);
        }
    }
}

TEST(Network, YbusIsSparse) {
    // A connected network cannot have more than a few nonzeros per row; the
    // whole point of the sparse storage is that this stays O(n) rather than n^2.
    const pf::Network net(pf::load_case_file(data_path("case118.m")));
    const std::size_t dense = net.bus_count() * net.bus_count();
    CHECK(net.ybus().nonzeros() < dense / 4);
    CHECK(net.ybus().nonzeros() >= 2 * net.branch_count());
}

TEST(Network, YbusRowSumIsTheBusCharging) {
    // No tap on any line at bus 1 of case14, so the off-diagonal terms cancel
    // the series part and each row sums to j * (half of the line charging).
    const pf::Network net(pf::load_case_file(data_path("case14.m")));
    pf::Complex expected{0.0, 0.0};
    for (std::size_t k = 0; k < net.branch_count(); ++k) {
        const pf::Branch& branch = net.branch(k);
        if (!net.branch_in_service(k) || branch.ratio != 0.0) {
            continue;
        }
        if (branch.from == 1 || branch.to == 1) {
            expected += pf::Complex{0.0, branch.b / 2.0};
        }
    }
    pf::Complex row_sum{0.0, 0.0};
    for (std::size_t j = 0; j < net.bus_count(); ++j) {
        row_sum += net.ybus()(0, j);
    }
    CHECK_NEAR(std::abs(row_sum - expected), 0.0, 1e-12);
}

TEST(Network, BuildsSeriesAdmittanceForASimpleLine) {
    const pf::CaseData data = make_two_bus_case();
    const pf::Network net(data);
    const pf::Complex series = 1.0 / pf::Complex{0.01, 0.05};
    const pf::Complex terminal = series + pf::Complex{0.0, 0.01 / 2.0};

    CHECK_NEAR(std::abs(net.ybus()(0, 0) - terminal), 0.0, 1e-12);
    CHECK_NEAR(std::abs(net.ybus()(1, 1) - terminal), 0.0, 1e-12);
    CHECK_NEAR(std::abs(net.ybus()(0, 1) + series), 0.0, 1e-12);
    CHECK_NEAR(std::abs(net.ybus()(1, 0) + series), 0.0, 1e-12);
}

TEST(Network, AppliesTransformerTapOnTheFromSide) {
    pf::CaseData data = make_two_bus_case();
    data.branches[0].ratio = 1.05;
    const pf::Network net(data);

    const pf::Complex series = 1.0 / pf::Complex{0.01, 0.05};
    const pf::Complex ytt = series + pf::Complex{0.0, 0.01 / 2.0};
    const pf::Complex yff = ytt / (1.05 * 1.05);
    const pf::Complex yft = -series / 1.05;

    CHECK_NEAR(std::abs(net.ybus()(0, 0) - yff), 0.0, 1e-12);
    CHECK_NEAR(std::abs(net.ybus()(1, 1) - ytt), 0.0, 1e-12);
    CHECK_NEAR(std::abs(net.ybus()(0, 1) - yft), 0.0, 1e-12);
}

TEST(Network, AddsBusShuntToTheDiagonal) {
    pf::CaseData data = make_two_bus_case();
    data.buses[1].bs = 30.0;  // 30 MVAr at 1 pu on a 100 MVA base
    const pf::Network net(data);

    const pf::Complex expected{0.0, 0.3};
    CHECK_NEAR(std::abs(net.bus_shunts()[1] - expected), 0.0, 1e-12);
    CHECK_NEAR(std::abs(net.ybus()(1, 0) - net.ybus()(0, 1)), 0.0, 1e-12);
}

TEST(Network, RejectsDuplicateBusIds) {
    pf::CaseData data = make_two_bus_case();
    data.buses[1].id = 1;
    CHECK_THROWS_AS(build(data), "duplicate bus id");
}

TEST(Network, RejectsUnknownBranchEndpoint) {
    pf::CaseData data = make_two_bus_case();
    data.branches[0].to = 7;
    CHECK_THROWS_AS(build(data), "unknown bus 7");
}

TEST(Network, RejectsSelfLoop) {
    pf::CaseData data = make_two_bus_case();
    data.branches[0].to = 1;
    CHECK_THROWS_AS(build(data), "connects a bus to itself");
}

TEST(Network, RejectsZeroImpedanceBranch) {
    pf::CaseData data = make_two_bus_case();
    data.branches[0].r = 0.0;
    data.branches[0].x = 0.0;
    CHECK_THROWS_AS(build(data), "zero series impedance");
}

TEST(Network, RejectsTwoSlackBuses) {
    pf::CaseData data = make_two_bus_case();
    data.buses[1].type = pf::BusType::Slack;
    CHECK_THROWS_AS(build(data), "exactly one slack bus");
}

TEST(Network, RejectsCaseWithoutSlack) {
    pf::CaseData data = make_two_bus_case();
    data.buses[0].type = pf::BusType::PQ;
    CHECK_THROWS_AS(build(data), "exactly one slack bus");
}

TEST(Network, MarksUnconnectedBusAsIsolated) {
    pf::CaseData data = make_two_bus_case();
    pf::Bus orphan;
    orphan.id = 3;
    orphan.type = pf::BusType::PQ;
    orphan.pd = 5.0;
    orphan.vm = 1.0;
    data.buses.push_back(orphan);

    const pf::Network net = build(data);
    CHECK_EQ(net.isolated_buses().size(), std::size_t{1});
    CHECK_EQ(net.isolated_buses()[0], pf::Index{2});
    CHECK_EQ(net.pq_buses().size(), std::size_t{1});
    CHECK_EQ(net.non_slack_buses().size(), std::size_t{1});
}

TEST(Network, IndexOfRejectsUnknownId) {
    const pf::Network net(pf::load_case_file(data_path("case14.m")));
    CHECK_EQ(net.index_of(14), pf::Index{13});
    CHECK(net.has_bus(7));
    CHECK_FALSE(net.has_bus(99));
    CHECK_THROWS_AS(net.index_of(99), "unknown bus id 99");
}

TEST(Network, TotalsLoadAndGeneration) {
    const pf::Network net(pf::load_case_file(data_path("case14.m")));
    CHECK_NEAR(net.total_load_mw(), 259.0, 1e-9);
    CHECK_NEAR(net.total_load_mvar(), 73.5, 1e-9);
    CHECK_NEAR(net.total_generation_mw(), 272.4, 1e-9);
    CHECK_NEAR(net.base_mva(), 100.0, 1e-12);
}

TEST(Network, IgnoresOutOfServiceBranch) {
    pf::CaseData data = make_two_bus_case();
    // Add a third bus joined by a branch that is switched off: the third bus
    // must end up isolated rather than breaking the solve.
    pf::Bus extra;
    extra.id = 3;
    extra.type = pf::BusType::PQ;
    extra.pd = 5.0;
    extra.vm = 1.0;
    data.buses.push_back(extra);

    pf::Branch switched_off;
    switched_off.from = 2;
    switched_off.to = 3;
    switched_off.r = 0.01;
    switched_off.x = 0.05;
    switched_off.in_service = false;
    data.branches.push_back(switched_off);

    const pf::Network net = build(data);
    CHECK_EQ(net.isolated_buses().size(), std::size_t{1});
    CHECK_FALSE(net.branch_in_service(1));
    CHECK_NEAR(std::abs(net.ybus()(2, 2)), 0.0, 1e-12);
}
