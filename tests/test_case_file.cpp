#include <string>
#include <vector>

#include "powerflow/case_file.hpp"
#include "test_framework.hpp"

namespace {

const std::string kDataDir = POWERFLOW_DATA_DIR;

std::string data_path(const std::string& name) { return kDataDir + "/" + name; }

// A minimal but structurally valid case. Used as the base for the error tests.
const char* kMinimalCase = R"(mpc.baseMVA = 100;
mpc.bus = [
    1 3 0   0 0 0 1 1.0  0  0 1 1.1 0.9;
    2 1 10  5 0 0 1 1.0  0  0 1 1.1 0.9;
];
mpc.gen = [
    1 10 0 100 -100 1.0 100 1 100 0;
];
mpc.branch = [
    1 2 0.01 0.1 0.02 0 0 0 0 0 1 -360 360;
];
)";

}  // namespace

TEST(CaseFile, LoadsIeee14FromDisk) {
    const pf::CaseData data = pf::load_case_file(data_path("case14.m"));
    CHECK_EQ(data.name, std::string("case14"));
    CHECK_EQ(data.buses.size(), std::size_t{14});
    CHECK_EQ(data.generators.size(), std::size_t{5});
    CHECK_EQ(data.branches.size(), std::size_t{20});
    CHECK_NEAR(data.base_mva, 100.0, 1e-12);
    CHECK(data.source_path.find("case14.m") != std::string::npos);
}

TEST(CaseFile, ReadsBusColumns) {
    const pf::CaseData data = pf::load_case_file(data_path("case14.m"));
    CHECK_EQ(data.buses[0].id, 1);
    CHECK(data.buses[0].type == pf::BusType::Slack);
    CHECK(data.buses[1].type == pf::BusType::PV);
    CHECK(data.buses[3].type == pf::BusType::PQ);
    CHECK_NEAR(data.buses[3].pd, 47.8, 1e-12);
    CHECK_NEAR(data.buses[3].qd, -3.9, 1e-12);
    CHECK_NEAR(data.buses[8].bs, 19.0, 1e-12);  // the only bus shunt in case14
    CHECK_NEAR(data.buses[13].vmin, 0.94, 1e-12);
}

TEST(CaseFile, ReadsBranchTapAndRating) {
    const pf::CaseData data = pf::load_case_file(data_path("case14.m"));
    CHECK_EQ(data.branches[7].from, 4);
    CHECK_EQ(data.branches[7].to, 7);
    CHECK_NEAR(data.branches[7].ratio, 0.978, 1e-12);
    CHECK_NEAR(data.branches[7].x, 0.20912, 1e-12);
    CHECK_NEAR(data.branches[7].rate_a, 0.0, 1e-12);
    CHECK(data.branches[7].in_service);
}

TEST(CaseFile, ReadsGeneratorLimits) {
    const pf::CaseData data = pf::load_case_file(data_path("case14.m"));
    CHECK_EQ(data.generators[0].bus, 1);
    CHECK_NEAR(data.generators[0].pg, 232.4, 1e-12);
    CHECK_NEAR(data.generators[0].qmax, 10.0, 1e-12);
    CHECK_NEAR(data.generators[0].vg, 1.06, 1e-12);
    CHECK(data.generators[0].online);
}

TEST(CaseFile, IgnoresCommentLines) {
    const pf::CaseData data = pf::parse_case_text(kMinimalCase, "minimal");
    CHECK_EQ(data.buses.size(), std::size_t{2});
    CHECK_EQ(data.branches.size(), std::size_t{1});
    CHECK_NEAR(data.branches[0].b, 0.02, 1e-12);
}

TEST(CaseFile, DoesNotConfuseBusNameWithBus) {
    const std::string text = R"(mpc.baseMVA = 100;
mpc.bus_name = { 'north'; 'south' };
mpc.bus = [
    1 3 0  0 0 0 1 1.0 0 0 1 1.1 0.9;
    2 1 10 5 0 0 1 1.0 0 0 1 1.1 0.9;
];
mpc.gen = [
    1 10 0 100 -100 1.0 100 1 100 0;
];
mpc.branch = [
    1 2 0.01 0.1 0 0 0 0 0 0 1 -360 360;
];
)";
    const pf::CaseData data = pf::parse_case_text(text, "names");
    CHECK_EQ(data.buses.size(), std::size_t{2});
}

TEST(CaseFile, RejectsMissingBranchSection) {
    const std::string text = R"(mpc.baseMVA = 100;
mpc.bus = [ 1 3 0 0 0 0 1 1.0 0 0 1 1.1 0.9; ];
mpc.gen = [ 1 0 0 10 -10 1.0 100 1 10 0 ];
)";
    CHECK_THROWS_AS(pf::parse_case_text(text, "broken"), "mpc.branch");
}

TEST(CaseFile, RejectsMissingBaseMva) {
    const std::string text = R"(mpc.bus = [ 1 3 0 0 0 0 1 1.0 0 0 1 1.1 0.9; ];
mpc.gen = [ 1 0 0 10 -10 1.0 100 1 10 0 ];
mpc.branch = [ 1 2 0.01 0.1 0 0 0 0 0 0 1 -360 360 ];
)";
    CHECK_THROWS_AS(pf::parse_case_text(text, "broken"), "mpc.baseMVA");
}

TEST(CaseFile, RejectsTooFewColumns) {
    const std::string text = R"(mpc.baseMVA = 100;
mpc.bus = [ 1 3 0 0 0 0 1 1.0 0 0; ];
mpc.gen = [ 1 0 0 10 -10 1.0 100 1 10 0 ];
mpc.branch = [ 1 2 0.01 0.1 0 0 0 0 0 0 1 -360 360 ];
)";
    CHECK_THROWS_AS(pf::parse_case_text(text, "broken"), "columns");
}

TEST(CaseFile, RejectsNonNumericValue) {
    const std::string text = R"(mpc.baseMVA = 100;
mpc.bus = [ 1 3 0 0 0 0 1 1.0 zero 0 1 1.1 0.9; ];
mpc.gen = [ 1 0 0 10 -10 1.0 100 1 10 0 ];
mpc.branch = [ 1 2 0.01 0.1 0 0 0 0 0 0 1 -360 360 ];
)";
    CHECK_THROWS_AS(pf::parse_case_text(text, "broken"), "non-numeric");
}

TEST(CaseFile, RejectsUnknownBusType) {
    const std::string text = R"(mpc.baseMVA = 100;
mpc.bus = [ 1 7 0 0 0 0 1 1.0 0 0 1 1.1 0.9; ];
mpc.gen = [ 1 0 0 10 -10 1.0 100 1 10 0 ];
mpc.branch = [ 1 2 0.01 0.1 0 0 0 0 0 0 1 -360 360 ];
)";
    CHECK_THROWS_AS(pf::parse_case_text(text, "broken"), "unsupported bus type");
}

TEST(CaseFile, RejectsUnreadablePath) {
    CHECK_THROWS_AS(pf::load_case_file("no/such/directory/case99.m"), "cannot open");
}

TEST(CaseFile, AcceptsCommaSeparatedRows) {
    const std::string text = R"(mpc.baseMVA = 100;
mpc.bus = [
    1, 3, 0, 0, 0, 0, 1, 1.0, 0, 0, 1, 1.1, 0.9;
    2, 1, 10, 5, 0, 0, 1, 1.0, 0, 0, 1, 1.1, 0.9;
];
mpc.gen = [ 1, 10, 0, 100, -100, 1.0, 100, 1, 100, 0 ];
mpc.branch = [ 1, 2, 0.01, 0.1, 0, 0, 0, 0, 0, 0, 1, -360, 360 ];
)";
    const pf::CaseData data = pf::parse_case_text(text, "commas");
    CHECK_EQ(data.buses.size(), std::size_t{2});
    CHECK_NEAR(data.buses[1].pd, 10.0, 1e-12);
}

TEST(CaseFile, ReadsEveryShippedCase) {
    const std::vector<std::string> names{"case9.m", "case14.m", "case30.m", "case57.m",
                                         "case118.m"};
    for (const std::string& name : names) {
        const pf::CaseData data = pf::load_case_file(data_path(name));
        CHECK(!data.buses.empty());
        CHECK(!data.generators.empty());
        CHECK(!data.branches.empty());
        CHECK_NEAR(data.base_mva, 100.0, 1e-12);
    }
}
