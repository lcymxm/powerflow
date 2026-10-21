#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "powerflow/case_file.hpp"
#include "powerflow/dense.hpp"
#include "powerflow/jacobian.hpp"
#include "powerflow/network.hpp"
#include "powerflow/sparse.hpp"
#include "powerflow/sparse_lu.hpp"
#include "test_framework.hpp"

namespace {

const std::string kDataDir = POWERFLOW_DATA_DIR;

std::string data_path(const std::string& name) { return kDataDir + "/" + name; }

// Sparsely banded, diagonally dominant matrix: easy to factorise but with real
// structure, so the ordering has something to work on.
pf::SparseMatrix<double> make_banded(std::size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    pf::SparseMatrix<double>::Builder builder(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        builder.add(i, i, 10.0 + std::abs(dist(rng)) * 5.0);
        if (i + 1 < n) {
            const double v = dist(rng) * 2.0;
            builder.add(i, i + 1, v);
            builder.add(i + 1, i, v);
        }
        if (i + 9 < n) {
            const double v = dist(rng) * 2.0;
            builder.add(i, i + 9, v);
            builder.add(i + 9, i, v);
        }
    }
    return builder.build();
}

double relative_residual(const pf::SparseMatrix<double>& a, const std::vector<double>& x,
                         const std::vector<double>& b) {
    const std::vector<double> ax = a.multiply(x);
    double norm_a = 0.0;
    for (std::size_t r = 0; r < a.rows(); ++r) {
        double row = 0.0;
        for (std::size_t k = a.row_begin(r); k < a.row_end(r); ++k) {
            row += std::abs(a.value(k));
        }
        norm_a = std::max(norm_a, row);
    }
    double norm_x = 0.0;
    double norm_b = 0.0;
    double residual = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        norm_x = std::max(norm_x, std::abs(x[i]));
        norm_b = std::max(norm_b, std::abs(b[i]));
        residual = std::max(residual, std::abs(ax[i] - b[i]));
    }
    return residual / (norm_a * norm_x + norm_b + 1e-300);
}

}  // namespace

TEST(Sparse, BuilderSumsDuplicatesAndSortsColumns) {
    pf::SparseMatrix<double>::Builder builder(3, 3);
    builder.add(0, 1, 1.0);
    builder.add(0, 1, 2.0);
    builder.add(0, 0, 5.0);
    builder.add(1, 0, 3.0);
    const pf::SparseMatrix<double> matrix = builder.build();

    CHECK_EQ(matrix.rows(), std::size_t{3});
    CHECK_EQ(matrix.nonzeros(), std::size_t{3});
    CHECK_NEAR(matrix.at(0, 1), 3.0, 1e-15);
    CHECK_NEAR(matrix.at(0, 0), 5.0, 1e-15);
    CHECK_EQ(matrix.column(matrix.row_begin(0)), std::size_t{0});
    CHECK_EQ(matrix.column(matrix.row_begin(0) + 1), std::size_t{1});
    CHECK_NEAR(matrix.at(2, 2), 0.0, 1e-15);
}

TEST(Sparse, RejectsOutOfRangeIndex) {
    pf::SparseMatrix<double>::Builder builder(2, 2);
    const auto add = [&builder]() { builder.add(0, 5, 1.0); };
    CHECK_THROWS_AS(add(), "out of range");
}

TEST(Sparse, MultiplyMatchesDenseMatrix) {
    const pf::SparseMatrix<double> a = make_banded(40, 3);
    const pf::DenseMatrix<double> dense = pf::to_dense(a);
    std::vector<double> x(40, 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] = 0.5 - 0.01 * static_cast<double>(i);
    }
    const std::vector<double> sparse_product = a.multiply(x);
    const std::vector<double> dense_product = pf::matvec(dense, x);
    for (std::size_t i = 0; i < x.size(); ++i) {
        CHECK_NEAR(sparse_product[i], dense_product[i], 1e-12);
    }
}

TEST(Sparse, SymmetricPatternSeesBothDirections) {
    pf::SparseMatrix<double>::Builder builder(3, 3);
    builder.add(0, 1, 1.0);  // only the upper triangle is stored
    builder.add(2, 2, 1.0);
    const auto adjacency = pf::symmetric_pattern(builder.build());
    CHECK_EQ(adjacency[0].size(), std::size_t{1});
    CHECK_EQ(adjacency[0][0], std::size_t{1});
    CHECK_EQ(adjacency[1].size(), std::size_t{1});
    CHECK_EQ(adjacency[1][0], std::size_t{0});
    CHECK(adjacency[2].empty());
}

TEST(MinimumDegree, ReturnsAPermutation) {
    const pf::SparseMatrix<double> a = make_banded(30, 11);
    const std::vector<std::size_t> order = pf::minimum_degree_order(pf::symmetric_pattern(a));
    CHECK_EQ(order.size(), std::size_t{30});
    std::vector<std::size_t> sorted = order;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 0; i < sorted.size(); ++i) {
        CHECK_EQ(sorted[i], i);
    }
}

TEST(MinimumDegree, ReducesTheFillIn) {
    const std::size_t n = 200;
    const pf::SparseMatrix<double> a = make_banded(n, 5);

    std::vector<double> truth(n, 1.0);
    const std::vector<double> rhs = a.multiply(truth);

    const pf::SparseLu plain(a, false);
    const pf::SparseLu ordered(a, true);

    CHECK(plain.stats().nnz_factors < n * n);
    CHECK(ordered.stats().nnz_factors < plain.stats().nnz_factors);
    CHECK_NEAR(relative_residual(a, plain.solve(rhs), rhs), 0.0, 1e-12);
    CHECK_NEAR(relative_residual(a, ordered.solve(rhs), rhs), 0.0, 1e-12);
}

TEST(SparseLu, SolvesDiagonallyDominantSystem) {
    const std::size_t n = 60;
    const pf::SparseMatrix<double> a = make_banded(n, 2);
    std::vector<double> truth(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        truth[i] = 0.25 * static_cast<double>(i) - 3.0;
    }
    const std::vector<double> rhs = a.multiply(truth);

    const pf::SparseLu solver(a, true);
    CHECK_FALSE(solver.is_singular());
    CHECK_NEAR(relative_residual(a, solver.solve(rhs), rhs), 0.0, 1e-12);
}

TEST(SparseLu, MatchesTheDenseSolver) {
    for (const std::size_t n : {7u, 33u, 120u}) {
        const pf::SparseMatrix<double> a = make_banded(n, 17u + static_cast<unsigned>(n));
        std::vector<double> truth(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            truth[i] = std::sin(0.1 * static_cast<double>(i)) + 0.3;
        }
        const std::vector<double> rhs = a.multiply(truth);

        const pf::SparseLu sparse(a, true);
        const std::vector<double> x_sparse = sparse.solve(rhs);

        const pf::DenseMatrix<double> dense = pf::to_dense(a);
        const pf::LuSolver<double> reference(dense);
        const std::vector<double> x_dense = reference.solve(rhs);

        for (std::size_t i = 0; i < n; ++i) {
            CHECK_NEAR(x_sparse[i], x_dense[i], 1e-10);
        }
    }
}

TEST(SparseLu, DetectsASingularMatrix) {
    pf::SparseMatrix<double>::Builder builder(3, 3);
    builder.add(0, 0, 1.0);
    builder.add(0, 1, 2.0);
    builder.add(1, 0, 2.0);
    builder.add(1, 1, 4.0);
    builder.add(2, 2, 1.0);
    const pf::SparseLu solver(builder.build(), true);
    CHECK(solver.is_singular());
    const std::vector<double> rhs{1.0, 2.0, 3.0};
    CHECK_THROWS_AS(solver.solve(rhs), "singular");
}

TEST(SparseJacobian, MatchesTheDenseAssembly) {
    for (const char* name : {"case14.m", "case57.m"}) {
        const pf::Network net(pf::load_case_file(data_path(name)));
        const std::size_t nb = net.bus_count();

        std::vector<double> vm(nb, 1.0);
        std::vector<double> va(nb, 0.0);
        for (std::size_t i = 0; i < nb; ++i) {
            vm[i] = net.bus(i).vm > 0.0 ? net.bus(i).vm : 1.0;
            va[i] = net.bus(i).va * pf::kDegToRad;
        }

        const pf::SparseMatrix<double> sparse =
            pf::assemble_sparse_jacobian(net.ybus(), vm, va, net.non_slack_buses(), net.pq_buses());
        const pf::DenseMatrix<double> dense =
            pf::build_power_flow_jacobian(net.ybus(), vm, va, net.non_slack_buses(),
                                          net.pq_buses());

        CHECK_EQ(sparse.rows(), dense.rows());
        for (std::size_t r = 0; r < dense.rows(); ++r) {
            for (std::size_t c = 0; c < dense.cols(); ++c) {
                CHECK_NEAR(sparse.at(r, c), dense(r, c), 1e-12);
            }
        }
    }
}

TEST(SparseJacobian, UsesFarLessStorageThanADenseOne) {
    const pf::Network net(pf::load_case_file(data_path("case118.m")));
    const std::size_t nb = net.bus_count();
    std::vector<double> vm(nb, 1.0);
    std::vector<double> va(nb, 0.0);
    for (std::size_t i = 0; i < nb; ++i) {
        vm[i] = net.bus(i).vm > 0.0 ? net.bus(i).vm : 1.0;
        va[i] = net.bus(i).va * pf::kDegToRad;
    }
    const pf::SparseMatrix<double> jac =
        pf::assemble_sparse_jacobian(net.ybus(), vm, va, net.non_slack_buses(), net.pq_buses());
    const double dense_entries = static_cast<double>(jac.rows()) * static_cast<double>(jac.cols());
    CHECK(static_cast<double>(jac.nonzeros()) < dense_entries / 8.0);
}
