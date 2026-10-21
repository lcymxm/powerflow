#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "powerflow/dense.hpp"
#include "test_framework.hpp"

namespace {

pf::RealMatrix make_matrix(const std::vector<std::vector<double>>& rows) {
    pf::RealMatrix m(rows.size(), rows.empty() ? 0 : rows.front().size(), 0.0);
    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (std::size_t c = 0; c < rows[r].size(); ++c) {
            m(r, c) = rows[r][c];
        }
    }
    return m;
}

}  // namespace

TEST(Dense, MatvecMultipliesRows) {
    const pf::RealMatrix m = make_matrix({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
    const std::vector<double> out = pf::matvec(m, std::vector<double>{1.0, 1.0, 1.0});
    CHECK_EQ(out.size(), std::size_t{2});
    CHECK_NEAR(out[0], 6.0, 1e-12);
    CHECK_NEAR(out[1], 15.0, 1e-12);
}

TEST(Dense, MatvecRejectsMismatchedLength) {
    const pf::RealMatrix m(2, 2, 1.0);
    const std::vector<double> wrong{1.0, 2.0, 3.0};
    CHECK_THROWS_AS(pf::matvec(m, wrong), "vector length");
}

TEST(Dense, IdentitySolvesToTheRightHandSide) {
    const pf::RealMatrix identity = pf::RealMatrix::identity(4);
    const pf::LuSolver<double> solver(identity);
    const std::vector<double> rhs{1.0, -2.0, 3.5, 0.0};
    const std::vector<double> x = solver.solve(rhs);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        CHECK_NEAR(x[i], rhs[i], 1e-12);
    }
}

TEST(Dense, SolvesWhenTheLeadingPivotIsZero) {
    // [[0, 2], [1, 3]] x = [4, 5]  =>  x = [-1, 2], only reachable after a row swap.
    const pf::RealMatrix a = make_matrix({{0.0, 2.0}, {1.0, 3.0}});
    const pf::LuSolver<double> solver(a);
    CHECK_FALSE(solver.is_singular());
    const std::vector<double> x = solver.solve(std::vector<double>{4.0, 5.0});
    CHECK_NEAR(x[0], -1.0, 1e-12);
    CHECK_NEAR(x[1], 2.0, 1e-12);
}

TEST(Dense, SolvesThreeByThreeToSmallResidual) {
    const pf::RealMatrix a = make_matrix({{4.0, 1.0, 0.0}, {1.0, 3.0, 1.0}, {0.0, 1.0, 2.0}});
    const std::vector<double> rhs{1.0, 2.0, 3.0};
    const pf::LuSolver<double> solver(a);
    const std::vector<double> x = solver.solve(rhs);
    const std::vector<double> check = pf::matvec(a, x);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        CHECK_NEAR(check[i], rhs[i], 1e-10);
    }
}

TEST(Dense, RefusesSingularMatrix) {
    // Second row is twice the first: rank one.
    const pf::RealMatrix a = make_matrix({{1.0, 2.0}, {2.0, 4.0}});
    const pf::LuSolver<double> solver(a);
    CHECK(solver.is_singular());
    const std::vector<double> rhs{1.0, 2.0};
    CHECK_THROWS_AS(solver.solve(rhs), "singular");
}

TEST(Dense, RejectsNonSquareMatrix) {
    const pf::RealMatrix a(2, 3, 1.0);
    const auto factorise = [&a]() {
        const pf::LuSolver<double> solver(a);
        CHECK_NEAR(solver.size(), 2.0, 1e-12);
    };
    CHECK_THROWS_AS(factorise(), "square");
}

TEST(Dense, SolvesPermutedSystemThatNeedsManyPivots) {
    // Regression test. A diagonally dominant matrix stays well conditioned under
    // any row permutation, but a permuted one forces several row interchanges
    // during the factorisation. Applying those interchanges column by column
    // instead of row by row still converges on small systems and then quietly
    // loses most of the accuracy here, so this case has to stay covered.
    const std::size_t n = 37;
    std::uint32_t seed = 12345u;
    const auto next = [&seed]() {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<double>(seed >> 8) / static_cast<double>(1u << 24);
    };

    pf::RealMatrix base(n, n, 0.0);
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t c = 0; c < n; ++c) {
            base(r, c) = (next() - 0.5) * 20.0;
        }
        base(r, r) += 200.0;
    }

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) {
        order[i] = i;
    }
    for (std::size_t i = n; i-- > 1;) {
        const std::size_t j = static_cast<std::size_t>(next() * static_cast<double>(i + 1));
        std::swap(order[i], order[j]);
    }

    pf::RealMatrix a(n, n, 0.0);
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t c = 0; c < n; ++c) {
            a(r, c) = base(order[r], c);
        }
    }

    std::vector<double> truth(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        truth[i] = next() - 0.5;
    }
    const std::vector<double> rhs = pf::matvec(a, truth);

    const pf::LuSolver<double> solver(a);
    CHECK_FALSE(solver.is_singular());
    const std::vector<double> x = solver.solve(rhs);
    const std::vector<double> check = pf::matvec(a, x);

    double norm_a = 0.0;
    double norm_x = 0.0;
    double norm_rhs = 0.0;
    double residual = 0.0;
    for (std::size_t r = 0; r < n; ++r) {
        double row_sum = 0.0;
        for (std::size_t c = 0; c < n; ++c) {
            row_sum += std::abs(a(r, c));
        }
        norm_a = std::max(norm_a, row_sum);
        residual = std::max(residual, std::abs(check[r] - rhs[r]));
    }
    for (std::size_t i = 0; i < n; ++i) {
        norm_x = std::max(norm_x, std::abs(x[i]));
        norm_rhs = std::max(norm_rhs, std::abs(rhs[i]));
    }
    CHECK(residual <= 1e-10 * (norm_a * norm_x + norm_rhs));
    for (std::size_t i = 0; i < n; ++i) {
        CHECK_NEAR(x[i], truth[i], 1e-9);
    }
}

TEST(Dense, SolvesComplexSystem) {
    const pf::DenseMatrix<pf::Complex> a = pf::DenseMatrix<pf::Complex>::identity(2);
    const pf::LuSolver<pf::Complex> solver(a);
    const std::vector<pf::Complex> rhs{pf::Complex{1.0, 2.0}, pf::Complex{-3.0, 0.5}};
    const std::vector<pf::Complex> x = solver.solve(rhs);
    CHECK_NEAR(std::abs(x[0] - rhs[0]), 0.0, 1e-12);
    CHECK_NEAR(std::abs(x[1] - rhs[1]), 0.0, 1e-12);
}
