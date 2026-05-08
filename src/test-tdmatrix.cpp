#include "algorithms.h"
#include "tdmatrix.h"

#include <testthat.h>

#include <cmath>
#include <vector>

using namespace sc;

namespace
{
    bool approxEqual(double a, double b, double tol = 1e-9) noexcept
    {
        return std::abs(a - b) <= tol;
    }
}

context("TDMatrix")
{
    test_that("default-constructed matrix has zero size")
    {
        TDMatrix m;
        expect_true(m.size() == 0);
        expect_false(m.isPrepared());
    }

    test_that("sized constructor allocates the bands")
    {
        TDMatrix m(5);
        expect_true(m.size() == 5);
        expect_true(static_cast<int>(m.fullDiag().size()) == 5);
        expect_true(static_cast<int>(m.fullLower().size()) == 4);
        expect_true(static_cast<int>(m.fullUpper().size()) == 4);
    }

    test_that("absMax picks the largest absolute element across all bands")
    {
        TDMatrix m(3);
        m.diag(0)  = 1.0;
        m.diag(1)  = -7.0;
        m.diag(2)  = 2.0;
        m.lower(0) = 3.0;
        m.lower(1) = 4.0;
        m.upper(0) = -5.0;
        m.upper(1) = 6.0;
        expect_true(approxEqual(m.absMax(), 7.0));
    }

    test_that("multiplyBy scales every band")
    {
        TDMatrix m(3);
        m.diag(0)  = 1.0;
        m.diag(1)  = 2.0;
        m.lower(0) = 3.0;
        m.upper(0) = 4.0;
        m.multiplyBy(2.0);
        expect_true(approxEqual(m.diag(0), 2.0));
        expect_true(approxEqual(m.diag(1), 4.0));
        expect_true(approxEqual(m.lower(0), 6.0));
        expect_true(approxEqual(m.upper(0), 8.0));
    }

    test_that("matrix-vector product is correct on a small case")
    {
        // | 1 2 0 |   | 1 |   | 5  |
        // | 3 4 5 | * | 2 | = | 26 |
        // | 0 6 7 |   | 3 |   | 33 |
        TDMatrix m(3);
        m.diag(0)  = 1.0;  m.upper(0) = 2.0;
        m.lower(0) = 3.0;  m.diag(1)  = 4.0;  m.upper(1) = 5.0;
        m.lower(1) = 6.0;  m.diag(2)  = 7.0;
        std::vector<double> x{1.0, 2.0, 3.0};
        const auto y = m * x;
        expect_true(approxEqual(y[0], 5.0));
        expect_true(approxEqual(y[1], 26.0));
        expect_true(approxEqual(y[2], 33.0));
    }

    test_that("inlineMultiply matches operator*")
    {
        TDMatrix m(3);
        m.diag(0)  = 1.0;  m.upper(0) = 2.0;
        m.lower(0) = 3.0;  m.diag(1)  = 4.0;  m.upper(1) = 5.0;
        m.lower(1) = 6.0;  m.diag(2)  = 7.0;
        std::vector<double> x{1.0, 2.0, 3.0};
        const auto y = m * x;
        m.inlineMultiply(x);
        expect_true(approxEqual(x[0], y[0]));
        expect_true(approxEqual(x[1], y[1]));
        expect_true(approxEqual(x[2], y[2]));
    }

    test_that("isDiagonalDominant returns true for SPD-shaped matrices")
    {
        TDMatrix m(3);
        m.diag(0)  = 5.0;  m.diag(1)  = 5.0;  m.diag(2)  = 5.0;
        m.lower(0) = 1.0;  m.lower(1) = 1.0;
        m.upper(0) = 1.0;  m.upper(1) = 1.0;
        expect_true(m.isDiagonalDominant());
    }
}

context("Thomas solver")
{
    // Reference test: A * x = b, solve gives back x.
    auto build = [] {
        TDMatrix m(5);
        for (int i = 0; i < 5; ++i) m.diag(i) = static_cast<double>(i + 1);
        for (int i = 0; i < 4; ++i)
        {
            m.lower(i) = static_cast<double>(i + 2);
            m.upper(i) = static_cast<double>(i + 2);
        }
        return m;
    };
    // Expected x = (1, 2, 3, 4, 5), so b = A*x = (5, 15, 31, 53, 45).

    test_that("thomasIP recovers the known solution")
    {
        auto m = build();
        std::vector<double> b{5, 15, 31, 53, 45};
        algorithm::thomasIP(m, b);
        for (int i = 0; i < 5; ++i)
        {
            expect_true(approxEqual(b[static_cast<std::size_t>(i)],
                                    static_cast<double>(i + 1)));
        }
    }

    test_that("thomasReUseIP recovers the solution and is reusable")
    {
        auto m = build();

        std::vector<double> b1{5, 15, 31, 53, 45};
        algorithm::thomasReUseIP(m, b1);
        for (int i = 0; i < 5; ++i)
        {
            expect_true(approxEqual(b1[static_cast<std::size_t>(i)],
                                    static_cast<double>(i + 1)));
        }
        expect_true(m.isPrepared());

        // A second solve with the same matrix and a fresh rhs.
        std::vector<double> b2{5, 15, 31, 53, 45};
        algorithm::thomasReUseIP(m, b2);
        for (int i = 0; i < 5; ++i)
        {
            expect_true(approxEqual(b2[static_cast<std::size_t>(i)],
                                    static_cast<double>(i + 1)));
        }
    }
}
