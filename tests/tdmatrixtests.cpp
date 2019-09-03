#include "../../tdmatrix.h"
#include "catch.hpp"

using namespace sc;

TEST_CASE("Init test 1", "[TDMatrix]")
{
    // init test
    TDMatrix matrix(3);
    REQUIRE(matrix(0, 0) == 0.0);
    REQUIRE(matrix(1, 1) == 0.0);
    REQUIRE(matrix(2, 2) == 0.0);
    REQUIRE(matrix(0, 1) == 0.0);
    REQUIRE(matrix(1, 2) == 0.0);
    REQUIRE(matrix(1, 0) == 0.0);
    REQUIRE(matrix(2, 1) == 0.0);
}

TEST_CASE("Init test 2", "[TDMatrix]")
{
    // init test
    TDMatrix matrix(3);
    REQUIRE(matrix.diag(0) == 0.0);
    REQUIRE(matrix.diag(1) == 0.0);
    REQUIRE(matrix.diag(2) == 0.0);
    REQUIRE(matrix.lower(0) == 0.0);
    REQUIRE(matrix.lower(1) == 0.0);
    REQUIRE(matrix.upper(0) == 0.0);
    REQUIRE(matrix.upper(1) == 0.0);
}

TEST_CASE("Size test", "[TDMatrix]")
{
    TDMatrix matrix(3);
    REQUIRE(matrix.size() == 3);

    TDMatrix matrix2;
    REQUIRE(matrix2.size() == 0);
}

TEST_CASE("Assign test 1", "[TDMatrix]")
{
    // init test
    TDMatrix matrix(3);
    matrix(0, 0) = 1.0;
    matrix(1, 1) = 2.0;
    matrix(2, 2) = 3.0;
    matrix(0, 1) = 4.0;
    matrix(1, 2) = 5.0;
    matrix(1, 0) = 6.0;
    matrix(2, 1) = 7.0;

    REQUIRE(matrix(0, 0) == 1.0);
    REQUIRE(matrix(1, 1) == 2.0);
    REQUIRE(matrix(2, 2) == 3.0);
    REQUIRE(matrix(0, 1) == 4.0);
    REQUIRE(matrix(1, 2) == 5.0);
    REQUIRE(matrix(1, 0) == 6.0);
    REQUIRE(matrix(2, 1) == 7.0);
}

TEST_CASE("Assign test 2", "[TDMatrix]")
{
    // init test
    TDMatrix matrix(3);
    matrix.diag(0)  = 1.0;
    matrix.diag(1)  = 2.0;
    matrix.diag(2)  = 3.0;
    matrix.lower(0) = 4.0;
    matrix.lower(1) = 5.0;
    matrix.upper(0) = 6.0;
    matrix.upper(1) = 7.0;

    REQUIRE(matrix.diag(0) == 1.0);
    REQUIRE(matrix.diag(1) == 2.0);
    REQUIRE(matrix.diag(2) == 3.0);
    REQUIRE(matrix.lower(0) == 4.0);
    REQUIRE(matrix.lower(1) == 5.0);
    REQUIRE(matrix.upper(0) == 6.0);
    REQUIRE(matrix.upper(1) == 7.0);
}

TEST_CASE("Clear test", "[TDMatrix]")
{
    TDMatrix matrix(3);
    matrix(0, 0) = 1.0;
    matrix(1, 1) = 2.0;
    matrix(2, 2) = 3.0;
    matrix(0, 1) = 4.0;
    matrix(1, 2) = 5.0;
    matrix(1, 0) = 6.0;
    matrix(2, 1) = 7.0;
    matrix.clear();
    REQUIRE(matrix.diag(0) == 0.0);
    REQUIRE(matrix.diag(1) == 0.0);
    REQUIRE(matrix.diag(2) == 0.0);
    REQUIRE(matrix.lower(0) == 0.0);
    REQUIRE(matrix.lower(1) == 0.0);
    REQUIRE(matrix.upper(0) == 0.0);
    REQUIRE(matrix.upper(1) == 0.0);
}

TEST_CASE("Copy test", "[TDMatrix]")
{
    TDMatrix matrix(3);
    matrix.diag(0)  = 1.0;
    matrix.diag(1)  = 2.0;
    matrix.diag(2)  = 3.0;
    matrix.lower(0) = 4.0;
    matrix.lower(1) = 5.0;
    matrix.upper(0) = 6.0;
    matrix.upper(1) = 7.0;

    auto matrix2 = matrix;
    REQUIRE(matrix.diag(0) == 1.0);
    REQUIRE(matrix.diag(1) == 2.0);
    REQUIRE(matrix.diag(2) == 3.0);
    REQUIRE(matrix.lower(0) == 4.0);
    REQUIRE(matrix.lower(1) == 5.0);
    REQUIRE(matrix.upper(0) == 6.0);
    REQUIRE(matrix.upper(1) == 7.0);
}

TEST_CASE("Multi test", "[TDMatrix]")
{
    TDMatrix matrix(3);
    matrix.diag(0)  = 1.0;
    matrix.diag(1)  = 2.0;
    matrix.diag(2)  = 3.0;
    matrix.lower(0) = 4.0;
    matrix.lower(1) = 5.0;
    matrix.upper(0) = 6.0;
    matrix.upper(1) = 7.0;
    matrix.multiplyBy(2.0);

    REQUIRE(matrix.diag(0) == 2.0);
    REQUIRE(matrix.diag(1) == 4.0);
    REQUIRE(matrix.diag(2) == 6.0);
    REQUIRE(matrix.lower(0) == 8.0);
    REQUIRE(matrix.lower(1) == 10.0);
    REQUIRE(matrix.upper(0) == 12.0);
    REQUIRE(matrix.upper(1) == 14.0);
}

TEST_CASE("Max test 1", "[TDMatrix]")
{
    TDMatrix matrix(3);
    matrix.diag(0)  = 1.0;
    matrix.diag(1)  = 2.0;
    matrix.diag(2)  = 3.0;
    matrix.lower(0) = 4.0;
    matrix.lower(1) = 5.0;
    matrix.upper(0) = 6.0;
    matrix.upper(1) = 7.0;

    REQUIRE(matrix.max() == 7.0);
    matrix.lower(0) = 19.0;
    REQUIRE(matrix.max() == 19.0);
    matrix.lower(1) = 125.0;
    REQUIRE(matrix.max() == 125.0);
}

TEST_CASE("M-V-Multi test", "[TDMatrix]")
{
    TDMatrix matrix(3);
    matrix.diag(0)  = 1.0;
    matrix.diag(1)  = 2.0;
    matrix.diag(2)  = 3.0;
    matrix.lower(0) = 4.0;
    matrix.lower(1) = 15.0;
    matrix.upper(0) = 6.0;
    matrix.upper(1) = 7.0;

    std::vector<double> vec1({4,2,9});
    auto res1 = matrix * vec1;
    REQUIRE(res1[0] == 16);
    REQUIRE(res1[1] == 83);
    REQUIRE(res1[2] == 57);
}

TEST_CASE("M-V-Multi test 2", "[TDMatrix]")
{
    TDMatrix matrix(3);
    matrix.diag(0)  = 1.0;
    matrix.diag(1)  = 2.0;
    matrix.diag(2)  = 3.0;
    matrix.lower(0) = 4.0;
    matrix.lower(1) = 15.0;
    matrix.upper(0) = 6.0;
    matrix.upper(1) = 7.0;

    std::vector<double> vec1({4,2,9});
    matrix.inlineMultiply(vec1);
    REQUIRE(vec1[0] == 16);
    REQUIRE(vec1[1] == 83);
    REQUIRE(vec1[2] == 57);
}
