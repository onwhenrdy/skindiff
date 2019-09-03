#include "../../algorithms.h"
#include "../../tdmatrix.h"
#include "catch.hpp"

using namespace sc;

TEST_CASE("ThomasIP test 1", "[ThomasAlgorithm]")
{
    TDMatrix matrix(5);
    matrix.diag(0) = 1.0;
    matrix.diag(1) = 2.0;
    matrix.diag(2) = 3.0;
    matrix.diag(3) = 4.0;
    matrix.diag(4) = 5.0;

    matrix.lower(0) = 2.0;
    matrix.lower(1) = 3.0;
    matrix.lower(2) = 4.0;
    matrix.lower(3) = 5.0;

    matrix.upper(0) = 2.0;
    matrix.upper(1) = 3.0;
    matrix.upper(2) = 4.0;
    matrix.upper(3) = 5.0;

    std::vector<double> b({5, 15, 31, 53, 45});

    algorithm::thomasIP(matrix, b);

    REQUIRE(b[0] == Approx(1.0));
    REQUIRE(b[1] == Approx(2.0));
    REQUIRE(b[2] == Approx(3.0));
    REQUIRE(b[3] == Approx(4.0));
    REQUIRE(b[4] == Approx(5.0));
}

TEST_CASE("GaussPivot test 1", "[GaussPivot]")
{
    TDMatrix matrix(5);
    matrix.diag(0) = 1.0;
    matrix.diag(1) = 2.0;
    matrix.diag(2) = 3.0;
    matrix.diag(3) = 4.0;
    matrix.diag(4) = 5.0;

    matrix.lower(0) = 2.0;
    matrix.lower(1) = 3.0;
    matrix.lower(2) = 4.0;
    matrix.lower(3) = 5.0;

    matrix.upper(0) = 2.0;
    matrix.upper(1) = 3.0;
    matrix.upper(2) = 4.0;
    matrix.upper(3) = 5.0;

    std::vector<double> b({5, 15, 31, 53, 45});

    algorithm::gaussPivotIP(matrix, b);

    REQUIRE(b[0] == Approx(1.0));
    REQUIRE(b[1] == Approx(2.0));
    REQUIRE(b[2] == Approx(3.0));
    REQUIRE(b[3] == Approx(4.0));
    REQUIRE(b[4] == Approx(5.0));
}

TEST_CASE("GaussReUsePivot test 1", "[GaussReUsePivot]")
{
    TDMatrix matrix(5);
    matrix.diag(0) = 1.0;
    matrix.diag(1) = 2.0;
    matrix.diag(2) = 3.0;
    matrix.diag(3) = 4.0;
    matrix.diag(4) = 5.0;

    matrix.lower(0) = 2.0;
    matrix.lower(1) = 3.0;
    matrix.lower(2) = 4.0;
    matrix.lower(3) = 5.0;

    matrix.upper(0) = 2.0;
    matrix.upper(1) = 3.0;
    matrix.upper(2) = 4.0;
    matrix.upper(3) = 5.0;

    std::vector<double> b({5, 15, 31, 53, 45});

    algorithm::gaussReUsePivotIP(matrix, b);

    REQUIRE(b[0] == Approx(1.0));
    REQUIRE(b[1] == Approx(2.0));
    REQUIRE(b[2] == Approx(3.0));
    REQUIRE(b[3] == Approx(4.0));
    REQUIRE(b[4] == Approx(5.0));

    std::vector<double> c({5, 15, 31, 53, 45});

    algorithm::gaussReUsePivotIP(matrix, c);

    REQUIRE(c[0] == Approx(1.0));
    REQUIRE(c[1] == Approx(2.0));
    REQUIRE(c[2] == Approx(3.0));
    REQUIRE(c[3] == Approx(4.0));
    REQUIRE(c[4] == Approx(5.0));
}

TEST_CASE("ThomasReUseIP test 1", "[ThomasReUseIP]")
{
    TDMatrix matrix(5);
    matrix.diag(0) = 1.0;
    matrix.diag(1) = 2.0;
    matrix.diag(2) = 3.0;
    matrix.diag(3) = 4.0;
    matrix.diag(4) = 5.0;

    matrix.lower(0) = 2.0;
    matrix.lower(1) = 3.0;
    matrix.lower(2) = 4.0;
    matrix.lower(3) = 5.0;

    matrix.upper(0) = 2.0;
    matrix.upper(1) = 3.0;
    matrix.upper(2) = 4.0;
    matrix.upper(3) = 5.0;

    std::vector<double> b({5, 15, 31, 53, 45});

    algorithm::thomasReUseIP(matrix, b);

    REQUIRE(b[0] == Approx(1.0));
    REQUIRE(b[1] == Approx(2.0));
    REQUIRE(b[2] == Approx(3.0));
    REQUIRE(b[3] == Approx(4.0));
    REQUIRE(b[4] == Approx(5.0));

    std::vector<double> c({5, 15, 31, 53, 45});

    algorithm::thomasReUseIP(matrix, c);

    REQUIRE(c[0] == Approx(1.0));
    REQUIRE(c[1] == Approx(2.0));
    REQUIRE(c[2] == Approx(3.0));
    REQUIRE(c[3] == Approx(4.0));
    REQUIRE(c[4] == Approx(5.0));
}
