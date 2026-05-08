// testthat-based C++ unit tests. Compiled into the package and exposed
// to R via .cpp_run_tests() (see RcppExports). Other test-*.cpp files
// register their own contexts; only this file emits the actual runner.
//
// Rcpp.h MUST be included before testthat.h here: testthat.h drags in
// Rinternals.h, which defines `isNull` etc. as macros. If Rcpp.h comes
// after, those macros eat Rcpp's member functions of the same name.
#include <Rcpp.h>

#define TESTTHAT_TEST_RUNNER
#include <testthat.h>

extern "C" SEXP run_testthat_tests(SEXP use_xml_sxp);

// [[Rcpp::export(name = ".cpp_run_tests", rng = false)]]
Rcpp::RObject cpp_run_tests()
{
    return run_testthat_tests(Rf_ScalarLogical(0));
}
