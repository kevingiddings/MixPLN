// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EM_loop
Rcpp::List EM_loop(Rcpp::List& start, arma::vec& lib_mat, Rcpp::List& m, Rcpp::List mu, Rcpp::List sigma, Rcpp::List isigma, Rcpp::List iOsigO, Rcpp::List& S, Rcpp::List& O_list, arma::mat& Y, arma::mat z, arma::mat& O_mat, Rcpp::List pi_g, const int& d, const int& N, const int& G);
RcppExport SEXP _DangLoop_EM_loop(SEXP startSEXP, SEXP lib_matSEXP, SEXP mSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP isigmaSEXP, SEXP iOsigOSEXP, SEXP SSEXP, SEXP O_listSEXP, SEXP YSEXP, SEXP zSEXP, SEXP O_matSEXP, SEXP pi_gSEXP, SEXP dSEXP, SEXP NSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type start(startSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lib_mat(lib_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type m(mSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type isigma(isigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type iOsigO(iOsigOSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type S(SSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type O_list(O_listSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type O_mat(O_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pi_g(pi_gSEXP);
    Rcpp::traits::input_parameter< const int& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(EM_loop(start, lib_mat, m, mu, sigma, isigma, iOsigO, S, O_list, Y, z, O_mat, pi_g, d, N, G));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _DangLoop_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _DangLoop_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _DangLoop_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _DangLoop_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
