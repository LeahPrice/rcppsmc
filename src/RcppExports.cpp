// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// blockpfGaussianOpt_cpp
Rcpp::List blockpfGaussianOpt_cpp(Rcpp::NumericVector data, long inlNumber, long inlLag);
RcppExport SEXP RcppSMC_blockpfGaussianOpt_cpp(SEXP dataSEXP, SEXP inlNumberSEXP, SEXP inlLagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< long >::type inlNumber(inlNumberSEXP);
    Rcpp::traits::input_parameter< long >::type inlLag(inlLagSEXP);
    rcpp_result_gen = Rcpp::wrap(blockpfGaussianOpt_cpp(data, inlNumber, inlLag));
    return rcpp_result_gen;
END_RCPP
}
// pfLineartBS_cpp
Rcpp::List pfLineartBS_cpp(arma::mat data, unsigned long inlNumber, bool useF, Rcpp::Function f);
RcppExport SEXP RcppSMC_pfLineartBS_cpp(SEXP dataSEXP, SEXP inlNumberSEXP, SEXP useFSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type inlNumber(inlNumberSEXP);
    Rcpp::traits::input_parameter< bool >::type useF(useFSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(pfLineartBS_cpp(data, inlNumber, useF, f));
    return rcpp_result_gen;
END_RCPP
}
// pfNonlinBS_cpp
Rcpp::List pfNonlinBS_cpp(arma::vec data, long inlNumber);
RcppExport SEXP RcppSMC_pfNonlinBS_cpp(SEXP dataSEXP, SEXP inlNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type data(dataSEXP);
    Rcpp::traits::input_parameter< long >::type inlNumber(inlNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(pfNonlinBS_cpp(data, inlNumber));
    return rcpp_result_gen;
END_RCPP
}
// radiataBS_cpp
Rcpp::List radiataBS_cpp(arma::mat data, unsigned long inlNumber);
RcppExport SEXP RcppSMC_radiataBS_cpp(SEXP dataSEXP, SEXP inlNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type inlNumber(inlNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(radiataBS_cpp(data, inlNumber));
    return rcpp_result_gen;
END_RCPP
}
// radiataPPBS_cpp
Rcpp::List radiataPPBS_cpp(arma::mat data, arma::vec intemps, unsigned long inlNumber);
RcppExport SEXP RcppSMC_radiataPPBS_cpp(SEXP dataSEXP, SEXP intempsSEXP, SEXP inlNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type intemps(intempsSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type inlNumber(inlNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(radiataPPBS_cpp(data, intemps, inlNumber));
    return rcpp_result_gen;
END_RCPP
}