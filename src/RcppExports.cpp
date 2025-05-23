// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// updateCorrectCoefCpp
List updateCorrectCoefCpp(List X0_glist, List X1_glist, List ThetaList, NumericVector a_i, NumericVector b_i, double penal_ksi, double penal_gamma);
RcppExport SEXP _CordBat_updateCorrectCoefCpp(SEXP X0_glistSEXP, SEXP X1_glistSEXP, SEXP ThetaListSEXP, SEXP a_iSEXP, SEXP b_iSEXP, SEXP penal_ksiSEXP, SEXP penal_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X0_glist(X0_glistSEXP);
    Rcpp::traits::input_parameter< List >::type X1_glist(X1_glistSEXP);
    Rcpp::traits::input_parameter< List >::type ThetaList(ThetaListSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a_i(a_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b_i(b_iSEXP);
    Rcpp::traits::input_parameter< double >::type penal_ksi(penal_ksiSEXP);
    Rcpp::traits::input_parameter< double >::type penal_gamma(penal_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(updateCorrectCoefCpp(X0_glist, X1_glist, ThetaList, a_i, b_i, penal_ksi, penal_gamma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CordBat_updateCorrectCoefCpp", (DL_FUNC) &_CordBat_updateCorrectCoefCpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_CordBat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
