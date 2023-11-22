// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// prior_loglik
double prior_loglik(NumericVector para);
RcppExport SEXP _seroreconstruct_prior_loglik(SEXP paraSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP);
    rcpp_result_gen = Rcpp::wrap(prior_loglik(para));
    return rcpp_result_gen;
END_RCPP
}
// sim_data
List sim_data(NumericMatrix data1, NumericMatrix ILI, NumericVector para, NumericVector para2);
RcppExport SEXP _seroreconstruct_sim_data(SEXP data1SEXP, SEXP ILISEXP, SEXP paraSEXP, SEXP para2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ILI(ILISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para2(para2SEXP);
    rcpp_result_gen = Rcpp::wrap(sim_data(data1, ILI, para, para2));
    return rcpp_result_gen;
END_RCPP
}
// loglik
List loglik(NumericMatrix data11, NumericMatrix data111, NumericMatrix data21, NumericMatrix ILI, NumericVector para, NumericVector para2, int level1, int level2, int level3, int season, NumericMatrix blankmatrix);
RcppExport SEXP _seroreconstruct_loglik(SEXP data11SEXP, SEXP data111SEXP, SEXP data21SEXP, SEXP ILISEXP, SEXP paraSEXP, SEXP para2SEXP, SEXP level1SEXP, SEXP level2SEXP, SEXP level3SEXP, SEXP seasonSEXP, SEXP blankmatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data11(data11SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data111(data111SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data21(data21SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ILI(ILISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para2(para2SEXP);
    Rcpp::traits::input_parameter< int >::type level1(level1SEXP);
    Rcpp::traits::input_parameter< int >::type level2(level2SEXP);
    Rcpp::traits::input_parameter< int >::type level3(level3SEXP);
    Rcpp::traits::input_parameter< int >::type season(seasonSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type blankmatrix(blankmatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik(data11, data111, data21, ILI, para, para2, level1, level2, level3, season, blankmatrix));
    return rcpp_result_gen;
END_RCPP
}
// all_update
List all_update(NumericMatrix data11, NumericMatrix data111, NumericMatrix data21, NumericMatrix ILI, NumericVector para, NumericVector para2, NumericMatrix loglik1, NumericMatrix loglik2, NumericMatrix loglik3);
RcppExport SEXP _seroreconstruct_all_update(SEXP data11SEXP, SEXP data111SEXP, SEXP data21SEXP, SEXP ILISEXP, SEXP paraSEXP, SEXP para2SEXP, SEXP loglik1SEXP, SEXP loglik2SEXP, SEXP loglik3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data11(data11SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data111(data111SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data21(data21SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ILI(ILISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para2(para2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loglik1(loglik1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loglik2(loglik2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loglik3(loglik3SEXP);
    rcpp_result_gen = Rcpp::wrap(all_update(data11, data111, data21, ILI, para, para2, loglik1, loglik2, loglik3));
    return rcpp_result_gen;
END_RCPP
}
// add_remove_infection
List add_remove_infection(NumericMatrix data11, NumericMatrix data111, NumericMatrix data21, NumericMatrix ILI, NumericVector para, NumericVector para2, NumericMatrix loglik1, NumericMatrix loglik2, NumericMatrix loglik3);
RcppExport SEXP _seroreconstruct_add_remove_infection(SEXP data11SEXP, SEXP data111SEXP, SEXP data21SEXP, SEXP ILISEXP, SEXP paraSEXP, SEXP para2SEXP, SEXP loglik1SEXP, SEXP loglik2SEXP, SEXP loglik3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data11(data11SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data111(data111SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data21(data21SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ILI(ILISEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para2(para2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loglik1(loglik1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loglik2(loglik2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type loglik3(loglik3SEXP);
    rcpp_result_gen = Rcpp::wrap(add_remove_infection(data11, data111, data21, ILI, para, para2, loglik1, loglik2, loglik3));
    return rcpp_result_gen;
END_RCPP
}
// mcmc
List mcmc(NumericMatrix input1, NumericMatrix input2, NumericMatrix input3, NumericMatrix ILI, int mcmc_n, NumericVector int_para, NumericVector int_para2, NumericVector int_para3, NumericVector paraseason, NumericVector move, NumericVector sigma, NumericVector sigma3, int burnin, int thinning);
RcppExport SEXP _seroreconstruct_mcmc(SEXP input1SEXP, SEXP input2SEXP, SEXP input3SEXP, SEXP ILISEXP, SEXP mcmc_nSEXP, SEXP int_paraSEXP, SEXP int_para2SEXP, SEXP int_para3SEXP, SEXP paraseasonSEXP, SEXP moveSEXP, SEXP sigmaSEXP, SEXP sigma3SEXP, SEXP burninSEXP, SEXP thinningSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type input1(input1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type input2(input2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type input3(input3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ILI(ILISEXP);
    Rcpp::traits::input_parameter< int >::type mcmc_n(mcmc_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type int_para(int_paraSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type int_para2(int_para2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type int_para3(int_para3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type paraseason(paraseasonSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type move(moveSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma3(sigma3SEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thinning(thinningSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc(input1, input2, input3, ILI, mcmc_n, int_para, int_para2, int_para3, paraseason, move, sigma, sigma3, burnin, thinning));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_seroreconstruct_prior_loglik", (DL_FUNC) &_seroreconstruct_prior_loglik, 1},
    {"_seroreconstruct_sim_data", (DL_FUNC) &_seroreconstruct_sim_data, 4},
    {"_seroreconstruct_loglik", (DL_FUNC) &_seroreconstruct_loglik, 11},
    {"_seroreconstruct_all_update", (DL_FUNC) &_seroreconstruct_all_update, 9},
    {"_seroreconstruct_add_remove_infection", (DL_FUNC) &_seroreconstruct_add_remove_infection, 9},
    {"_seroreconstruct_mcmc", (DL_FUNC) &_seroreconstruct_mcmc, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_seroreconstruct(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
