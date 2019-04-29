// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// RandomStart_cpp
List RandomStart_cpp(const arma::vec& input_vector);
RcppExport SEXP _ClusterSSCA_RandomStart_cpp(SEXP input_vectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type input_vector(input_vectorSEXP);
    rcpp_result_gen = Rcpp::wrap(RandomStart_cpp(input_vector));
    return rcpp_result_gen;
END_RCPP
}
// MatrixCenter_cpp
arma::mat MatrixCenter_cpp(const arma::mat& input_matrix, const int& center, const int& scale);
RcppExport SEXP _ClusterSSCA_MatrixCenter_cpp(SEXP input_matrixSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type input_matrix(input_matrixSEXP);
    Rcpp::traits::input_parameter< const int& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const int& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(MatrixCenter_cpp(input_matrix, center, scale));
    return rcpp_result_gen;
END_RCPP
}
// sca_common_cpp
List sca_common_cpp(arma::mat data_con, int common);
RcppExport SEXP _ClusterSSCA_sca_common_cpp(SEXP data_conSEXP, SEXP commonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_con(data_conSEXP);
    Rcpp::traits::input_parameter< int >::type common(commonSEXP);
    rcpp_result_gen = Rcpp::wrap(sca_common_cpp(data_con, common));
    return rcpp_result_gen;
END_RCPP
}
// csca_cpp
List csca_cpp(arma::mat xc, arma::vec nvar, int nblock, int all_components, int ncluster, int iteration);
RcppExport SEXP _ClusterSSCA_csca_cpp(SEXP xcSEXP, SEXP nvarSEXP, SEXP nblockSEXP, SEXP all_componentsSEXP, SEXP nclusterSEXP, SEXP iterationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xc(xcSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< int >::type all_components(all_componentsSEXP);
    Rcpp::traits::input_parameter< int >::type ncluster(nclusterSEXP);
    Rcpp::traits::input_parameter< int >::type iteration(iterationSEXP);
    rcpp_result_gen = Rcpp::wrap(csca_cpp(xc, nvar, nblock, all_components, ncluster, iteration));
    return rcpp_result_gen;
END_RCPP
}
// sparsedisco_cpp
List sparsedisco_cpp(arma::mat data_con, arma::vec nvar, int all_components, arma::mat p, arma::uvec fixed_zeros, int n_zeros);
RcppExport SEXP _ClusterSSCA_sparsedisco_cpp(SEXP data_conSEXP, SEXP nvarSEXP, SEXP all_componentsSEXP, SEXP pSEXP, SEXP fixed_zerosSEXP, SEXP n_zerosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_con(data_conSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type all_components(all_componentsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fixed_zeros(fixed_zerosSEXP);
    Rcpp::traits::input_parameter< int >::type n_zeros(n_zerosSEXP);
    rcpp_result_gen = Rcpp::wrap(sparsedisco_cpp(data_con, nvar, all_components, p, fixed_zeros, n_zeros));
    return rcpp_result_gen;
END_RCPP
}
// IntSparseSca_rational_full_cpp
List IntSparseSca_rational_full_cpp(arma::mat data_con, arma::vec nvar, int nblock, int all_components, arma::mat initial, arma::uvec fixed_zeros, double sparsity);
RcppExport SEXP _ClusterSSCA_IntSparseSca_rational_full_cpp(SEXP data_conSEXP, SEXP nvarSEXP, SEXP nblockSEXP, SEXP all_componentsSEXP, SEXP initialSEXP, SEXP fixed_zerosSEXP, SEXP sparsitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_con(data_conSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< int >::type all_components(all_componentsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fixed_zeros(fixed_zerosSEXP);
    Rcpp::traits::input_parameter< double >::type sparsity(sparsitySEXP);
    rcpp_result_gen = Rcpp::wrap(IntSparseSca_rational_full_cpp(data_con, nvar, nblock, all_components, initial, fixed_zeros, sparsity));
    return rcpp_result_gen;
END_RCPP
}
// IntSparseSca_random_cpp
List IntSparseSca_random_cpp(arma::mat data_con, arma::vec nvar, int nblock, int all_components, arma::uvec fixed_zeros, double sparsity);
RcppExport SEXP _ClusterSSCA_IntSparseSca_random_cpp(SEXP data_conSEXP, SEXP nvarSEXP, SEXP nblockSEXP, SEXP all_componentsSEXP, SEXP fixed_zerosSEXP, SEXP sparsitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_con(data_conSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< int >::type all_components(all_componentsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fixed_zeros(fixed_zerosSEXP);
    Rcpp::traits::input_parameter< double >::type sparsity(sparsitySEXP);
    rcpp_result_gen = Rcpp::wrap(IntSparseSca_random_cpp(data_con, nvar, nblock, all_components, fixed_zeros, sparsity));
    return rcpp_result_gen;
END_RCPP
}
// cssca_quick_cpp
List cssca_quick_cpp(arma::mat data_con, arma::vec nvar, int nblock, int common, arma::vec distinct, int ncluster, int nrespondents, double sparsity, arma::vec feed, double cutoff_prop);
RcppExport SEXP _ClusterSSCA_cssca_quick_cpp(SEXP data_conSEXP, SEXP nvarSEXP, SEXP nblockSEXP, SEXP commonSEXP, SEXP distinctSEXP, SEXP nclusterSEXP, SEXP nrespondentsSEXP, SEXP sparsitySEXP, SEXP feedSEXP, SEXP cutoff_propSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_con(data_conSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type nblock(nblockSEXP);
    Rcpp::traits::input_parameter< int >::type common(commonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type distinct(distinctSEXP);
    Rcpp::traits::input_parameter< int >::type ncluster(nclusterSEXP);
    Rcpp::traits::input_parameter< int >::type nrespondents(nrespondentsSEXP);
    Rcpp::traits::input_parameter< double >::type sparsity(sparsitySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type feed(feedSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff_prop(cutoff_propSEXP);
    rcpp_result_gen = Rcpp::wrap(cssca_quick_cpp(data_con, nvar, nblock, common, distinct, ncluster, nrespondents, sparsity, feed, cutoff_prop));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ClusterSSCA_RandomStart_cpp", (DL_FUNC) &_ClusterSSCA_RandomStart_cpp, 1},
    {"_ClusterSSCA_MatrixCenter_cpp", (DL_FUNC) &_ClusterSSCA_MatrixCenter_cpp, 3},
    {"_ClusterSSCA_sca_common_cpp", (DL_FUNC) &_ClusterSSCA_sca_common_cpp, 2},
    {"_ClusterSSCA_csca_cpp", (DL_FUNC) &_ClusterSSCA_csca_cpp, 6},
    {"_ClusterSSCA_sparsedisco_cpp", (DL_FUNC) &_ClusterSSCA_sparsedisco_cpp, 6},
    {"_ClusterSSCA_IntSparseSca_rational_full_cpp", (DL_FUNC) &_ClusterSSCA_IntSparseSca_rational_full_cpp, 7},
    {"_ClusterSSCA_IntSparseSca_random_cpp", (DL_FUNC) &_ClusterSSCA_IntSparseSca_random_cpp, 6},
    {"_ClusterSSCA_cssca_quick_cpp", (DL_FUNC) &_ClusterSSCA_cssca_quick_cpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_ClusterSSCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}