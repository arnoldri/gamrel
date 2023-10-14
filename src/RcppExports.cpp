// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// csqrt
NumericVector csqrt(NumericVector x);
RcppExport SEXP _gamrel_csqrt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(csqrt(x));
    return rcpp_result_gen;
END_RCPP
}
// order1_c
IntegerVector order1_c(NumericVector x);
RcppExport SEXP _gamrel_order1_c(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order1_c(x));
    return rcpp_result_gen;
END_RCPP
}
// mean_c
double mean_c(NumericVector x);
RcppExport SEXP _gamrel_mean_c(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_c(x));
    return rcpp_result_gen;
END_RCPP
}
// makev_c
NumericVector makev_c(NumericVector uvec);
RcppExport SEXP _gamrel_makev_c(SEXP uvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type uvec(uvecSEXP);
    rcpp_result_gen = Rcpp::wrap(makev_c(uvec));
    return rcpp_result_gen;
END_RCPP
}
// lprior_ifrdfr_c
NumericVector lprior_ifrdfr_c(double eta, double gamma, NumericVector thetavec, NumericVector vvec, double alpha, double beta, double phi, double nu, double a1, double a2, double b1, double b2, double f1, double f2);
RcppExport SEXP _gamrel_lprior_ifrdfr_c(SEXP etaSEXP, SEXP gammaSEXP, SEXP thetavecSEXP, SEXP vvecSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP phiSEXP, SEXP nuSEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP f1SEXP, SEXP f2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec(vvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< double >::type f2(f2SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_ifrdfr_c(eta, gamma, thetavec, vvec, alpha, beta, phi, nu, a1, a2, b1, b2, f1, f2));
    return rcpp_result_gen;
END_RCPP
}
// lprior_lwb_c
NumericVector lprior_lwb_c(double eta, double a, double gamma, NumericVector thetavec, NumericVector vvec, double alpha, double beta, double phi, double nu, double c1, double c2, double a1, double a2, double b1, double b2, double f1, double f2);
RcppExport SEXP _gamrel_lprior_lwb_c(SEXP etaSEXP, SEXP aSEXP, SEXP gammaSEXP, SEXP thetavecSEXP, SEXP vvecSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP phiSEXP, SEXP nuSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP f1SEXP, SEXP f2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec(vvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< double >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< double >::type f2(f2SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_lwb_c(eta, a, gamma, thetavec, vvec, alpha, beta, phi, nu, c1, c2, a1, a2, b1, b2, f1, f2));
    return rcpp_result_gen;
END_RCPP
}
// lprior_sbt_c
NumericVector lprior_sbt_c(double eta, double gamma1, NumericVector thetavec1, NumericVector vvec1, double alpha1, double beta1, double phi1, double gamma2, NumericVector thetavec2, NumericVector vvec2, double alpha2, double beta2, double phi2, double nu, double a1, double a2, double b1, double b2, double f11, double f21, double f12, double f22);
RcppExport SEXP _gamrel_lprior_sbt_c(SEXP etaSEXP, SEXP gamma1SEXP, SEXP thetavec1SEXP, SEXP vvec1SEXP, SEXP alpha1SEXP, SEXP beta1SEXP, SEXP phi1SEXP, SEXP gamma2SEXP, SEXP thetavec2SEXP, SEXP vvec2SEXP, SEXP alpha2SEXP, SEXP beta2SEXP, SEXP phi2SEXP, SEXP nuSEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP f11SEXP, SEXP f21SEXP, SEXP f12SEXP, SEXP f22SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec1(thetavec1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec1(vvec1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec2(thetavec2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec2(vvec2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< double >::type phi2(phi2SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type f11(f11SEXP);
    Rcpp::traits::input_parameter< double >::type f21(f21SEXP);
    Rcpp::traits::input_parameter< double >::type f12(f12SEXP);
    Rcpp::traits::input_parameter< double >::type f22(f22SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_sbt_c(eta, gamma1, thetavec1, vvec1, alpha1, beta1, phi1, gamma2, thetavec2, vvec2, alpha2, beta2, phi2, nu, a1, a2, b1, b2, f11, f21, f12, f22));
    return rcpp_result_gen;
END_RCPP
}
// lprior_mbt_c
NumericVector lprior_mbt_c(double pival, double eta1, double gamma1, NumericVector thetavec1, NumericVector vvec1, double alpha1, double beta1, double phi1, double eta2, double gamma2, NumericVector thetavec2, NumericVector vvec2, double alpha2, double beta2, double phi2, double nu, double a1, double a2, double b1, double b2, double f11, double f21, double f12, double f22);
RcppExport SEXP _gamrel_lprior_mbt_c(SEXP pivalSEXP, SEXP eta1SEXP, SEXP gamma1SEXP, SEXP thetavec1SEXP, SEXP vvec1SEXP, SEXP alpha1SEXP, SEXP beta1SEXP, SEXP phi1SEXP, SEXP eta2SEXP, SEXP gamma2SEXP, SEXP thetavec2SEXP, SEXP vvec2SEXP, SEXP alpha2SEXP, SEXP beta2SEXP, SEXP phi2SEXP, SEXP nuSEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP f11SEXP, SEXP f21SEXP, SEXP f12SEXP, SEXP f22SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type pival(pivalSEXP);
    Rcpp::traits::input_parameter< double >::type eta1(eta1SEXP);
    Rcpp::traits::input_parameter< double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec1(thetavec1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec1(vvec1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< double >::type eta2(eta2SEXP);
    Rcpp::traits::input_parameter< double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec2(thetavec2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec2(vvec2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< double >::type phi2(phi2SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type f11(f11SEXP);
    Rcpp::traits::input_parameter< double >::type f21(f21SEXP);
    Rcpp::traits::input_parameter< double >::type f12(f12SEXP);
    Rcpp::traits::input_parameter< double >::type f22(f22SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_mbt_c(pival, eta1, gamma1, thetavec1, vvec1, alpha1, beta1, phi1, eta2, gamma2, thetavec2, vvec2, alpha2, beta2, phi2, nu, a1, a2, b1, b2, f11, f21, f12, f22));
    return rcpp_result_gen;
END_RCPP
}
// lprior_lcv_c
NumericVector lprior_lcv_c(double lambda0, double w0, double gamma, NumericVector thetavec, NumericVector vvec, double alpha, double beta, double phi, double s1, double s2, double sigmap_w0, double a1, double a2, double b1, double b2, double f1, double f2);
RcppExport SEXP _gamrel_lprior_lcv_c(SEXP lambda0SEXP, SEXP w0SEXP, SEXP gammaSEXP, SEXP thetavecSEXP, SEXP vvecSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP phiSEXP, SEXP s1SEXP, SEXP s2SEXP, SEXP sigmap_w0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP f1SEXP, SEXP f2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec(vvecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type sigmap_w0(sigmap_w0SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< double >::type f2(f2SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_lcv_c(lambda0, w0, gamma, thetavec, vvec, alpha, beta, phi, s1, s2, sigmap_w0, a1, a2, b1, b2, f1, f2));
    return rcpp_result_gen;
END_RCPP
}
// lprior_cvx_c
NumericVector lprior_cvx_c(double eta, double tau, double gamma1, NumericVector thetavec1, NumericVector vvec1, double alpha1, double beta1, double phi1, double gamma2, NumericVector thetavec2, NumericVector vvec2, double alpha2, double beta2, double phi2, double nu, double c1, double c2, double a1, double a2, double b1, double b2, double f1, double f2);
RcppExport SEXP _gamrel_lprior_cvx_c(SEXP etaSEXP, SEXP tauSEXP, SEXP gamma1SEXP, SEXP thetavec1SEXP, SEXP vvec1SEXP, SEXP alpha1SEXP, SEXP beta1SEXP, SEXP phi1SEXP, SEXP gamma2SEXP, SEXP thetavec2SEXP, SEXP vvec2SEXP, SEXP alpha2SEXP, SEXP beta2SEXP, SEXP phi2SEXP, SEXP nuSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP f1SEXP, SEXP f2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec1(thetavec1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec1(vvec1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type phi1(phi1SEXP);
    Rcpp::traits::input_parameter< double >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec2(thetavec2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vvec2(vvec2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< double >::type phi2(phi2SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< double >::type c2(c2SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< double >::type f2(f2SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_cvx_c(eta, tau, gamma1, thetavec1, vvec1, alpha1, beta1, phi1, gamma2, thetavec2, vvec2, alpha2, beta2, phi2, nu, c1, c2, a1, a2, b1, b2, f1, f2));
    return rcpp_result_gen;
END_RCPP
}
// lprior_mew_c
NumericVector lprior_mew_c(double lambda, double alpha, double theta, double gamma, double s1, double s2, double a1, double a2, double t1, double t2, double g1, double g2);
RcppExport SEXP _gamrel_lprior_mew_c(SEXP lambdaSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP gammaSEXP, SEXP s1SEXP, SEXP s2SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP t1SEXP, SEXP t2SEXP, SEXP g1SEXP, SEXP g2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< double >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< double >::type g2(g2SEXP);
    rcpp_result_gen = Rcpp::wrap(lprior_mew_c(lambda, alpha, theta, gamma, s1, s2, a1, a2, t1, t2, g1, g2));
    return rcpp_result_gen;
END_RCPP
}
// lambda_func_ifr_c
NumericVector lambda_func_ifr_c(NumericVector tvec, double lambda0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_lambda_func_ifr_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_func_ifr_c(tvec, lambda0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// lambda_func_dfr_c
NumericVector lambda_func_dfr_c(NumericVector tvec, double lambda0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_lambda_func_dfr_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_func_dfr_c(tvec, lambda0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// lambda_func_lwb_c
NumericVector lambda_func_lwb_c(NumericVector tvec, double lambda0, double a, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_lambda_func_lwb_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP aSEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_func_lwb_c(tvec, lambda0, a, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// lambda_func_lcv_c
NumericVector lambda_func_lcv_c(NumericVector tvec, double lambda0, double w0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_lambda_func_lcv_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP w0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_func_lcv_c(tvec, lambda0, w0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// lambda_func_cvx_c
NumericVector lambda_func_cvx_c(NumericVector tvec, double lambda0, double tau, NumericVector thetavec1, NumericVector wvec1, NumericVector thetavec2, NumericVector wvec2);
RcppExport SEXP _gamrel_lambda_func_cvx_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP tauSEXP, SEXP thetavec1SEXP, SEXP wvec1SEXP, SEXP thetavec2SEXP, SEXP wvec2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec1(thetavec1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec1(wvec1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec2(thetavec2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec2(wvec2SEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_func_cvx_c(tvec, lambda0, tau, thetavec1, wvec1, thetavec2, wvec2));
    return rcpp_result_gen;
END_RCPP
}
// int_lambda_func_ifr_c
NumericVector int_lambda_func_ifr_c(NumericVector tvec, double lambda0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_int_lambda_func_ifr_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(int_lambda_func_ifr_c(tvec, lambda0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// int_lambda_func_dfr_c
NumericVector int_lambda_func_dfr_c(NumericVector tvec, double lambda0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_int_lambda_func_dfr_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(int_lambda_func_dfr_c(tvec, lambda0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// int_lambda_func_lwb_c
NumericVector int_lambda_func_lwb_c(NumericVector tvec, double lambda0, double a, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_int_lambda_func_lwb_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP aSEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(int_lambda_func_lwb_c(tvec, lambda0, a, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// int_lambda_func_lcv_c
NumericVector int_lambda_func_lcv_c(NumericVector tvec, double lambda0, double w0, NumericVector thetavec, NumericVector wvec, double epsilon);
RcppExport SEXP _gamrel_int_lambda_func_lcv_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP w0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(int_lambda_func_lcv_c(tvec, lambda0, w0, thetavec, wvec, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// both_lambda_func_ifr_c
NumericMatrix both_lambda_func_ifr_c(NumericVector tvec, double lambda0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_both_lambda_func_ifr_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(both_lambda_func_ifr_c(tvec, lambda0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// both_lambda_func_dfr_c
NumericMatrix both_lambda_func_dfr_c(NumericVector tvec, double lambda0, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_both_lambda_func_dfr_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(both_lambda_func_dfr_c(tvec, lambda0, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// both_lambda_func_lwb_c
NumericMatrix both_lambda_func_lwb_c(NumericVector tvec, double lambda0, double a, NumericVector thetavec, NumericVector wvec);
RcppExport SEXP _gamrel_both_lambda_func_lwb_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP aSEXP, SEXP thetavecSEXP, SEXP wvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    rcpp_result_gen = Rcpp::wrap(both_lambda_func_lwb_c(tvec, lambda0, a, thetavec, wvec));
    return rcpp_result_gen;
END_RCPP
}
// both_lambda_func_lcv_c
NumericMatrix both_lambda_func_lcv_c(NumericVector tvec, double lambda0, double w0, NumericVector thetavec, NumericVector wvec, double epsilon);
RcppExport SEXP _gamrel_both_lambda_func_lcv_c(SEXP tvecSEXP, SEXP lambda0SEXP, SEXP w0SEXP, SEXP thetavecSEXP, SEXP wvecSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tvec(tvecSEXP);
    Rcpp::traits::input_parameter< double >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetavec(thetavecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wvec(wvecSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(both_lambda_func_lcv_c(tvec, lambda0, w0, thetavec, wvec, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// shellSort_c
NumericVector shellSort_c(NumericVector x);
RcppExport SEXP _gamrel_shellSort_c(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(shellSort_c(x));
    return rcpp_result_gen;
END_RCPP
}
// order_c
IntegerVector order_c(NumericVector x);
RcppExport SEXP _gamrel_order_c(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_c(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _gamrel_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gamrel_csqrt", (DL_FUNC) &_gamrel_csqrt, 1},
    {"_gamrel_order1_c", (DL_FUNC) &_gamrel_order1_c, 1},
    {"_gamrel_mean_c", (DL_FUNC) &_gamrel_mean_c, 1},
    {"_gamrel_makev_c", (DL_FUNC) &_gamrel_makev_c, 1},
    {"_gamrel_lprior_ifrdfr_c", (DL_FUNC) &_gamrel_lprior_ifrdfr_c, 14},
    {"_gamrel_lprior_lwb_c", (DL_FUNC) &_gamrel_lprior_lwb_c, 17},
    {"_gamrel_lprior_sbt_c", (DL_FUNC) &_gamrel_lprior_sbt_c, 22},
    {"_gamrel_lprior_mbt_c", (DL_FUNC) &_gamrel_lprior_mbt_c, 24},
    {"_gamrel_lprior_lcv_c", (DL_FUNC) &_gamrel_lprior_lcv_c, 17},
    {"_gamrel_lprior_cvx_c", (DL_FUNC) &_gamrel_lprior_cvx_c, 23},
    {"_gamrel_lprior_mew_c", (DL_FUNC) &_gamrel_lprior_mew_c, 12},
    {"_gamrel_lambda_func_ifr_c", (DL_FUNC) &_gamrel_lambda_func_ifr_c, 4},
    {"_gamrel_lambda_func_dfr_c", (DL_FUNC) &_gamrel_lambda_func_dfr_c, 4},
    {"_gamrel_lambda_func_lwb_c", (DL_FUNC) &_gamrel_lambda_func_lwb_c, 5},
    {"_gamrel_lambda_func_lcv_c", (DL_FUNC) &_gamrel_lambda_func_lcv_c, 5},
    {"_gamrel_lambda_func_cvx_c", (DL_FUNC) &_gamrel_lambda_func_cvx_c, 7},
    {"_gamrel_int_lambda_func_ifr_c", (DL_FUNC) &_gamrel_int_lambda_func_ifr_c, 4},
    {"_gamrel_int_lambda_func_dfr_c", (DL_FUNC) &_gamrel_int_lambda_func_dfr_c, 4},
    {"_gamrel_int_lambda_func_lwb_c", (DL_FUNC) &_gamrel_int_lambda_func_lwb_c, 5},
    {"_gamrel_int_lambda_func_lcv_c", (DL_FUNC) &_gamrel_int_lambda_func_lcv_c, 6},
    {"_gamrel_both_lambda_func_ifr_c", (DL_FUNC) &_gamrel_both_lambda_func_ifr_c, 4},
    {"_gamrel_both_lambda_func_dfr_c", (DL_FUNC) &_gamrel_both_lambda_func_dfr_c, 4},
    {"_gamrel_both_lambda_func_lwb_c", (DL_FUNC) &_gamrel_both_lambda_func_lwb_c, 5},
    {"_gamrel_both_lambda_func_lcv_c", (DL_FUNC) &_gamrel_both_lambda_func_lcv_c, 6},
    {"_gamrel_shellSort_c", (DL_FUNC) &_gamrel_shellSort_c, 1},
    {"_gamrel_order_c", (DL_FUNC) &_gamrel_order_c, 1},
    {"_gamrel_rcpp_hello_world", (DL_FUNC) &_gamrel_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_gamrel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
