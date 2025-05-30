#include <Rcpp.h>
using namespace Rcpp;

//**********************************************************************
//' Hazard rate function - CON
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' 
//' @description Hazard rate function for constant hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_con_c(NumericVector tvec,
                        double lambda0 
) {
 int n = tvec.size();
 NumericVector lambdavec(n);
 int i;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
 }
 
 return lambdavec;
}
 
//**********************************************************************
//' Integrated Hazard rate function - CON
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' 
//' @description Integrated hazard rate function for constant hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_con_c(NumericVector tvec,
                         double lambda0 
) {
  int n = tvec.size();
  NumericVector clambdavec(n);
  int i;
  
  for(i=0; i<n; i++) {
    clambdavec[i] = lambda0*tvec[i];
  }
  
  return clambdavec;
}

//**********************************************************************
//' Inverse survival function - CON
//' 
//' @param uvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' 
//' @description Inverse survival function for constant hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector invsurvf_con_c(NumericVector uvec,
                             double lambda0
                             ) {
 int n = uvec.size();
 NumericVector tvec(n);
 int i;
 
 for(i=0; i<n; i++) {
   tvec[i] = -log(uvec[i])/lambda0;
 }
 
 return tvec;
}


 