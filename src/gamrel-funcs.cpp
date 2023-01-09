#include <Rcpp.h>
using namespace Rcpp;

//' C version of order() - bad performance with duplicates
//' 
//' @export
// [[Rcpp::export]]
IntegerVector order1_c(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}



//' C version of mean()
//' 
//' @export
// [[Rcpp::export]]
double mean_c(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

//' Hazard rate function - IFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' @param lambda0 Offset to be added to the hazard rate
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector lambda_func_ifr_c(NumericVector tvec,
                                double lambda0, 
                                NumericVector thetavec, 
                                NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    lambdavec[i] = lambda0;
    for(k=0; k<kmax; k++) {
      if(thetavec[k]<tvec[i]) lambdavec[i] += wvec[k];
    }
  }

  return lambdavec;
}

//' Hazard rate function - DFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector lambda_func_dfr_c(NumericVector tvec, 
                                double lambda0, 
                                NumericVector thetavec, 
                                NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    lambdavec[i] = lambda0;
    for(k=0; k<kmax; k++) {
      if(thetavec[k]>=tvec[i]) lambdavec[i] += wvec[k];
    }
  }
  
  return lambdavec;
}

//' Hazard rate function - LWB
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param a inflection point
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @export
// [[Rcpp::export]]
NumericVector lambda_func_lwb_c(NumericVector tvec, 
                                double lambda0, 
                                double a,
                                NumericVector thetavec, 
                                NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    lambdavec[i] = 0;
    for(k=0; k<kmax; k++) {
      if(thetavec[k]<=abs(tvec[i]-a)) lambdavec[i] += wvec[k];
    }
  }
  
  return lambdavec;
}

//' Integrated hazard rate function - IFR version
//' 
//' @param tvec Locations at which to evaluate the function
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' @param lambda0 Offset to be added to the hazard rate (ignored)
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector int_lambda_func_ifr_c(NumericVector tvec, 
                                    double lambda0, 
                                    NumericVector thetavec, 
                                    NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector int_lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    int_lambdavec[i] = lambda0*tvec[i];
    for(k=0; k<kmax; k++) {
      if(thetavec[k]<tvec[i]) {
        int_lambdavec[i] += wvec[k]*(tvec[i]-thetavec[k]); 
      }
    }
  }

  return int_lambdavec;
}

//' Integrated hazard rate function - DFR version
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate (ignored)
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector int_lambda_func_dfr_c(NumericVector tvec, 
                                    double lambda0, 
                                    NumericVector thetavec, 
                                    NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector int_lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    int_lambdavec[i] = lambda0*tvec[i];
    for(k=0; k<kmax; k++) {
      if(thetavec[k]>=tvec[i]) {
        int_lambdavec[i] += wvec[k]*tvec[i];
      } else {
        int_lambdavec[i] += wvec[k]*thetavec[k];
      }
    }
  }
  
  return int_lambdavec;
}

//' Integrated hazard rate function - LWB version
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate (ignored)
//' @param a inflection point
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector int_lambda_func_lwb_c(NumericVector tvec, 
                                    double lambda0,
                                    double a, 
                                    NumericVector thetavec, 
                                    NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector int_lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    int_lambdavec[i] = 0;
    for(k=0; k<kmax; k++) {
      if(tvec[i]<=a && thetavec[k]<a) {
        int_lambdavec[i] += wvec[k]*fmin(tvec[i],a-thetavec[k]);
      } else if(tvec[i]>a) {
        int_lambdavec[i] += wvec[k]*(fmax(0,a-thetavec[k])
                                    +fmax(0,tvec[i]-a-thetavec[k]));
      }
    }
  }
  
  return int_lambdavec;
}

//' Integrated hazard rate function - Log Convex version
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 hazard rate scaling
//' @param w0 hazard rate gradient offset
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Integrated hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector int_lambda_func_lcv_c(NumericVector tvec, 
                                    double lambda0, 
                                    double w0,
                                    NumericVector thetavec, 
                                    NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  int kmaxp2 = kmax+2;
  NumericVector int_lambdavec(n);
  IntegerVector odx(kmax);
  NumericVector othetavec(kmaxp2);
  NumericVector owvec(kmaxp2);
  NumericVector s01vec(kmaxp2); 
  NumericVector s2vec(kmaxp2);
  NumericVector ccvec(kmaxp2);
  double cs;
  NumericVector ssvec(kmaxp2);
  int i,k,k1;
  
  odx = order1_c(thetavec);
  othetavec[0] = 0;
  owvec[0] = 0;
  for(k=0; k<kmax; k++) {
    othetavec[k+1] = thetavec[odx[k]-1];
    owvec[k+1] = wvec[odx[k]-1];
  }
  othetavec[kmaxp2-1] = 1.1*othetavec[kmaxp2-2];
  owvec[kmaxp2-1] = 0;

  s01vec[0] = w0;
  s2vec[0] = 0;
  for(k=1; k<kmaxp2; k++) {
    s01vec[k] = s01vec[k-1] + owvec[k];
    s2vec[k] = s2vec[k-1] + owvec[k]*othetavec[k];
  }
  
  ssvec[0] = 0;
  for(k=0; k<kmaxp2-1; k++) {
    ccvec[k] = lambda0*exp(-s2vec[k])/s01vec[k];
    cs = ccvec[k]*(exp(s01vec[k]*othetavec[k+1])-exp(s01vec[k]*othetavec[k]));
    ssvec[k+1] = ssvec[k] + cs;
  }
  ccvec[kmaxp2-1] = ccvec[kmaxp2-2];
  
  for(i=0; i<n; i++) {
    k1 = 0;
    while( (k1<kmaxp2-1) && (othetavec[k1+1]<=tvec[i]) ) {
      k1++;
    }
    int_lambdavec[i] = ssvec[k1] + ccvec[k1]*(exp(s01vec[k1]*tvec[i])-exp(s01vec[k1]*othetavec[k1]));
  }
  
  return int_lambdavec;
}

