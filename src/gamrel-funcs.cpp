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

//' Make DPP beta draws given a set of unscaled stick breaking weights
//' 
//' @export
// [[Rcpp::export]]
NumericVector makev_c(NumericVector uvec) {

  int kmax = uvec.size();
  NumericVector vvec(kmax);
  int k;
  double cp;
  
  vvec[0] = uvec[0];
  if(kmax>1) {
    cp = 1.0;
    for(k=1; k<kmax-1; k++) {
      vvec[k] = uvec[k]/(cp-uvec[k-1]);
      cp = cp*(1-vvec[k-1]);
    }
    vvec[kmax-1] = 0.5;
  }
  
  return vvec;
}


//**********************************************************************
//' Hazard rate function - IFR
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
      if(thetavec[k]<=tvec[i]) lambdavec[i] += wvec[k];
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
      if(thetavec[k]>tvec[i]) lambdavec[i] += wvec[k];
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
    lambdavec[i] = lambda0;
    for(k=0; k<kmax; k++) {
      if(thetavec[k]<=fabs(tvec[i]-a)) lambdavec[i] += wvec[k];
    }
  }
  
  return lambdavec;
}

// SBT - not needed
// MBT - not needed

//**********************************************************************
//' Hazard rate function - LCV
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 scaling
//' @param Offset to be added to the log hazard rate derivatice
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericVector lambda_func_lcv_c(NumericVector tvec,
                                double lambda0, 
                                double w0,
                                NumericVector thetavec, 
                                NumericVector wvec) {
  int n = tvec.size();
  int kmax = thetavec.size();
  NumericVector lambdavec(n);
  int i,k;
  
  for(i=0; i<n; i++) {
    lambdavec[i] = w0;
    for(k=0; k<kmax; k++) {
      if(thetavec[k]<=tvec[i]) lambdavec[i] += wvec[k];
    }
    lambdavec[i] = lambda0*exp(lambdavec[i]);
  }
  
  return lambdavec;
}



//------------------------------------------------------------

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
    int_lambdavec[i] = lambda0*tvec[i];
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

// SBT not needed
// MBT not needed

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
                                    NumericVector wvec,
                                    double epsilon) {  // needs 100*.Machine$double.neg.eps
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
    s01vec[k] = s01vec[k-1] + owvec[k];   // C
    s2vec[k] = s2vec[k-1] + owvec[k]*othetavec[k];  // D
  }
  
  ssvec[0] = 0;
  for(k=0; k<kmaxp2-1; k++) {
    ccvec[k] = lambda0*exp(-s2vec[k]);
    if(s01vec[k]==0) {
      cs = ccvec[k]*(othetavec[k+1]-othetavec[k]);
    } else {
      cs = (ccvec[k]/s01vec[k])*(exp(s01vec[k]*othetavec[k+1])-exp(s01vec[k]*othetavec[k]));
    }
    ssvec[k+1] = ssvec[k] + cs;
  }
  ccvec[kmaxp2-1] = ccvec[kmaxp2-2];
  
  // NumericVector excess(9);
  
  for(i=0; i<n; i++) {
    k1 = 0;
    while( (k1<kmaxp2-1) && (othetavec[k1+1]<=tvec[i]) ) {   // try < rather than <=?
      k1++;
    }
    int_lambdavec[i] = ssvec[k1];
    if(std::abs(s01vec[k1])<epsilon) {
      //excess[0] = -i;
      //excess[1] = k1;
      //excess[2] = othetavec[k1];
      //excess[3] = tvec[i];
      //excess[4] = othetavec[k1+1];
      //excess[5] = int_lambdavec[i];
      //excess[6] = ccvec[k1]*(tvec[i]-othetavec[k1]);
      //excess[7] = s01vec[k1];
      //excess[8] = ccvec[k1];
      
      int_lambdavec[i] += ccvec[k1]*(tvec[i]-othetavec[k1]);
    } else {
      //excess[0] = i;
      //excess[1] = k1;
      //excess[2] = othetavec[k1];
      //excess[3] = tvec[i];
      //excess[4] = othetavec[k1+1];
      //excess[5] = int_lambdavec[i];
      //excess[6] = (ccvec[k1]/s01vec[k1])*(exp(s01vec[k1]*tvec[i])-exp(s01vec[k1]*othetavec[k1]));
      //excess[7] = s01vec[k1];
      //excess[8] = ccvec[k1];
      
      int_lambdavec[i] += (ccvec[k1]/s01vec[k1])*(exp(s01vec[k1]*tvec[i])-exp(s01vec[k1]*othetavec[k1]));
    }
    //Rcpp::Rcout << excess << std::endl; // !!==
  }
  
  //Rcpp::Rcout << "s01vec" << std::endl; // !!==
  //Rcpp::Rcout << s01vec << std::endl; // **!!==
  //NumericVector excess(n);
  //Rcpp::Rcout << std::numeric_limits::epsilon( ) << std::endl;

  return int_lambdavec;
}

//*************************************************************************
//' Hazard and Integrated Hazard rate functions - IFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Hazard rate and integrated hazard rate functions generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix both_lambda_func_ifr_c(NumericVector tvec,
                                     double lambda0, 
                                     NumericVector thetavec, 
                                     NumericVector wvec) {
   int n = tvec.size();
   int kmax = thetavec.size();
   NumericMatrix bothlambdamat(n,2);
   int i,k;
   
   for(i=0; i<n; i++) {
     bothlambdamat(i,0) = lambda0;
     bothlambdamat(i,1) = lambda0*tvec[i];
     for(k=0; k<kmax; k++) {
       if(thetavec[k]<=tvec[i]) {
         bothlambdamat(i,0) += wvec[k];
         bothlambdamat(i,1) += wvec[k]*(tvec[i]-thetavec[k]); 
       }
     }
   }

   return bothlambdamat;
}

//' Hazard rate and integrated hazard rate functions - DFR
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
NumericMatrix both_lambda_func_dfr_c(NumericVector tvec, 
                                     double lambda0, 
                                     NumericVector thetavec, 
                                     NumericVector wvec) {
   int n = tvec.size();
   int kmax = thetavec.size();
   NumericMatrix bothlambdamat(n,2);
   int i,k;
   
   for(i=0; i<n; i++) {
     bothlambdamat(i,0) = lambda0;
     bothlambdamat(i,1) = lambda0*tvec[i];
     for(k=0; k<kmax; k++) {
       if(thetavec[k]>=tvec[i]) {
         bothlambdamat(i,0) += wvec[k];
         bothlambdamat(i,1) += wvec[k]*tvec[i];
       } else {
         bothlambdamat(i,1) += wvec[k]*thetavec[k];
       }
     }
   }

   return bothlambdamat;
}

//' Hazard and integrated hazard rate functions - LWB
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param a inflection point
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix both_lambda_func_lwb_c(NumericVector tvec, 
                                     double lambda0, 
                                     double a,
                                     NumericVector thetavec, 
                                     NumericVector wvec) {
   int n = tvec.size();
   int kmax = thetavec.size();
   NumericMatrix bothlambdamat(n,2);
   int i,k;
   
   for(i=0; i<n; i++) {
     bothlambdamat(i,0) = lambda0;
     bothlambdamat(i,1) = lambda0*tvec[i];
     for(k=0; k<kmax; k++) {
       if(thetavec[k]<=fabs(tvec[i]-a)) {
         bothlambdamat(i,0) += wvec[k];
       } 
       if(tvec[i]<=a && thetavec[k]<a) {
         bothlambdamat(i,1) += wvec[k]*fmin(tvec[i],a-thetavec[k]);
       } else if(tvec[i]>a) {
         bothlambdamat(i,1) += wvec[k]*(fmax(0,a-thetavec[k])
                                       +fmax(0,tvec[i]-a-thetavec[k]));
       }
     }
   }

   return bothlambdamat;
}

//**********************************************************************
//' Hazard and integrated hazard rate functions - LCV
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 scaling
//' @param Offset to be added to the log hazard rate derivatice
//' @param thetavec Set of support locations
//' @param wvec Set of associated weights
//' 
//' @description Hazard rate function generated by 
//' a set of discrete locations and weights and integrated
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix both_lambda_func_lcv_c(NumericVector tvec,
                                     double lambda0, 
                                     double w0,
                                     NumericVector thetavec, 
                                     NumericVector wvec,
                                     double epsilon) {  // needs 100*.Machine$double.neg.eps
  int n = tvec.size();
  int kmax = thetavec.size();
  int kmaxp2 = kmax+2;
  NumericMatrix bothlambdamat(n,2);
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
    s01vec[k] = s01vec[k-1] + owvec[k];   // C
    s2vec[k] = s2vec[k-1] + owvec[k]*othetavec[k];  // D
  }
  
  ssvec[0] = 0;
  for(k=0; k<kmaxp2-1; k++) {
    ccvec[k] = lambda0*exp(-s2vec[k]);
    if(s01vec[k]==0) {
      cs = ccvec[k]*(othetavec[k+1]-othetavec[k]);
    } else {
      cs = (ccvec[k]/s01vec[k])*(exp(s01vec[k]*othetavec[k+1])-exp(s01vec[k]*othetavec[k]));
    }
    ssvec[k+1] = ssvec[k] + cs;
  }
  ccvec[kmaxp2-1] = ccvec[kmaxp2-2];

  
  for(i=0; i<n; i++) {
    bothlambdamat(i,0) = w0;
    k1 = 0;
    while( (k1<kmaxp2-1) && (othetavec[k1+1]<=tvec[i]) ) {   // try < rather than <=?
      bothlambdamat(i,0) += wvec[k1];
      k1++;
    }
    bothlambdamat(i,1) = ssvec[k1];
    if(std::abs(s01vec[k1])<epsilon) {
      bothlambdamat(i,1) += ccvec[k1]*(tvec[i]-othetavec[k1]);
    } else {
      bothlambdamat(i,1) += (ccvec[k1]/s01vec[k1])*(exp(s01vec[k1]*tvec[i])-exp(s01vec[k1]*othetavec[k1]));
    }
    bothlambdamat(i,0) = lambda0*exp(bothlambdamat(i,0));
  }
  
  return bothlambdamat;
}



