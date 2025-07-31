#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//--// [[Rcpp::depends(fntl)]]
//--#include "fntl.h"

//**********************************************************************
//* CON - Constant Hazard Rate
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


//**********************************************************************
//* IFR - Increasing hazard rate
//**********************************************************************
//' Hazard rate function - IFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Hazard rate function for IFR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_ifr_c(NumericVector tvec,
                         double lambda0, 
                         NumericVector thetavec, 
                         NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector lambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
   for(k=0; k<kmax; k++) {
     if(thetavec[k]<tvec[i]) lambdavec[i] += wvec[k];
   }
 }
 
 return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - IFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Integrated hazard rate function for IFR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_ifr_c(NumericVector tvec,
                         double lambda0,
                         NumericVector thetavec, 
                         NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector clambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   clambdavec[i] = lambda0*tvec[i];
   for(k=0; k<kmax; k++) {
     if(thetavec[k]<tvec[i]) clambdavec[i] += wvec[k]*(tvec[i]-thetavec[k]);
   }
 }
 
 return clambdavec;
}


//**********************************************************************
//* DFR - Decreasing hazard rate
//**********************************************************************
//' Hazard rate function - DFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Hazard rate function for DFR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_dfr_c(NumericVector tvec,
                        double lambda0, 
                        NumericVector thetavec, 
                        NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector lambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
   for(k=0; k<kmax; k++) {
     if(thetavec[k]>tvec[i]) lambdavec[i] += wvec[k];
   }
 }
 
 return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - DFR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Integrated hazard rate function for DFR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_dfr_c(NumericVector tvec,
                        double lambda0,
                        NumericVector thetavec, 
                        NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector clambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   clambdavec[i] = lambda0*tvec[i];
   for(k=0; k<kmax; k++) {
     if(thetavec[k]<tvec[i]) {
       clambdavec[i] += wvec[k]*thetavec[k]; 
     } else {
       clambdavec[i] += wvec[k]*tvec[i]; 
     }
   }
 }
 
 return clambdavec;
}


//**********************************************************************
//* CIR - Piecewise linear increasing hazard rate
//**********************************************************************
//' Hazard rate function - CIR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Hazard rate function for CIR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_cir_c(NumericVector tvec,
                        double lambda0, 
                        NumericVector thetavec, 
                        NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector lambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
   for(k=0; k<kmax; k++) {
     if(thetavec[k]<tvec[i]) lambdavec[i] += wvec[k]*(tvec[i]-thetavec[k]);
   }
 }
 
 return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - CIR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Integrated hazard rate function for CIR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_cir_c(NumericVector tvec,
                        double lambda0,
                        NumericVector thetavec, 
                        NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector clambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   clambdavec[i] = lambda0*tvec[i];
   for(k=0; k<kmax; k++) {
     if(thetavec[k]<tvec[i]) clambdavec[i] += 0.5*wvec[k]*pow(tvec[i]-thetavec[k],2);
   }
 }
 
 return clambdavec;
}




//**********************************************************************
//* CDR - Piecewise linear decreasing hazard rate
//**********************************************************************
//' Hazard rate function - CDR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Hazard rate function for CDR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_cdr_c(NumericVector tvec,
                        double lambda0, 
                        NumericVector thetavec, 
                        NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector lambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
   for(k=0; k<kmax; k++) {
     if(thetavec[k]>tvec[i]) lambdavec[i] += wvec[k]*(thetavec[k]-tvec[i]);
   }
 }
 
 return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - CDR
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Integrated hazard rate function for CIR hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_cdr_c(NumericVector tvec,
                        double lambda0,
                        NumericVector thetavec, 
                        NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector clambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   clambdavec[i] = lambda0*tvec[i];
   for(k=0; k<kmax; k++) {
     clambdavec[i] += 0.5*wvec[k]*pow(thetavec[k],2);
     if(thetavec[k]>tvec[i]) clambdavec[i] += -0.5*wvec[k]*pow(thetavec[k]-tvec[i],2);
   }
 }
 
 return clambdavec;
}


//**********************************************************************
//* LWB - Decreasing hazard rate
//**********************************************************************
//' Hazard rate function - LWB
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param a Change point
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Hazard rate function for LWB hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_lwb_c(NumericVector tvec,
                         double lambda0, 
                         double a, 
                         NumericVector thetavec, 
                         NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector lambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
   for(k=0; k<kmax; k++) {
     if(  (tvec[i]<a && tvec[i]<a-thetavec[k])
        ||(tvec[i]>a && tvec[i]>a+thetavec[k]) ) lambdavec[i] += wvec[k];
   }
 }
 
 return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - LWB
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param a Change point
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Integrated hazard rate function for LWB hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_lwb_c(NumericVector tvec,
                         double lambda0,
                         double a,
                         NumericVector thetavec, 
                         NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector clambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   clambdavec[i] = lambda0*tvec[i];
   for(k=0; k<kmax; k++) {
     if(tvec[i]<a) {
       if(thetavec[k]<a) {
          clambdavec[i] += wvec[k]*std::min(tvec[i],a-thetavec[k]);
       }
     } else {
       if(thetavec[k]<a) {
         clambdavec[i] += wvec[k]*(a-thetavec[k]);
       }
       if(thetavec[k]<tvec[i]-a) {
         clambdavec[i] += wvec[k]*(tvec[i]-a-thetavec[k]);
       }
     }
   }
 }
 
 return clambdavec;
}


//**********************************************************************
//* HBT - Ho Bathtub
//**********************************************************************
//' Hazard rate function - HBT
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param a Change point
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Hazard rate function for HBT hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_hbt_c(NumericVector tvec,
                         double lambda0, 
                         double a, 
                         NumericVector thetavec, 
                         NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector lambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   lambdavec[i] = lambda0;
   for(k=0; k<kmax; k++) {
     if(   (tvec[i]<a && tvec[i]<thetavec[k] && thetavec[k]<a)
        || (tvec[i]>a && a<thetavec[k] && thetavec[k]<tvec[i]) ) {
       lambdavec[i] += wvec[k]; 
     }
   }
 }
 
 return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - HBT
//' 
//' @param tvec Locations at which to evaluate the function
//' @param lambda0 Offset to be added to the hazard rate
//' @param a Change point
//' @param thetavec Locations
//' @param wvec Weights
//' 
//' @description Integrated hazard rate function for HBT hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_hbt_c(NumericVector tvec,
                         double lambda0,
                         double a,
                         NumericVector thetavec, 
                         NumericVector wvec
) {
 int n = tvec.size();
 int kmax = thetavec.size();
 NumericVector clambdavec(n);
 int i, k;
 
 for(i=0; i<n; i++) {
   clambdavec[i] = lambda0*tvec[i];
   for(k=0; k<kmax; k++) {
     if(tvec[i]<a) {
       if(thetavec[k]<a) {
          clambdavec[i] += wvec[k]*std::min(tvec[i],thetavec[k]);
       }
     } else {
       if(thetavec[k]<a) {
         clambdavec[i] += wvec[k]*thetavec[k]; 
       } else {
         clambdavec[i] += wvec[k]*std::max(tvec[i]-thetavec[k],0.0);
       }
     }
   }
 }
 
 return clambdavec;
}




//**********************************************************************
//* MEW - Modified Expontential Weibull
//**********************************************************************
//' Hazard rate function - MEW
//' 
//' @param tvec Locations at which to evaluate the function
//' @param alpha alpha
//' @param beta beta
//' @param nu nu
//' @param mu mu 
//' 
//' @description Hazard rate function for MEW hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector hazf_mew_c(NumericVector tvec,
                         double alpha, double beta,
                         double mu, double nu
) {
  int n = tvec.size();
  NumericVector lambdavec(n);
  double z, zd;
  int i;
  
  for(i=0; i<n; i++) {
    z = exp(pow(mu*tvec[i],beta))-1;
    zd = (1+z)*mu*beta*pow(mu*tvec[i],beta-1);
    lambdavec[i] = zd*(nu+alpha*pow(z,alpha-1));
  }
  
  return lambdavec;
}

//**********************************************************************
//' Integrated Hazard rate function - MEW
//' 
//' @param tvec Locations at which to evaluate the function
//' @param alpha alpha
//' @param beta beta
//' @param nu nu
//' @param mu mu 
//' 
//' @description Integrated hazard rate function for MEW hazard
//' 
//' @export
// [[Rcpp::export]]
NumericVector chzf_mew_c(NumericVector tvec,
                         double alpha, double beta,
                         double mu, double nu
) {
  int n = tvec.size();
  NumericVector clambdavec(n);
  double z;
  int i;
  
  for(i=0; i<n; i++) {
    z = exp(pow(mu*tvec[i],beta))-1;
    clambdavec[i] = nu*z + pow(z,alpha);
  }
  
  return clambdavec;
}


 
 