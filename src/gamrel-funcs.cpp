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
//' Log prior - IFR/DFR
//' 
//' 
//' @param eta eta
//' @param gamma gamma
//' 
//' @description log(prior) for the IFR and DFR models - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_ifrdfr_c(double eta, 
                              double gamma,
                              NumericVector thetavec, 
                              NumericVector vvec,
                              double alpha, 
                              double beta,
                              double phi,
                              double nu, double a1, double a2, double b1, double b2, double f1, double f2) {
  
    int kmax = vvec.size();
    NumericVector lpriorvec(7);
    int i;
  
    // eta
    lpriorvec[0] = -eta*nu;
    // gamma
    lpriorvec[1] = R::dgamma(gamma, alpha, 1./beta, true);
    // thetavec
    lpriorvec[2] = kmax*log(phi);
    for(i=0; i<kmax; i++) lpriorvec[2] += (-phi*thetavec[i]);
    // vvec
    lpriorvec[3] =  (kmax-1)*log(alpha);
    for(i=0; i<kmax-1; i++) lpriorvec[3] += (alpha-1)*log(1-vvec[i]);
    // alpha
    lpriorvec[4] = (a1-1)*log(alpha) - a2*alpha;
    // beta
    lpriorvec[5] = (b1-1)*log(beta) - b2*beta;
    // phi
    lpriorvec[6] = (f1-1)*log(phi) - f2*phi;
    
    return lpriorvec;
}   

//' Log prior - LWB
//' 
//' 
//' @param eta eta
//' @param gamma gamma
//' 
//' @description log(prior) for the LWB model - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_lwb_c(double eta, 
                           double a,
                           double gamma,
                           NumericVector thetavec, 
                           NumericVector vvec,
                           double alpha, 
                           double beta,
                           double phi,
                           double nu, double c1, double c2, 
                           double a1, double a2, double b1, double b2, double f1, double f2) {
   
   int kmax = vvec.size();
   NumericVector lpriorvec(8);
   int i;
   
   // eta
   lpriorvec[0] = -eta*nu;
   // a 
   lpriorvec[1] = (c1-1)*log(a) - c2*a;
   // gamma
   lpriorvec[2] = R::dgamma(gamma, alpha, 1./beta, true);
   // thetavec
   lpriorvec[3] = kmax*log(phi);
   for(i=0; i<kmax; i++) lpriorvec[3] += (-phi*thetavec[i]);
   // vvec
   lpriorvec[4] =  (kmax-1)*log(alpha);
   for(i=0; i<kmax-1; i++) lpriorvec[4] += (alpha-1)*log(1-vvec[i]);
   // alpha
   lpriorvec[5] = (a1-1)*log(alpha) - a2*alpha;
   // beta
   lpriorvec[6] = (b1-1)*log(beta) - b2*beta;
   // phi
   lpriorvec[7] = (f1-1)*log(phi) - f2*phi;
   
   return lpriorvec;
}

//' Log prior - SBT
//' 
//' 
//' @param eta eta
//' @param gamma gamma
//' 
//' @description log(prior) for the SBT model - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_sbt_c(double eta, 
                            double gamma1,
                            NumericVector thetavec1, 
                            NumericVector vvec1,
                            double alpha1, 
                            double beta1,
                            double phi1,
                            double gamma2,
                            NumericVector thetavec2, 
                            NumericVector vvec2,
                            double alpha2, 
                            double beta2,
                            double phi2,
                            double nu, double a1, double a2, double b1, double b2, 
                            double f11, double f21, double f12, double f22) {
   
   int kmax = vvec1.size();
   NumericVector lpriorvec(13);
   int i;
   
   // eta
   lpriorvec[0] = -eta*nu;
   
   // gamma1
   lpriorvec[1] = R::dgamma(gamma1, alpha1, 1./beta1, true);
   // thetavec1
   lpriorvec[2] = kmax*log(phi1);
   for(i=0; i<kmax; i++) lpriorvec[2] += (-phi1*thetavec1[i]);
   // vvec1
   lpriorvec[3] =  (kmax-1)*log(alpha1);
   for(i=0; i<kmax-1; i++) lpriorvec[3] += (alpha1-1)*log(1-vvec1[i]);
   // alpha1
   lpriorvec[4] = (a1-1)*log(alpha1) - a2*alpha1;
   // beta1
   lpriorvec[5] = (b1-1)*log(beta1) - b2*beta1;
   // phi1
   lpriorvec[6] = (f11-1)*log(phi1) - f21*phi1;
   
   // gamma2
   lpriorvec[7] = R::dgamma(gamma2, alpha2, 1./beta2, true);
   // thetavec2
   lpriorvec[8] = kmax*log(phi2);
   for(i=0; i<kmax; i++) lpriorvec[8] += (-phi2*thetavec2[i]);
   // vvec2
   lpriorvec[9] =  (kmax-1)*log(alpha2);
   for(i=0; i<kmax-1; i++) lpriorvec[9] += (alpha2-1)*log(1-vvec2[i]);
   // alpha2
   lpriorvec[10] = (a1-1)*log(alpha2) - a2*alpha2;
   // beta2
   lpriorvec[11] = (b1-1)*log(beta2) - b2*beta2;
   // phi2
   lpriorvec[12] = (f12-1)*log(phi2) - f22*phi2;
   
   return lpriorvec;
}


//' Log prior - MBT
//' 
//' 
//' @param eta eta
//' @param gamma gamma
//' 
//' @description log(prior) for the MBT model - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_mbt_c(double pival, 
                           double eta1, 
                           double gamma1,
                           NumericVector thetavec1, 
                           NumericVector vvec1,
                           double alpha1, 
                           double beta1,
                           double phi1,
                           double eta2, 
                           double gamma2,
                           NumericVector thetavec2, 
                           NumericVector vvec2,
                           double alpha2, 
                           double beta2,
                           double phi2,
                           double nu, double a1, double a2, double b1, double b2, 
                           double f11, double f21, double f12, double f22) {
   
   int kmax = vvec1.size();
   NumericVector lpriorvec(15);
   int i;

   // pival
   lpriorvec[0] = 0;
   
   // eta1
   lpriorvec[1] = 0; // always zero-eta*nu;
   // gamma1
   lpriorvec[2] = R::dgamma(gamma1, alpha1, 1./beta1, true);
   // thetavec1
   lpriorvec[3] = kmax*log(phi1);
   for(i=0; i<kmax; i++) lpriorvec[3] += (-phi1*thetavec1[i]);
   // vvec1
   lpriorvec[4] =  (kmax-1)*log(alpha1);
   for(i=0; i<kmax-1; i++) lpriorvec[4] += (alpha1-1)*log(1-vvec1[i]);
   // alpha1
   lpriorvec[5] = (a1-1)*log(alpha1) - a2*alpha1;
   // beta1
   lpriorvec[6] = (b1-1)*log(beta1) - b2*beta1;
   // phi1
   lpriorvec[7] = (f11-1)*log(phi1) - f21*phi1;
   
   // eta2
   lpriorvec[8] = -eta2*nu;
   // gamma2
   lpriorvec[9] = R::dgamma(gamma2, alpha2, 1./beta2, true);
   // thetavec2
   lpriorvec[10] = kmax*log(phi2);
   for(i=0; i<kmax; i++) lpriorvec[10] += (-phi2*thetavec2[i]);
   // vvec2
   lpriorvec[11] =  (kmax-1)*log(alpha2);
   for(i=0; i<kmax-1; i++) lpriorvec[11] += (alpha2-1)*log(1-vvec2[i]);
   // alpha2
   lpriorvec[12] = (a1-1)*log(alpha2) - a2*alpha2;
   // beta2
   lpriorvec[13] = (b1-1)*log(beta2) - b2*beta2;
   // phi2
   lpriorvec[14] = (f12-1)*log(phi2) - f22*phi2;
   
   return lpriorvec;
}

//' Log prior - LCV
//' 
//' 
//' @param eta eta
//' @param gamma gamma
//' 
//' @description log(prior) for the LCV model - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_lcv_c(double lambda0,
                           double w0,
                           double gamma,
                           NumericVector thetavec, 
                           NumericVector vvec,
                           double alpha, 
                           double beta,
                           double phi,
                           double s1, double s2, double sigmap_w0, 
                           double a1, double a2, double b1, double b2, double f1, double f2) {
   
   int kmax = vvec.size();
   NumericVector lpriorvec(8);
   int i;

   // lambda0
   lpriorvec[0] = (s1-1)*log(lambda0) - s2*lambda0;
   // w0
   lpriorvec[1] = -0.5*pow(w0/sigmap_w0,2);
   // gamma
   lpriorvec[2] = R::dgamma(gamma, alpha, 1./beta, true);
   // thetavec
   lpriorvec[3] = kmax*log(phi);
   for(i=0; i<kmax; i++) lpriorvec[3] += (-phi*thetavec[i]);
   // vvec
   lpriorvec[4] =  (kmax-1)*log(alpha);
   for(i=0; i<kmax-1; i++) lpriorvec[4] += (alpha-1)*log(1-vvec[i]);
   // alpha
   lpriorvec[5] = (a1-1)*log(alpha) - a2*alpha;
   // beta
   lpriorvec[6] = (b1-1)*log(beta) - b2*beta;
   // phi
   lpriorvec[7] = (f1-1)*log(phi) - f2*phi;
   
   return lpriorvec;
}

//' Log prior - CVX
//' 
//' 
//' @param eta eta
//' @param gamma gamma
//' 
//' @description log(prior) for the CVX model - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_cvx_c(double eta, 
                            double tau,
                            double gamma1,
                            NumericVector thetavec1, 
                            NumericVector vvec1,
                            double alpha1, 
                            double beta1,
                            double phi1,
                            double gamma2,
                            NumericVector thetavec2, 
                            NumericVector vvec2,
                            double alpha2, 
                            double beta2,
                            double phi2,
                            double nu, 
                            double c1, double c2, 
                            double a1, double a2, double b1, double b2, 
                            double f1, double f2) {
   
   int kmax = vvec1.size();
   NumericVector lpriorvec(13);
   int i;
   
   // eta
   lpriorvec[0] = -eta*nu;
   // tau
   lpriorvec[1] = (c1-1)*log(tau) - c2*tau;
   
   // gamma1
   lpriorvec[2] = R::dgamma(gamma1, alpha1, 1./beta1, true);
   // thetavec1
   lpriorvec[3] = -kmax*log(tau);
   // vvec1
   lpriorvec[4] =  (kmax-1)*log(alpha1);
   for(i=0; i<kmax-1; i++) lpriorvec[3] += (alpha1-1)*log(1-vvec1[i]);
   // alpha1
   lpriorvec[5] = (a1-1)*log(alpha1) - a2*alpha1;
   // beta1
   lpriorvec[6] = (b1-1)*log(beta1) - b2*beta1;
   
   // gamma2
   lpriorvec[7] = R::dgamma(gamma2, alpha2, 1./beta2, true);
   // thetavec2
   lpriorvec[8] = kmax*log(phi2);
   for(i=0; i<kmax; i++) lpriorvec[8] += (-phi2*thetavec2[i]);
   // vvec2
   lpriorvec[9] =  (kmax-1)*log(alpha2);
   for(i=0; i<kmax-1; i++) lpriorvec[9] += (alpha2-1)*log(1-vvec2[i]);
   // alpha2
   lpriorvec[10] = (a1-1)*log(alpha2) - a2*alpha2;
   // beta2
   lpriorvec[11] = (b1-1)*log(beta2) - b2*beta2;
   // phi
   lpriorvec[12] = (f1-1)*log(phi2) - f2*phi2;
   
   return lpriorvec;
}


//' Log prior - MEW
//' 
//' 
//' @param lambda lambda
//' @param alpha alpha
//' @param theta theta
//' @param gamma gamma
//' 
//' @description log(prior) for the MEW model - in vector form
//' 
//' @export
// [[Rcpp::export]]
NumericVector lprior_mew_c(double lambda,
                           double alpha,
                           double theta,
                           double gamma, 
                           double s1, double s2, 
                           double a1, double a2, double t1, double t2, double g1, double g2) {
   
   NumericVector lpriorvec(4);

   // lambda
   lpriorvec[0] = (s1-1)*log(lambda) - s2*lambda;
   // alpha
   lpriorvec[1] = (a1-1)*log(alpha) - a2*alpha;
   // theta
   lpriorvec[2] = (t1-1)*log(theta) - t2*theta;
   // gamma
   lpriorvec[3] = (g1-1)*log(gamma) - g2*gamma;
   
   return lpriorvec;
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
  double s;
  int i,k;
  
  for(i=0; i<n; i++) {
    s = w0*tvec[i];
    for(k=0; k<kmax; k++) {
      if(thetavec[k]<=tvec[i]) s += wvec[k]*(tvec[i]-thetavec[k]);
    }
    lambdavec[i] = lambda0*exp(s);
  }
  
  return lambdavec;
}

//' Hazard rate function - CVX
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
NumericVector lambda_func_cvx_c(NumericVector tvec,
                                 double lambda0, 
                                 double tau,
                                 NumericVector thetavec1, 
                                 NumericVector wvec1,
                                 NumericVector thetavec2,
                                 NumericVector wvec2) {
   int n = tvec.size();
   int kmax = thetavec1.size();
   NumericVector lambdavec(n);
   int i,k;
   
   for(i=0; i<n; i++) {
     lambdavec[i] = lambda0;
     for(k=0; k<kmax; k++) {
       if(thetavec1[k]-tvec[i]>0) lambdavec[i] += wvec1[k]*(thetavec1[k]-tvec[i]);
     }
     for(k=0; k<kmax; k++) {
       if(tvec[i]-thetavec2[k]>0) lambdavec[i] += wvec2[k]*(tvec[i]-thetavec1[k]);
     }
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
    k1 = 0;
    while( (k1<kmaxp2-1) && (othetavec[k1+1]<=tvec[i]) ) {   // try < rather than <=?
      k1++;
    }
    bothlambdamat(i,0) = lambda0*exp(s01vec[k1]*tvec[i]-s2vec[k1]);
    bothlambdamat(i,1) = ssvec[k1];
    if(std::abs(s01vec[k1])<epsilon) {
      bothlambdamat(i,1) += ccvec[k1]*(tvec[i]-othetavec[k1]);
    } else {
      bothlambdamat(i,1) += (ccvec[k1]/s01vec[k1])*(exp(s01vec[k1]*tvec[i])-exp(s01vec[k1]*othetavec[k1]));
    }
  }
  
  return bothlambdamat;
}



