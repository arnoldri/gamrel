
#' Log likelihood of the current state
#'
llikef.old <- function(state, datlist, fpar, model) {
  # log likelihood of observations tvec
  lambda.vec <- lambda.func(tvec=datlist$tvec[datlist$obs],
                            model.list=c(list(model=model,kmax=fpar$kmax),
                                         state),
                            use.Cpp=fpar$use.Cpp)
  int.lambda.vec <- int.lambda.func(tvec=datlist$tvec,
                                    model.list=c(list(model=model,kmax=fpar$kmax),state),
                                    use.Cpp=fpar$use.Cpp)
  
  retval <- ( sum(log(lambda.vec)) - sum(int.lambda.vec) )
  return(retval)
}


#' Components of the log prior of the current state - return as a vector
#'
#' @export
lpriorf.vector.old <- function(state, fpar, model) {
  # log prior of the state
  lprior.vec <- rep(NA,length=length(fpar$parnames))
  names(lprior.vec) <- fpar$parnames
  if(model%in%c("IFR","DFR")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_ifrdfr_c(state$eta, 
                                    state$gamma,
                                    state$thetavec, 
                                    state$vvec,
                                    state$alpha, 
                                    state$beta,
                                    state$phi,
                                    fpar$nu, fpar$a1, fpar$a2, fpar$b1, fpar$b2, fpar$f1, fpar$f2)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for eta
      ## lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
      lprior.vec["eta"] <- -state$eta*fpar$nu 
      # prior for gamma
      lprior.vec["gamma"] <- dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
      # prior for thetavec
      lprior.vec["thetavec"] <- sum( (log(state$phi)-state$phi*state$thetavec) )
      # prior for vvec
      lprior.vec["vvec"] <-  sum( (log(state$alpha) + (state$alpha-1)*log(1-state$vvec[-fpar$kmax])) )
      # prior for alpha
      #lprior.vec["alpha"] <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha"] <- (fpar$a1-1)*log(state$alpha) - fpar$a2*state$alpha
      # prior for beta
      #lprior.vec["beta"] <- dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta"] <- (fpar$b1-1)*log(state$beta) - fpar$b2*state$beta
      # prior for phi
      #lprior.vec["phi"] <- dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
      lprior.vec["phi"] <- (fpar$f1-1)*log(state$phi) - fpar$f2*state$phi
    }
  } else if(model%in%c("LWB")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_lwb_c(state$eta, 
                                 state$a,
                                 state$gamma,
                                 state$thetavec, 
                                 state$vvec,
                                 state$alpha, 
                                 state$beta,
                                 state$phi,
                                 fpar$nu, fpar$c1, fpar$c2, 
                                 fpar$a1, fpar$a2, fpar$b1, fpar$b2, fpar$f1, fpar$f2)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for eta
      ## lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
      lprior.vec["eta"] <- -state$eta*fpar$nu 
      # prior for a
      ## lprior.vec["a"] <- dgamma(state$a, fpar$c1, fpar$c2)
      lprior.vec["a"] <- (fpar$c1-1)*log(state$a) - fpar$c2*state$a
      # prior for gamma
      lprior.vec["gamma"] <- dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
      # prior for thetavec
      lprior.vec["thetavec"] <- sum( (log(state$phi)-state$phi*state$thetavec) )
      # prior for vvec
      lprior.vec["vvec"] <-  sum( (log(state$alpha) + (state$alpha-1)*log(1-state$vvec[-fpar$kmax])) )
      # prior for alpha
      #lprior.vec["alpha"] <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha"] <- (fpar$a1-1)*log(state$alpha) - fpar$a2*state$alpha
      # prior for beta
      #lprior.vec["beta"] <- dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta"] <- (fpar$b1-1)*log(state$beta) - fpar$b2*state$beta
      # prior for phi
      #lprior.vec["phi"] <- dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
      lprior.vec["phi"] <- (fpar$f1-1)*log(state$phi) - fpar$f2*state$phi
    }
  } else if(model%in%c("SBT")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_sbt_c(state$eta, 
                                 state$gamma1,
                                 state$thetavec1, 
                                 state$vvec1,
                                 state$alpha1, 
                                 state$beta1,
                                 state$phi1,
                                 state$gamma2,
                                 state$thetavec2, 
                                 state$vvec2,
                                 state$alpha2, 
                                 state$beta2,
                                 state$phi2,
                                 fpar$nu, fpar$a1, fpar$a2, fpar$b1, fpar$b2, 
                                 fpar$f11, fpar$f21, fpar$f12, fpar$f22)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for eta
      ## lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
      lprior.vec["eta"] <- -state$eta*fpar$nu 
      
      # prior for gamma1
      lprior.vec["gamma1"] <- dgamma(state$gamma1, state$alpha1, state$beta1, log=TRUE)
      # prior for thetavec1
      lprior.vec["thetavec1"] <- sum( (log(state$phi1)-state$phi1*state$thetavec1) )
      # prior for vvec1
      lprior.vec["vvec1"] <-  sum( (log(state$alpha1) + (state$alpha1-1)*log(1-state$vvec1[-fpar$kmax])) )
      # prior for alpha1
      #lprior.vec["alpha1"] <- dgamma(state$alpha1, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha1"] <- (fpar$a1-1)*log(state$alpha1) - fpar$a2*state$alpha1
      # prior for beta1
      #lprior.vec["beta1"] <- dgamma(state$beta1, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta1"] <- (fpar$b1-1)*log(state$beta1) - fpar$b2*state$beta1
      # prior for phi1
      #lprior.vec["phi1"] <- dgamma(state$phi1, fpar$f11, fpar$f21, log=TRUE)
      lprior.vec["phi1"] <- (fpar$f11-1)*log(state$phi1) - fpar$f21*state$phi1
      
      # prior for gamma2
      lprior.vec["gamma2"] <- dgamma(state$gamma2, state$alpha2, state$beta2, log=TRUE)
      # prior for thetavec2
      lprior.vec["thetavec2"] <- sum( (log(state$phi2)-state$phi2*state$thetavec2) )
      # prior for vvec2
      lprior.vec["vvec2"] <-  sum( (log(state$alpha2) + (state$alpha2-1)*log(1-state$vvec2[-fpar$kmax])) )
      # prior for alpha2
      #lprior.vec["alpha2"] <- dgamma(state$alpha2, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha2"] <- (fpar$a1-1)*log(state$alpha2) - fpar$a2*state$alpha2
      # prior for beta2
      #lprior.vec["beta2"] <- dgamma(state$beta2, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta2"] <- (fpar$b1-1)*log(state$beta2) - fpar$b2*state$beta2
      # prior for phi2
      #lprior.vec["phi2"] <- dgamma(state$phi2, fpar$f12, fpar$f22, log=TRUE)
      lprior.vec["phi2"] <- (fpar$f12-1)*log(state$phi2) - fpar$f22*state$phi2
    }
  } else if(model%in%c("MBT")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_mbt_c(state$pival,
                                 state$eta1, 
                                 state$gamma1,
                                 state$thetavec1, 
                                 state$vvec1,
                                 state$alpha1, 
                                 state$beta1,
                                 state$phi1,
                                 state$eta2, 
                                 state$gamma2,
                                 state$thetavec2, 
                                 state$vvec2,
                                 state$alpha2, 
                                 state$beta2,
                                 state$phi2,
                                 fpar$nu, fpar$a1, fpar$a2, fpar$b1, fpar$b2, 
                                 fpar$f11, fpar$f21, fpar$f12, fpar$f22)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for pival
      lprior.vec["pival"] <- 0
      
      # prior for eta1
      lprior.vec["eta1"] <- 0 ##dexp(state$eta1, fpar$nu)
      # prior for gamma1
      lprior.vec["gamma1"] <- dgamma(state$gamma1, state$alpha1, state$beta1, log=TRUE)
      # prior for thetavec1
      lprior.vec["thetavec1"] <- sum( (log(state$phi1)-state$phi1*state$thetavec1) )
      # prior for vvec1
      lprior.vec["vvec1"] <-  sum( (log(state$alpha1) + (state$alpha1-1)*log(1-state$vvec1[-fpar$kmax])) )
      # prior for alpha1
      #lprior.vec["alpha1"] <- dgamma(state$alpha1, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha1"] <- (fpar$a1-1)*log(state$alpha1) - fpar$a2*state$alpha1
      # prior for beta1
      #lprior.vec["beta1"] <- dgamma(state$beta1, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta1"] <- (fpar$b1-1)*log(state$beta1) - fpar$b2*state$beta1
      # prior for phi1
      #lprior.vec["phi1"] <- dgamma(state$phi1, fpar$f11, fpar$f21, log=TRUE)
      lprior.vec["phi1"] <- (fpar$f11-1)*log(state$phi1) - fpar$f21*state$phi1
      
      # prior for eta2
      ## lprior.vec["eta2"] <- dexp(state$eta2, fpar$nu)
      lprior.vec["eta2"] <- -state$eta2*fpar$nu 
      # prior for gamma2
      lprior.vec["gamma2"] <- dgamma(state$gamma2, state$alpha2, state$beta2, log=TRUE)
      # prior for thetavec2
      lprior.vec["thetavec2"] <- sum( (log(state$phi2)-state$phi2*state$thetavec2) )
      # prior for vvec2
      lprior.vec["vvec2"] <-  sum( (log(state$alpha2) + (state$alpha2-1)*log(1-state$vvec2[-fpar$kmax])) )
      # prior for alpha2
      #lprior.vec["alpha2"] <- dgamma(state$alpha2, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha2"] <- (fpar$a1-1)*log(state$alpha2) - fpar$a2*state$alpha2
      # prior for beta2
      #lprior.vec["beta2"] <- dgamma(state$beta2, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta2"] <- (fpar$b1-1)*log(state$beta2) - fpar$b2*state$beta2
      # prior for phi2
      #lprior.vec["phi2"] <- dgamma(state$phi2, fpar$f12, fpar$f22, log=TRUE)
      lprior.vec["phi2"] <- (fpar$f12-1)*log(state$phi2) - fpar$f22*state$phi2
    }
  } else if(model%in%c("LCV")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_lcv_c(state$lambda0,
                                 state$w0,
                                 state$gamma,
                                 state$thetavec, 
                                 state$vvec,
                                 state$alpha, 
                                 state$beta,
                                 state$phi,
                                 fpar$s1, fpar$s2, fpar$sigmap.w0, 
                                 fpar$a1, fpar$a2, fpar$b1, fpar$b2, fpar$f1, fpar$f2)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for lambda0
      ## lprior.vec["lambda0"] <- dgamma(state$lambda0, fpar$s1, fpar$s2)
      lprior.vec["lambda0"] <- (fpar$s1-1)*log(state$lambda0) - fpar$s2*state$lambda0
      # prior for w0
      #lprior.vec["w0"] <- dnorm(state$w0, 0, fpar$sigmap.w0)
      lprior.vec["w0"] <- -0.5*(state$w0/fpar$sigmap.w0)^2
      # prior for gamma
      lprior.vec["gamma"] <- dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
      # prior for thetavec
      lprior.vec["thetavec"] <- sum( (log(state$phi)-state$phi*state$thetavec) )
      # prior for vvec
      lprior.vec["vvec"] <-  sum( (log(state$alpha) + (state$alpha-1)*log(1-state$vvec[-fpar$kmax])) )
      # prior for alpha
      #lprior.vec["alpha"] <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha"] <- (fpar$a1-1)*log(state$alpha) - fpar$a2*state$alpha
      # prior for beta
      #lprior.vec["beta"] <- dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta"] <- (fpar$b1-1)*log(state$beta) - fpar$b2*state$beta
      # prior for phi
      #lprior.vec["phi"] <- dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
      lprior.vec["phi"] <- (fpar$f1-1)*log(state$phi) - fpar$f2*state$phi
    }
  } else if(model%in%c("CIR","CDR")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_circdr_c(state$eta, 
                                    state$gamma,
                                    state$thetavec, 
                                    state$vvec,
                                    state$alpha, 
                                    state$beta,
                                    state$phi,
                                    fpar$nu, fpar$a1, fpar$a2, fpar$b1, fpar$b2, fpar$f1, fpar$f2)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for eta
      ## lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
      lprior.vec["eta"] <- -state$eta*fpar$nu 
      # prior for gamma
      lprior.vec["gamma"] <- dgamma(state$gamma, state$alpha, state$beta, log=TRUE)
      # prior for thetavec
      lprior.vec["thetavec"] <- sum( (log(state$phi)-state$phi*state$thetavec) )
      # prior for vvec
      lprior.vec["vvec"] <-  sum( (log(state$alpha) + (state$alpha-1)*log(1-state$vvec[-fpar$kmax])) )
      # prior for alpha
      #lprior.vec["alpha"] <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha"] <- (fpar$a1-1)*log(state$alpha) - fpar$a2*state$alpha
      # prior for beta
      #lprior.vec["beta"] <- dgamma(state$beta, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta"] <- (fpar$b1-1)*log(state$beta) - fpar$b2*state$beta
      # prior for phi
      #lprior.vec["phi"] <- dgamma(state$phi, fpar$f1, fpar$f2, log=TRUE)
      lprior.vec["phi"] <- (fpar$f1-1)*log(state$phi) - fpar$f2*state$phi
    }
  } else if(model%in%c("CVX")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_cvx_c(state$eta, 
                                 state$gamma1,
                                 state$thetavec1, 
                                 state$vvec1,
                                 state$alpha1, 
                                 state$beta1,
                                 state$phi1,
                                 state$gamma2,
                                 state$thetavec2, 
                                 state$vvec2,
                                 state$alpha2, 
                                 state$beta2,
                                 state$phi2,
                                 fpar$nu, 
                                 fpar$a1, fpar$a2, fpar$b1, fpar$b2, 
                                 fpar$f11, fpar$f21,
                                 fpar$f12, fpar$f22)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for eta
      ## lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
      lprior.vec["eta"] <- -state$eta*fpar$nu 
      
      # prior for gamma1
      lprior.vec["gamma1"] <- dgamma(state$gamma1, state$alpha1, state$beta1, log=TRUE)
      # prior for thetavec1
      lprior.vec["thetavec1"] <- sum( (log(state$phi1)-state$phi1*state$thetavec1) )
      # prior for vvec1
      lprior.vec["vvec1"] <-  sum( (log(state$alpha1) + (state$alpha1-1)*log(1-state$vvec1[-fpar$kmax])) )
      # prior for alpha1
      #lprior.vec["alpha1"] <- dgamma(state$alpha1, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha1"] <- (fpar$a1-1)*log(state$alpha1) - fpar$a2*state$alpha1
      # prior for beta1
      #lprior.vec["beta1"] <- dgamma(state$beta1, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta1"] <- (fpar$b1-1)*log(state$beta1) - fpar$b2*state$beta1
      # prior for phi1
      #lprior.vec["phi1"] <- dgamma(state$phi1, fpar$f11, fpar$f21, log=TRUE)
      lprior.vec["phi1"] <- (fpar$f11-1)*log(state$phi) - fpar$f21*state$phi
      
      # prior for gamma2
      lprior.vec["gamma2"] <- dgamma(state$gamma2, state$alpha2, state$beta2, log=TRUE)
      # prior for thetavec2
      lprior.vec["thetavec2"] <- sum( (log(state$phi2)-state$phi2*state$thetavec2) )
      # prior for vvec2
      lprior.vec["vvec2"] <-  sum( (log(state$alpha2) + (state$alpha2-1)*log(1-state$vvec2[-fpar$kmax])) )
      # prior for alpha2
      #lprior.vec["alpha2"] <- dgamma(state$alpha2, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha2"] <- (fpar$a1-1)*log(state$alpha2) - fpar$a2*state$alpha2
      # prior for beta2
      #lprior.vec["beta2"] <- dgamma(state$beta2, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta2"] <- (fpar$b1-1)*log(state$beta2) - fpar$b2*state$beta2
      # prior for phi2
      #lprior.vec["phi2"] <- dgamma(state$phi2, fpar$f12, fpar$f22, log=TRUE)
      lprior.vec["phi2"] <- (fpar$f12-1)*log(state$phi) - fpar$f22*state$phi
    }
  } else if(model%in%c("MEW")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_mew_c(state$lambda,
                                 state$alpha,
                                 state$theta,
                                 state$gamma, 
                                 fpar$s1, fpar$s2, 
                                 fpar$a1, fpar$a2, fpar$t1, fpar$t2, fpar$g1, fpar$g2)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for lambda
      ## lprior.vec["lambda"] <- dgamma(state$lambda, fpar$s1, fpar$s2)
      lprior.vec["lambda"] <- (fpar$s1-1)*log(state$lambda) - fpar$s2*state$lambda
      # prior for alpha
      #lprior.vec["alpha"] <- dgamma(state$alpha, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha"] <- (fpar$a1-1)*log(state$alpha) - fpar$a2*state$alpha
      # prior for theta
      #lprior.vec["theta"] <- dgamma(state$theta, fpar$t1, fpar$t2, log=TRUE)
      lprior.vec["theta"] <- (fpar$t1-1)*log(state$theta) - fpar$t2*state$theta
      # prior for gamma
      #lprior.vec["gamma"] <- dgamma(state$gamma, fpar$g1, fpar$g2, log=TRUE)
      lprior.vec["gamma"] <- (fpar$g1-1)*log(state$gamma) - fpar$g2*state$gamma
    }
    
  } else {
    stop("Specified model has not been implemented")
  }
  return(lprior.vec)
}
