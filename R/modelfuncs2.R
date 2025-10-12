####################################################
# Utils
unitquadratic <- function(u,a) {
  # -1 < a < 1 for function increasing at both 0 and 1
  # may have vector u and scalar a, or both vectors of the same length
  a*u^2 + (1-a)*u
}
solveunitquadratic <- function(d,a) {
  # -1 < a < 1 for function increasing at both 0 and 1
  # return solution in 0<u<1, assuming 0<d<1
  # may have vector u and scalar a, or both vectors of the same length
  a <- rep_len(a,length.out=length(d))
  u <- ifelse(a==0,d,
              -(1/(2.*a))*( 1-a - sqrt( (1-a)^2+4*a*d) ))
  return(u)
}
ptquadratic <- function(t,t1,t2,f1,f2,f1dash) {
  # quadratic passing through (t1,f1) and (t2,f2)
  # with gradient f1 at t1
  dfdt <- (f2-f1)/(t2-t1)
  a <- ifelse(f1==f2 | dfdt==f1dash,0,
              1 - f1dash/dfdt)
  f1 + (f2-f1)*unitquadratic((t-t1)/(t2-t1),a)
}
solveptquadratic <- function(f,t1,t2,f1,f2,f1dash) {
  # find t where f(t)=f
  dfdt <- (f2-f1)/(t2-t1)
  a <- ifelse(f1==f2 | dfdt==f1dash,0,
              1 - f1dash/dfdt)
  t1 + (t2-t1)*solveunitquadratic((f-f1)/(f2-f1),a)
}
ptlinear <- function(t,t1,t2,f1,f2) {
  # linear function passing through (t1,f1) and (t2,f2)
  dfdt <- (f2-f1)/(t2-t1)
  f1 + (t-t1)*dfdt
}
solveptlinear <- function(f,t1,t2,f1,f2) {
  # find t where f(t)=f
  dfdt <- (f2-f1)/(t2-t1)
  t1 + (f-f1)/dfdt
}

####################################################
#' Make data list object
#' 
#' @param tvec Vector of failure times
#' @param obs Vector of TRUE (observed) or FALSE (right censored) indicators
#' 
#' @export
make.datlist <- function(tvec, obs=NULL) {
  n <- length(tvec)
  if(is.null(obs)) {
    obs <- rep(TRUE,n)
  } else if(length(obs==1)) {
    obs <- rep(obs,n)
  } else if(length(obs)!=n) {
    stop("obs must be the same length as tvec")
  }
  datlist <- list(n=n, nobs=sum(obs), tvec=tvec, obs=obs)
  return(datlist)
}

####################################################
# Generic
# Models are CON/ IFR/DFR/ CIR/CDR/ LWB/ HBT/HCV/ SBT/SCV/ MBT/LCV/ MEW

#' Hazard function
#' 
#' @param tvec Vector of failure times
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' 
#' @export
hazf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=hazf.CON(tvec, model.list, use.Cpp),
         IFR=hazf.IFR(tvec, model.list, use.Cpp),
         DFR=hazf.DFR(tvec, model.list, use.Cpp),
         CIR=hazf.CIR(tvec, model.list, use.Cpp),
         CDR=hazf.CDR(tvec, model.list, use.Cpp),
         LWB=hazf.LWB(tvec, model.list, use.Cpp),
         HBT=hazf.HBT(tvec, model.list, use.Cpp),
         HCV=hazf.HCV(tvec, model.list, use.Cpp),
         SBT=hazf.SBT(tvec, model.list, use.Cpp),
         SCV=hazf.SCV(tvec, model.list, use.Cpp),
         MBT=hazf.MBT(tvec, model.list, use.Cpp),
         LCV=hazf.LCV(tvec, model.list, use.Cpp),
         MEW=hazf.MEW(tvec, model.list, use.Cpp)
  )
}

#' Cumulative hazard function
#' 
#' @param tvec Vector of failure times
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' 
#' @export
chzf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=chzf.CON(tvec, model.list, use.Cpp),
         IFR=chzf.IFR(tvec, model.list, use.Cpp),
         DFR=chzf.DFR(tvec, model.list, use.Cpp),
         CIR=chzf.CIR(tvec, model.list, use.Cpp),
         CDR=chzf.CDR(tvec, model.list, use.Cpp),
         LWB=chzf.LWB(tvec, model.list, use.Cpp),
         HBT=chzf.HBT(tvec, model.list, use.Cpp),
         HCV=chzf.HCV(tvec, model.list, use.Cpp),
         SBT=chzf.SBT(tvec, model.list, use.Cpp),
         SCV=chzf.SCV(tvec, model.list, use.Cpp),
         MBT=chzf.MBT(tvec, model.list, use.Cpp),
         LCV=chzf.LCV(tvec, model.list, use.Cpp),
         MEW=chzf.MEW(tvec, model.list, use.Cpp)
  )
}

#' Survival function
#' 
#' @param tvec Vector of failure times
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' @param log Return log of survival function?
#' 
#' @export
survf <- function(tvec, model.list, use.Cpp=FALSE, log=FALSE) {
  survf <- -chzf(tvec, model.list, use.Cpp)
  if(!log) survf <- exp(survf)
  return(survf)
}

#' Inverse Survival function
#' 
#' @param uvec Vector of survival function values
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' 
#' @export
invsurvf <- function(uvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=invsurvf.CON(uvec, model.list, use.Cpp),
         IFR=invsurvf.IFR(uvec, model.list, use.Cpp),
         DFR=invsurvf.DFR(uvec, model.list, use.Cpp),
         CIR=invsurvf.CIR(uvec, model.list, use.Cpp),
         CDR=invsurvf.CDR(uvec, model.list, use.Cpp),
         LWB=invsurvf.LWB(uvec, model.list, use.Cpp),
         HBT=invsurvf.HBT(uvec, model.list, use.Cpp),
         HCV=invsurvf.HCV(uvec, model.list, use.Cpp),
         SBT=invsurvf.SBT(uvec, model.list, use.Cpp),
         SCV=invsurvf.SCV(uvec, model.list, use.Cpp),
         MBT=invsurvf.MBT(uvec, model.list, use.Cpp),
         LCV=invsurvf.LCV(uvec, model.list, use.Cpp),
         MEW=invsurvf.MEW(uvec, model.list, use.Cpp)
  )
}


#' Hazard and cumulative hazard function values
#' 
#' @param tvec Vector of failure times
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' @param epsilon Small nonzero value
#' 
#' @export
hazf_chzf <- function(tvec, model.list, use.Cpp, 
                      epsilon=100*.Machine$double.neg.eps) {
  if(use.Cpp) {
    retval <- switch(model.list$model,
         CON=hazf_chzf_con_c(tvec, model.list$lambda0),
         IFR=hazf_chzf_ifr_c(tvec, model.list$lambda0,
                             model.list$thetavec, model.list$wvec),
         DFR=hazf_chzf_dfr_c(tvec, model.list$lambda0,
                             model.list$thetavec, model.list$wvec),
         CIR=hazf_chzf_cir_c(tvec, model.list$lambda0,
                             model.list$thetavec, model.list$wvec),
         CDR=hazf_chzf_cdr_c(tvec, model.list$lambda0,
                             model.list$thetavec, model.list$wvec),
         LWB=hazf_chzf_lwb_c(tvec, model.list$lambda0, model.list$a, 
                             model.list$thetavec, model.list$wvec),
         HBT=hazf_chzf_hbt_c(tvec, model.list$lambda0, model.list$a, 
                             model.list$thetavec, model.list$wvec),
         HCV=hazf_chzf_hcv_c(tvec, model.list$lambda0, model.list$a, 
                             model.list$thetavec, model.list$wvec),
         SBT=hazf_chzf_sbt_c(tvec, model.list$lambda0,
                             model.list$thetavec1, model.list$wvec1,
                             model.list$thetavec2, model.list$wvec2),
         SCV=hazf_chzf_scv_c(tvec, model.list$lambda0,
                             model.list$thetavec1, model.list$wvec1,
                             model.list$thetavec2, model.list$wvec2),
         MBT=hazf_chzf_mbt_c(tvec, model.list$pival, model.list$lambda0,
                             model.list$thetavec1, model.list$wvec1,
                             model.list$thetavec2, model.list$wvec2),
         LCV=hazf_chzf_lcv_c(tvec, model.list$lambda0, model.list$w0,
                             model.list$thetavec, model.list$wvec,
                             epsilon),
         MEW=hazf_chzf_mew_c(tvec, model.list$alpha, model.list$beta,
                             model.list$mu, model.list$nu))
  } else {
    retval <- cbind(hazf(tvec, model.list, use.Cpp),
                    chzf(tvec, model.list, use.Cpp))
  } 
  return(retval)
}

#' Probability density function
#' 
#' @param tvec Vector of failure times
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' 
#' @export
densf <- function(tvec, model.list, use.Cpp=FALSE, log=FALSE) {
  if(log) {
     retval <- log(hazf(tvec, model.list, use.Cpp)) - chzf(tvec, model.list, use.Cpp)
  } else {
     retval <- hazf(tvec, model.list, use.Cpp)*exp(-chzf(tvec, model.list, use.Cpp))
  }
  return(retval)
}

#' Generate failure times from a model
#' 
#' @param n Number of failure times to simulate
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?  TRUE or FALSE
#' @param seed Random number seed
#' 
#' @export
rgen <- function(n, model.list, use.Cpp=FALSE, seed=NULL) {
  if(model.list$model=="MBT") { # special case of the Mixture Bathtub
    tvec <- rgen.MBT(n, model.list, use.Cpp=use.Cpp, seed=seed)
  } else {
    if(!is.null(seed)) set.seed(seed)
    uvec <- runif(n)
    tvec <- invsurvf(uvec, model.list=model.list, use.Cpp=use.Cpp)
  }
  return(tvec)
}

#' Log likelihood function
#' 
#' @param tvec Vector of failure times
#' @param parvec Parameter vector
#' @param model Model code (three letters)
#' @param verbose Display output in iterative evaluations?  TRUE or FALSE
#' 
#' @export
llfunc <- function(datlist, model.list, verbose=FALSE) {
  retval <- sum(densf(datlist$tvec[datlist$obs], model.list, use.Cpp=TRUE, log=TRUE))
  retval <- retval + sum(densf(datlist$tvec[!datlist$obs], model.list, use.Cpp=TRUE, log=TRUE))
  return(retval)
}
#' Log likelihood of the current state
#'
#' @export
llikef <- function(state, datlist, fpar, model) {
  # log likelihood of observations datlist
  #hazfuncs <- both.lambda.funcs(datlist$tvec,
  #                              model.list=c(list(model=model,kmax=fpar$kmax),
  #                                           state),
  #                              use.Cpp=fpar$use.Cpp, epsilon=fpar$epsilon)
  hazfuncs <- hazf_chzf(datlist$tvec, 
                        model.list=c(list(model=model,kmax=fpar$kmax),
                                          state), 
                        use.Cpp=fpar$use.Cpp, 
                        epsilon=100*.Machine$double.neg.eps) 
  if(model=="LCV" && all(!is.nan(hazfuncs[,1])) 
     && any(is.nan(hazfuncs[,2]))) {
    # Catch the case where the LCV integrated hazard function overflows
    cat("*")
    retval <- -Inf
  } else {
    retval <- sum(log(hazfuncs[datlist$obs,1]) - hazfuncs[,2])
  }
  return(retval)
}

#' Log prior of the current state - scalar valued
#'
#' @export
lpriorf <- function(state, fpar, model, use.Cpp=TRUE) {
  # log prior of the state
  retval <- sum(lpriorfunc.vector(state, fpar, model, use.Cpp=use.Cpp))
  return(retval)
}

#' Log prior function - vector valued
#' 
#' @param model.list Model specification
#' @param use.Cpp Use C++ functions?
#' 
#' @export
lpriorfunc.vector <- function(state, fpar, model, use.Cpp=TRUE) {
    retval <- switch(model, 
         CON=logprior.CON(state, fpar, use.Cpp),
         IFR=logprior.IFR(state, fpar, use.Cpp),
         DFR=logprior.DFR(state, fpar, use.Cpp),
         CIR=logprior.CIR(state, fpar, use.Cpp),
         CDR=logprior.CDR(state, fpar, use.Cpp),
         LWB=logprior.LWB(state, fpar, use.Cpp),
         HBT=logprior.HBT(state, fpar, use.Cpp),
         HCV=logprior.HCV(state, fpar, use.Cpp),
         SBT=logprior.SBT(state, fpar, use.Cpp),
         SCV=logprior.SCV(state, fpar, use.Cpp),
         MBT=logprior.MBT(state, fpar, use.Cpp),
         LCV=logprior.LCV(state, fpar, use.Cpp),
         MEW=logprior.MEW(state, fpar, use.Cpp)
  )
  return(retval)
}

#' Convert model list object to parameter vector (incomplete)
#' 
#' @param model.list Model specification
#' 
#' @export
as.parvec <- function(model.list) { ##!!!== incomplete
  if(model.list$model=="CON") {
    parvec <- log(model.list$lambda0)
  } else if(model.list$model=="MEW") {
    parvec <- log(c(model.list$alpha, model.list$beta, model.list$mu, model.list$nu))
  } else {
    stop(paste("Model",model.list$model,"not recognised"))
  }
}
#' Convert parameter vector to model list object (incomplete)
#' 
#' @param model.list Model specification
#' 
#' @export
as.model.list <- function(parvec,model) {  ##!!== incomplete
  if(model=="CON") {
    model.list <- list(model=model, lambda0=exp(parvec[1]))
  } else if(model=="MEW") {
    model.list <- list(model=model, 
                       alpha=exp(parvec[1]), beta=exp(parvec[2]), 
                       mu=exp(parvec[3]), nu=exp(parvec[4]))
  } else {
    stop(paste("Model",model,"not recognised"))
  }
}

#' Plot a hazard function
#' 
#' @export
plot.hazf <- function(model.list, xlim, use.Cpp=FALSE, n=101, add=FALSE, 
                      xlab="t", ylab=expression(lambda(t)), scale=1, ...) {
  tvec <- seq(from=xlim[1], to=xlim[2], length=n)
  fvec <- scale*hazf(tvec,model.list,use.Cpp) 
  ylim <- c(0,1.1*max(fvec))
  if(!add) {
    plot(tvec, fvec, ylim=ylim, xlab=xlab, ylab=ylab, type="l", ...)
  } else {
    lines(tvec, fvec, ...)
  }
  if(!add) mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

#' Plot a cumulative hazard function
#' 
#' @export
plot.chzf <- function(model.list, xlim, use.Cpp=FALSE, n=101, add=FALSE, 
                      xlab="t", ylab=expression(Lambda(t)), scale=1, ...) {
  curve(scale*chzf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  if(!add) mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

#' Plot a survival function
#' 
#' @export
plot.survf <- function(model.list, xlim, ylim=c(0,1), use.Cpp=FALSE, n=101, add=FALSE, 
                       xlab="t", ylab=expression(bar(F)(t)), scale=1, ...) {
  curve(scale*survf(x,model.list,use.Cpp), xlim=xlim, ylim=ylim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  if(!add) mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

#' Plot a density function
#' 
#' @export
plot.densf <- function(model.list, xlim, use.Cpp=FALSE, n=101, add=FALSE, 
                       xlab="t", ylab=expression(f(t)), scale=1, ...) {
  curve(scale*densf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  if(!add) mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

#' Plot an inverse survival function
#' 
#' @export
plot.invsurvf <- function(model.list, xlim=c(0,1)+0.0001*c(1,-1), use.Cpp=FALSE, n=101, add=FALSE, 
                          xlab=expression(bar(F)(t)), ylab="t", scale=1, ...) {
  curve(scale*invsurvf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  if(!add) mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

#' Maximum likelihood fitting for a restricted set of models
#' 
#' @param model Chosen model
#' @param tvec Vector of failure times
#' @param obs Vector of logical values indicating observation or right-censoring
#' 
#' @export
mlfit <- function(model, tvec, obs=TRUE) {
   if(model=="CON") {
     retval <- mlfit.CON(tvec, obs)
   } else {
     stop("Model not recognised")
   }
  return(retval)
}

####################################################
# CON
# lambda(t) = lambda0
# model.list = list(model="CON", lambda0)
hazf.CON <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_con_c(tvec, model.list$lambda0)
  } else {
    lambdavec <- rep(model.list$lambda0, length=length(tvec))
  }
  return(lambdavec)
}
chzf.CON <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_con_c(tvec, model.list$lambda0)
  } else {
    clambdavec <- model.list$lambda0*tvec
  }
  return(clambdavec)
}
invsurvf.CON <- function(uvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    tvec <- invsurvf_con_c(uvec, model.list$lambda0)
  } else {
    tvec <- -log(uvec)/model.list$lambda0
  }
  return(tvec)
}
logprior.CON <- function(state, fpar, use.Cpp) {
  if(use.Cpp) {
    retval <- logprior_con_c(state$lambda0, fpar$nu)
  } else {
    retval <- dexp(state$lambda0, fpar$nu, log=TRUE)
  }
  return(retval)
}
mlfit.CON <- function(tvec, obs=TRUE) {
  datlist <- make.datlist(tvec, obs)
  ff <- function(par, tvec, obs, verbose=FALSE) {
    lambda0 <- exp(par)
    n0 <- sum(obs)
    retval <- n0*par - sum(lambda0*tvec)
  }
  lambda0 <- 3./mean(tvec[datlist$obs])
  par0 <- log(lambda0)
  opt <- optim(par, ff, 
               lower=log(1/max(tvec)), upper=log(1/min(tvec)),
               hessian=TRUE, method="Brent",
               control=list(fnscale=-1),
               tvec=datlist$tvec, obs=datlist$obs)
  lambda0 <- exp(opt$par)
  model.list <- list(model="CON", lambda0=lambda0)
  se.lambda0 <- as.vector(lambda0/sqrt(-opt$hessian))
  return(c(model.list,list(se.lambda0=se.lambda0)))
}

####################################################
# IFR
# lambda(t) = lambda0 + int_0^t G(du)
# model.list = list(model="IFR", kmax, lambda0, thetavec, wvec)
hazf.IFR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_ifr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- model.list$lambda0 + as.vector(outer(tvec,model.list$thetavec,">")%*%model.list$wvec)
  }
  return(lambdavec)
}
chzf.IFR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_ifr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- ( (model.list$lambda0 + sum(model.list$wvec))*tvec 
                    - as.vector(outer(tvec,model.list$thetavec,pmin)%*%model.list$wvec)
                  )
    clambdavec <- pmax(0,clambdavec)
  }
  return(clambdavec)
}
invsurvf.IFR <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  odx <- order(model.list$thetavec)
  # integrated hazard at these values
  ivec <- c(0,chzf(model.list$thetavec[odx],model.list,use.Cpp))
  # cumulative sum of weights
  cvec <- c(0,cumsum(model.list$wvec[odx]))
  # cumulative sum of weights*theta
  dvec <- c(0,cumsum((model.list$wvec*model.list$thetavec)[odx]))
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # linear interpolation
  tvec <- (dvec[kvec]-log(uvec))/(model.list$lambda0+cvec[kvec])
  return(tvec)
}

####################################################
# DFR
# lambda(t) = lambda0 + int_t^infty G(du)
# model.list = list(model="DFR", kmax, lambda0, thetavec, wvec)
hazf.DFR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_dfr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- model.list$lambda0 + as.vector(outer(tvec,model.list$thetavec,"<")%*%model.list$wvec)
  }
  return(lambdavec)
}
chzf.DFR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_dfr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- model.list$lambda0*tvec + as.vector(outer(tvec,model.list$thetavec,pmin)%*%model.list$wvec)
  }
  return(clambdavec)
}
invsurvf.DFR <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  othetavec <- sort(model.list$thetavec)
  othetavec <- c(0,othetavec,1.2*othetavec[length(othetavec)])
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  kvecp1 <- kvec+1
  # linear interpolation
  tvec <- ( othetavec[kvec]
            + (othetavec[kvecp1]-othetavec[kvec])/(ivec[kvecp1]-ivec[kvec])*(-log(uvec)-ivec[kvec])
  )
  return(tvec)
}

####################################################
# CIR
# lambda(t) = lambda0 + int_0^t int_0^u G(dv) du
# model.list = list(model="CIR", kmax, lambda0, thetavec, wvec)
hazf.CIR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_cir_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- model.list$lambda0 + as.vector(outer(tvec,model.list$thetavec,function(t,theta) pmax(0,t-theta))%*%model.list$wvec)
  }
  return(lambdavec)
}
chzf.CIR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_cir_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- ( model.list$lambda0*tvec
                    + 0.5*as.vector(outer(tvec,model.list$thetavec,function(t,theta) (pmax(0,t-theta))^2)%*%model.list$wvec)
    )
  }
  return(clambdavec)
}
invsurvf.CIR <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  odx <- order(model.list$thetavec)
  othetavec <- model.list$thetavec[odx]
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # hazard at these values
  hvec <- hazf(othetavec,model.list,use.Cpp)
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # quadratic interpolation
  tvec <- solveptquadratic(-log(uvec),
                           othetavec[kvec],othetavec[kvec+1],
                           ivec[kvec],ivec[kvec+1],
                           hvec[kvec]) 
  return(tvec)
}


####################################################
# CDR
# lambda(t) = lambda0 + int_t^infty int_u^infty G(dv) du
# model.list = list(model="CIR", kmax, lambda0, thetavec, wvec)
hazf.CDR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_cdr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- model.list$lambda0 + as.vector(outer(tvec,model.list$thetavec,function(t,theta) pmax(0,theta-t))%*%model.list$wvec)
  }
  return(lambdavec)
}
chzf.CDR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_cdr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- ( model.list$lambda0*tvec
                    + 0.5*sum(model.list$wvec*model.list$thetavec^2)
                    - 0.5*as.vector(outer(tvec,model.list$thetavec,function(t,theta) (pmax(0,theta-t))^2)%*%model.list$wvec)
    )
  }
  return(clambdavec)
}
invsurvf.CDR <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  odx <- order(model.list$thetavec)
  othetavec <- model.list$thetavec[odx]
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # hazard at these values
  hvec <- hazf(othetavec,model.list,use.Cpp)
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # quadratic interpolation
  tvec <- solveptquadratic(-log(uvec),
                           othetavec[kvec],othetavec[kvec+1],
                           ivec[kvec],ivec[kvec+1],
                           hvec[kvec]) 
  return(tvec)
}

####################################################
# LWB
# lambda(t) = lambda0 + int_0^|t-a| G(du)
# model.list = list(model="LWB", kmax, lambda0, a, thetavec, wvec)
hazf.LWB <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_lwb_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- model.list$lambda0 + ifelse(tvec<model.list$a,
                                  as.vector(outer(tvec,model.list$a-model.list$thetavec,"<")%*%model.list$wvec),
                                  as.vector(outer(tvec,model.list$a+model.list$thetavec,">")%*%model.list$wvec))
  }
  return(lambdavec)
}
chzf.LWB <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_lwb_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- model.list$lambda0*tvec
    clambdavec <- clambdavec + as.vector(outer(pmin(model.list$a,tvec),(model.list$a-model.list$thetavec),pmin)
                                         %*%
                                         ((model.list$thetavec<model.list$a)*model.list$wvec)
                                        )
    clambdavec <- clambdavec + ifelse(tvec>model.list$a,
                                      as.vector(outer(tvec,model.list$thetavec,
                                                      function(t,theta) (t-model.list$a-theta)*(theta<t-a))
                                                %*%model.list$wvec),
                                      0)
  }
  return(clambdavec)
}
invsurvf.LWB <- function(uvec, model.list, use.Cpp=FALSE) {
  # knots at which the function changes behaviour
  othetavec <- model.list$a + c(0,model.list$thetavec,-model.list$thetavec)
  othetavec <- othetavec[othetavec>0]
  othetavec <- sort(othetavec)
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # linear interpolation
  tvec <- solveptlinear(-log(uvec),
                        othetavec[kvec],othetavec[kvec+1],
                        ivec[kvec],ivec[kvec+1])
  return(tvec)
}

####################################################
# HBT
# lambda(t) = lambda0 + int_min(a,t)^a G(du) + int_a^max(a,t) G(du)
# model.list = list(model="CIR", kmax, lambda0, a, thetavec, wvec)
hazf.HBT <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_hbt_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- (model.list$lambda0 + 
                            ifelse(tvec<model.list$a,
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) (t<theta & theta<model.list$a))%*%model.list$wvec),
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) (model.list$a<theta & theta<t))%*%model.list$wvec)
                                   ) )
  }
  return(lambdavec)
}
chzf.HBT <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_hbt_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- model.list$lambda0*tvec
    clambdavec <- clambdavec + ifelse(tvec<model.list$a,
                                      as.vector(outer(tvec,model.list$thetavec,function(t,theta) ((theta<model.list$a)*pmin(t,theta)))%*%model.list$wvec),
                                      as.vector(outer(tvec,model.list$thetavec,function(t,theta) ((theta<model.list$a)*pmin(model.list$a,theta)))%*%model.list$wvec)
                                      +as.vector(outer(tvec,model.list$thetavec,function(t,theta) ((theta>model.list$a)*pmax(t-theta,0)))%*%model.list$wvec)
    )
  }
  return(clambdavec)
}
invsurvf.HBT <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  othetavec <- c(model.list$a, model.list$thetavec)
  othetavec <- sort(othetavec)
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # quadratic interpolation
  tvec <- solveptlinear(-log(uvec),
                        othetavec[kvec],othetavec[kvec+1],
                        ivec[kvec],ivec[kvec+1])
  return(tvec)
}

####################################################
# HCV
# lambda(t) = lambda0 + int_t^infty int_min(a,u)^a G(dv)du + int_0^t int_a^max(a,u) G(dv)dt
# model.list = list(model="HCV", kmax, lambda0, a, thetavec, wvec)
hazf.HCV <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_hcv_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- ( model.list$lambda0 + 
                           ifelse(tvec<model.list$a,
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) 
                                    (t<theta & theta<model.list$a)*pmax(theta-t,0))%*%model.list$wvec),
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) 
                                    (model.list$a<theta & theta<t)*pmax(t-theta,0))%*%model.list$wvec)
                                  ) )
  }
  return(lambdavec)
}
chzf.HCV <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_hcv_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- model.list$lambda0*tvec
    clambdavec <- clambdavec + ifelse(tvec<model.list$a,
                                      as.vector(outer(tvec,model.list$thetavec,function(t,theta) 
                                        ((theta<model.list$a)*0.5*(theta^2-(pmax(0,theta-t))^2)))%*%model.list$wvec),
                                      as.vector(outer(tvec,model.list$thetavec,function(t,theta) 
                                        ((theta<model.list$a)*0.5*theta^2))%*%model.list$wvec)
                                      +as.vector(outer(tvec,model.list$thetavec,function(t,theta) 
                                        ((model.list$a<theta)*(theta<t)*0.5*(t-theta)^2))%*%model.list$wvec)
                                      )
  }
  return(clambdavec)
}
invsurvf.HCV <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  othetavec <- c(model.list$a, model.list$thetavec)
  othetavec <- sort(othetavec)
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # hazard at these values
  hvec <- hazf(othetavec,model.list,use.Cpp)
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # quadratic interpolation
  tvec <- solveptquadratic(-log(uvec),
                           othetavec[kvec],othetavec[kvec+1],
                           ivec[kvec],ivec[kvec+1],
                           hvec[kvec]) 
  return(tvec)
}

####################################################
# SBT
# lambda(t) = lambda0 + int_t^infty G1(du) + int_0^t G2(du) = DFR+IFR
# model.list = list(model="SBT", kmax, lambda0, thetavec1, wvec1, thetavec2, wvec2)
hazf.SBT <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_sbt_c(tvec, model.list$lambda0, 
                            model.list$thetavec1, model.list$wvec1,
                            model.list$thetavec2, model.list$wvec2)
  } else {
    lambdavec <- ( model.list$lambda0 
         + as.vector(outer(tvec,model.list$thetavec1,"<")%*%model.list$wvec1)
         + as.vector(outer(tvec,model.list$thetavec2,">")%*%model.list$wvec2) )
  }
  return(lambdavec)
}
chzf.SBT <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_sbt_c(tvec, model.list$lambda0, 
                             model.list$thetavec1, model.list$wvec1,
                             model.list$thetavec2, model.list$wvec2)
  } else {
    clambdavec <- ( model.list$lambda0*tvec
                    + as.vector(outer(tvec,model.list$thetavec1,pmin)%*%model.list$wvec1)
                    + sum(model.list$wvec2)*tvec 
                      - as.vector(outer(tvec,model.list$thetavec2,pmin)%*%model.list$wvec2)
                  )
  }
  return(clambdavec)
}
invsurvf.SBT <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  othetavec <- c(model.list$thetavec1, model.list$thetavec2)
  othetavec <- sort(othetavec)
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # quadratic interpolation
  tvec <- solveptlinear(-log(uvec),
                        othetavec[kvec],othetavec[kvec+1],
                        ivec[kvec],ivec[kvec+1])
  return(tvec)
}

####################################################
# SCV
# lambda(t) = lambda0 + int_t^infty int_u^infty G1(dv)du + int_0^t int_0^u G(dv)du
# model.list = list(model="SCV", kmax1, kmax2, lambda0, thetavec1, wvec1, thetavec2, wvec2)
hazf.SCV <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_scv_c(tvec, model.list$lambda0, 
                            model.list$thetavec1, model.list$wvec1,
                            model.list$thetavec2, model.list$wvec2)
  } else {
    lambdavec <- (model.list$lambda0 
                  + as.vector(outer(tvec,model.list$thetavec1,function(t,theta) pmax(0,theta-t))%*%model.list$wvec1)
                  + as.vector(outer(tvec,model.list$thetavec2,function(t,theta) pmax(0,t-theta))%*%model.list$wvec2)
    )
  }
  return(lambdavec)
}
chzf.SCV <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_scv_c(tvec, model.list$lambda0, 
                             model.list$thetavec1, model.list$wvec1,
                             model.list$thetavec2, model.list$wvec2)
  } else {
    clambdavec <- ( model.list$lambda0*tvec
                    + 0.5*sum(model.list$wvec1*model.list$thetavec1^2)
                    - 0.5*as.vector(outer(tvec,model.list$thetavec1,function(t,theta) (pmax(0,theta-t))^2)%*%model.list$wvec1)
                    + 0.5*as.vector(outer(tvec,model.list$thetavec2,function(t,theta) (pmax(0,t-theta))^2)%*%model.list$wvec2)
                   )
  }
  return(clambdavec)
}
invsurvf.SCV <- function(uvec, model.list, use.Cpp=FALSE) {
  # ordering of theta values
  othetavec <- c(model.list$thetavec1, model.list$thetavec2)
  othetavec <- sort(othetavec)
  othetavec <- c(0,othetavec,2.*othetavec[length(othetavec)])
  # hazard at these values
  hvec <- hazf(othetavec,model.list,use.Cpp)
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # location of -log(uvec) in the integrated hazards
  kvec <- findInterval(-log(uvec),ivec)
  # quadratic interpolation
  tvec <- solveptquadratic(-log(uvec),
                           othetavec[kvec],othetavec[kvec+1],
                           ivec[kvec],ivec[kvec+1],
                           hvec[kvec]) 
  return(tvec)
}

####################################################
# MBT
# lambda(t) = pi f1(t) + (1-pi) f2(t) / (pi Fbar1(t) + (1-pi) Fbar2(t))
# model.list = list(model="MBT", kmax1, kmax2, pival lambda0, thetavec1, wvec1, thetavec2, wvec2)
hazf.MBT <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_mbt_c(tvec, model.list$pival, model.list$lambda0, 
                            model.list$thetavec1, model.list$wvec1,
                            model.list$thetavec2, model.list$wvec2)
  } else {
    model.list1 <- list(model="DFR", kmax=model.list$kmax1,
                        lambda0=model.list$lambda0, 
                        thetavec=model.list$thetavec1, wvec=model.list$wvec1)
    model.list2 <- list(model="IFR", kmax=model.list$kmax2,
                        lambda0=0.0, 
                        thetavec=model.list$thetavec2, wvec=model.list$wvec2)
    hvec1 <- hazf(tvec, model.list1)
    hvec2 <- hazf(tvec, model.list2)
    ivec1 <- chzf(tvec, model.list1)
    ivec2 <- chzf(tvec, model.list2)
    fbar1 <- exp(-ivec1)
    fbar2 <- exp(-ivec2)
    lambdavec <- ( model.list$pival*hvec1*fbar1 + (1-model.list$pival)*hvec2*fbar2 )/(
                   model.list$pival*fbar1       + (1-model.list$pival)*fbar2 )
  }
  return(lambdavec)
}
chzf.MBT <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_mbt_c(tvec, model.list$pival, model.list$lambda0, 
                             model.list$thetavec1, model.list$wvec1,
                             model.list$thetavec2, model.list$wvec2)
  } else {
    model.list1 <- list(model="DFR", kmax=model.list$kmax1,
                        lambda0=model.list$lambda0, 
                        thetavec=model.list$thetavec1, wvec=model.list$wvec1)
    model.list2 <- list(model="IFR", kmax=model.list$kmax2,
                        lambda0=0.0, 
                        thetavec=model.list$thetavec2, wvec=model.list$wvec2)
    ivec1 <- chzf(tvec, model.list1)
    ivec2 <- chzf(tvec, model.list2)
    fbar1 <- exp(-ivec1)
    fbar2 <- exp(-ivec2)
    clambdavec <- -log( model.list$pival*fbar1 + (1-model.list$pival)*fbar2 )
  }
  return(clambdavec)
}
invsurvf.MBT <- function(uvec, model.list, use.Cpp=FALSE) {
  n <- length(uvec)
  # ordering of theta values
  othetavec <- c(model.list$thetavec1, model.list$thetavec2)
  othetavec <- sort(othetavec)
  ntt <- length(othetavec)
  ttmax <- 1e6*othetavec[ntt]
  othetavec <- c(0,othetavec,ttmax)
  ntt <- ntt+2
  # integrated hazard at these values
  ivec <- chzf(othetavec,model.list,use.Cpp)
  # largest non-infinite value
  imax <- min(c(ntt,which(ivec==Inf)-1))
  # location of -log(uvec) in the integrated hazards
  nloguvec <- -log(uvec)
  kvec <- findInterval(nloguvec,ivec)
  # root finding
  ff <- function(tt,f0) {
    chzf(tt,model.list,use.Cpp)-f0
  }
  #print(cbind(uvec,-log(uvec),
  #            kvec,
  #            ivec[kvec],ivec[kvec+1],
  #            ff(othetavec[kvec],-log(uvec)),
  #            ff(othetavec[kvec+1],-log(uvec)))) ##!!==
  tvec <- sapply(1:n,
                 function(i) {
                   if(kvec[i]>=imax) {
                     ttmax 
                   } else if(ivec[kvec[i]]==nloguvec[i]) {
                     ivec[kvec[i]]
                   } else {
                     uniroot(ff, 
                             interval=othetavec[kvec[i]+c(0,1)],
                             f0=nloguvec[i])$root
                   }
                 })
  return(tvec)
}
rgen.MBT <- function(n, model.list, use.Cpp=FALSE, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  model.list1 <- list(model="DFR", kmax=model.list$kmax1,
                      lambda0=model.list$lambda0, 
                      thetavec=model.list$thetavec1, wvec=model.list$wvec1)
  model.list2 <- list(model="IFR", kmax=model.list$kmax1,
                      lambda0=0.0, 
                      thetavec=model.list$thetavec2, wvec=model.list$wvec2)
  n1 <- rbinom(1,n,model.list$pival)
  n2 <- n-n1
  t1vec <- rgen(n1, model.list1, use.Cpp=use.Cpp)
  t2vec <- rgen(n2, model.list2, use.Cpp=use.Cpp)
  tvec <- sample(c(t1vec,t2vec))
  return(tvec)
}


####################################################
# LCV
# lambda(t) = lambda0 * exp( w0*t + int_0^t G(du) )
# model.list = list(model="LCV", lambda0, w0, thetavec, wvec)
hazf.LCV <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_lcv_c(tvec, model.list$lambda0, model.list$w0, model.list$thetavec, model.list$wvec)
  } else {
    loglambdavec <- log(model.list$lambda0) + model.list$w0*tvec
    loglambdavec <- loglambdavec + as.vector(outer(tvec,model.list$thetavec,function(t,theta) pmax(0,t-theta))%*%model.list$wvec)
    lambdavec <- exp(loglambdavec)
  }
  return(lambdavec)
}
chzf.LCV <- function(tvec, model.list, use.Cpp=FALSE, epsilon=100*.Machine$double.neg.eps) {
  if(use.Cpp) {
    clambdavec <- chzf_lcv_c(tvec, model.list$lambda0, model.list$w0, model.list$thetavec, model.list$wvec, epsilon)
  } else {
    odx <- order(model.list$thetavec) 
    othetavec <- model.list$thetavec[odx]  # theta*_k (1...K)
    owvec <- model.list$wvec[odx]          # w*_k     (1...K)
    kmax <- length(othetavec)   # K
    othetavec <- c(0,othetavec,1.1*othetavec[kmax]) # 0, thetavec, 1.1*last
    owvec <- c(0,owvec,0)                           # 0, wvec,     0
    kmaxp2 <- kmax+2            # K+2
    k1vec <- apply(outer(tvec, othetavec, function(t,theta) theta<=t),1,sum)
    k1vec <- pmin(pmax(1,k1vec),kmaxp2-1)
    
    s1vec <- cumsum(owvec)  # C
    s01vec <- model.list$w0 + s1vec    # C0
    s2vec <- cumsum(owvec*othetavec)    # D
    ccvec <- model.list$lambda0*exp(-s2vec)/s01vec # lambda0.exp(-D)/C0
    ccvec1 <- ccvec[-kmaxp2]
    s01vec1 <- s01vec[-kmaxp2]
    csvec1 <- ccvec1*(exp(s01vec1*othetavec[-1])-exp(s01vec1*othetavec[-kmaxp2]))
    csvec1 <- ifelse(s01vec1==0, 
                     model.list$lambda0*exp(-s2vec[-kmaxp2])*(othetavec[-1]-othetavec[-kmaxp2]), 
                     csvec1)
    ssvec <- c(0,cumsum(csvec1))
    
    ctvec <- ifelse(s01vec1[k1vec]==0,
                    model.list$lambda0*exp(-s2vec[k1vec])*(tvec-othetavec[k1vec]),
                    ccvec[k1vec]*( exp(s01vec[k1vec]*tvec)
                                   -exp(s01vec[k1vec]*othetavec[k1vec])))
    
    clambdavec <- ( ssvec[k1vec] + ctvec )
  }
  return(clambdavec)
}
invsurvf.LCV <- function(uvec, model.list, use.Cpp=FALSE) {
  n <- length(uvec)
  nloguvec <- -log(uvec)
  
  odx <- order(model.list$thetavec) 
  othetavec <- model.list$thetavec[odx]  # theta*_k (1...K)
  owvec <- model.list$wvec[odx]          # w*_k     (1...K)
  kmax <- length(othetavec)   # K
  othetavec <- c(0,othetavec,1.1*othetavec[kmax]) # 0, thetavec, 1.1*last
  owvec <- c(0,owvec,0)                           # 0, wvec,     0
  kmaxp2 <- kmax+2            # K+2
  
  s1vec <- cumsum(owvec)
  s01vec <- model.list$w0 + s1vec
  s2vec <- cumsum(owvec*othetavec)
  ccvec <- model.list$lambda0*exp(-s2vec)/s01vec
  ccvec1 <- ccvec[-kmaxp2]
  s01vec1 <- s01vec[-kmaxp2]
  csvec1 <- ccvec1*(exp(s01vec1*othetavec[-1])-exp(s01vec1*othetavec[-kmaxp2]))
  csvec1 <- ifelse(s01vec1==0, 
                   model.list$lambda0*exp(-s2vec[-kmaxp2])*(othetavec[-1]-othetavec[-kmaxp2]), 
                   csvec1)
  ssvec <- c(0,cumsum(csvec1))
  
  k1vec <- apply(outer(nloguvec, ssvec, function(nlu,ilamval) ilamval<=nlu),1,sum)
  k1vec <- pmax(1,k1vec)
  
  tvec <- ifelse(s01vec[k1vec]==0,
                 othetavec[k1vec] + (nloguvec-ssvec[k1vec])/(model.list$lambda0*exp(-s2vec[k1vec])),
                 (1/s01vec[k1vec])*log(
                   exp(s01vec[k1vec]*othetavec[k1vec])
                   + (nloguvec-ssvec[k1vec])/ccvec[k1vec]) 
                )
  
  return(tvec)
}  


####################################################
# MEW
# lambda(t) = Z'(t)[ nu + alpha*Z(t)^(alpha-1) ]
# Z(t) = exp((mu*t)^beta)-1
# model.list = list(model="MEW", alpha, beta, mu, nu)
hazf.MEW <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_mew_c(tvec, 
                            model.list$alpha,
                            model.list$beta,
                            model.list$mu,
                            model.list$nu)
  } else {
    zvec <- exp((model.list$mu*tvec)^model.list$beta)-1
    zdvec <- (1+zvec)*(model.list$mu*model.list$beta)*((model.list$mu*tvec)^(model.list$beta-1))
    lambdavec <- zdvec*(model.list$nu + model.list$alpha*zvec^(model.list$alpha-1))
  }
  return(lambdavec)
}
chzf.MEW <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_mew_c(tvec,
                             model.list$alpha,
                             model.list$beta,
                             model.list$mu,
                             model.list$nu)
  } else {
    zvec <- exp((model.list$mu*tvec)^model.list$beta)-1
    clambdavec <- model.list$nu*zvec + (zvec)^(model.list$alpha)
  }
  return(clambdavec)
}
invsurvf.MEW <- function(uvec, model.list, use.Cpp=FALSE) {
  fz <- function(expnz,logu,nu,alpha) {
    # expnz is in [0,1]
    z <- -log(expnz)
    nu*z + z^alpha + logu
  }
  expnz <- sapply(uvec, function(u) {
                     uniroot(fz,interval=c(0,1),
                             logu=log(u),
                             nu=model.list$nu,alpha=model.list$alpha)$root
                  })
  tvec <- (1/model.list$mu)*( log(1-log(expnz)) )^(1/model.list$beta)
  return(tvec)
}

####################################################
