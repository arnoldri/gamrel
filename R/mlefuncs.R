#' Fit the MLE of the MEW model
#' 
#' @description Parameters are epar=(lambda,alpha,theta,gamma)
#' 
#' @export
mlefit.mew <- function(epar, datlist, fpar, ppar) {
  parvec <- log(epar)
  names(parvec) <- c("lambda","alpha","theta","gamma")
  o1 <- optim(par=parvec, fn=ofunc.mew, gr=ofunc.gr.mew,
              control=list(fnscale=-1), 
              method="L-BFGS-B",
              datlist=datlist, fpar=fpar, ppar=ppar,
              hessian=TRUE,
              verbose=FALSE)
  o1$est <- exp(o1$par)
  o1$vcov.par <- -solve(o1$hessian)
  o1$vcov.est <- diag(o1$est)%*%o1$vcov.par%*%diag(o1$est)
  o1$se.par <- sqrt(diag(o1$vcov.par))
  o1$se.est <- sqrt(diag(o1$vcov.est))
  o1$ci.par <- cbind(estimate=o1$par, lower=o1$par-1.96*o1$se.par, upper=o1$par+1.96*o1$se.par)
  o1$ci.est <- exp(o1$ci.par)
  return(o1)
}

#' Log Likelihood of the MEW function for use by optim()
#' Temporary version
#' 
#' @description Parameterisation is log(lambda),log(alpha),log(theta),log(gamma)
#' 
#' @export
ofunc1.mew <- function(parvec, datlist, fpar, ppar, verbose=FALSE) {
  epar <- as.list(exp(parvec))
  names(epar) <- c("lambda","alpha","theta","gamma")
  
  state <- make.state(epar, datlist, fpar, ppar, model="MEW")
  retval <- llikef(state, datlist, fpar, model="MEW") 
  if(verbose) print(c(as.vector(unlist(epar)),retval))
  return(retval)
}

#' Log Likelihood of the MEW function for use by optim()
#' 
#' @description Parameterisation is log(lambda),log(alpha),log(theta),log(gamma)
#' 
#' @export
ofunc.mew <- function(parvec, datlist, fpar, ppar, verbose=FALSE) {
  epar <- as.list(exp(parvec))
  names(epar) <- c("lambda","alpha","theta","gamma")
  lambda <- epar$lambda
  alpha <- epar$alpha
  theta <- epar$theta
  gamma <- epar$gamma
  
  evec <- exp( (datlist$tvec/theta)^gamma )
  qvec <- lambda + alpha*(evec-1)^(alpha-1)
  ihvec <- lambda*(evec-1) + (evec-1)^alpha
  hvec <- (gamma/theta)*(datlist$tvec/theta)^(gamma-1)*evec*qvec
  retval <- sum(log(hvec)) - sum(ihvec)
  
  if(verbose) print(c(as.vector(unlist(epar)),retval))
  return(retval)
}

#' Gradient of Log Likelihood of the MEW function for use by optim()
#' 
#' @description Parameterisation is log(lambda),log(alpha),log(theta),log(gamma)
#' 
#' @export
ofunc.gr.mew <- function(parvec, datlist, fpar, ppar, verbose=FALSE) {
  epar <- as.list(exp(parvec))
  names(epar) <- c("lambda","alpha","theta","gamma")
  lambda <- epar$lambda
  alpha <- epar$alpha
  theta <- epar$theta
  gamma <- epar$gamma
  
  evec <- exp( (datlist$tvec/theta)^gamma )
  avec <- evec - 1
  bvec <- (datlist$tvec/theta)^gamma * evec
  cvec <- lambda + alpha*avec^(alpha-1)
  ihvec <- lambda*avec + avec^alpha
  hvec <- cvec * gamma/datlist$tvec * bvec
  
  llval <- sum(log(hvec)) - sum(ihvec)
  
  dllval.dlambda <- sum( 1/cvec - avec )
  dllval.dalpha  <- sum( avec^(alpha-1)*(1+alpha*log(avec))/cvec  - log(avec)*avec^alpha )
  dllval.dtheta  <- (gamma/theta)*sum(
    - alpha*(alpha-1)*avec^(alpha-2)*bvec/cvec # note: this term has wrong sign in original paper
    + bvec*cvec   # note: this term is incorrect in the original paper
    - (datlist$tvec/theta)^gamma 
    - 1
  ) 
  dllval.dgamma <- sum(
    1/gamma 
    + log(datlist$tvec/theta)*bvec*( alpha*(alpha-1)*avec^(alpha-2)/cvec - cvec )
    + log(datlist$tvec/theta)*(1 + (datlist$tvec/theta)^gamma)
  )
  retval <- c(lambda*dllval.dlambda, alpha*dllval.dalpha, theta*dllval.dtheta, gamma*dllval.dgamma)
  
  if(verbose) print(c(as.vector(unlist(epar)),retval))
  return(retval)
}

#' Gradient of integrated hazard function the MEW function
#' 
#' @description Parameterisation is lambda,alpha,theta,gamma
#' 
#' @export
inthaz.gr.mew <- function(tvec, lambda, alpha, theta, gamma) {

  evec <- exp( (tvec/theta)^gamma )
  avec <- evec - 1
  bvec <- (tvec/theta)^gamma * evec
  cvec <- lambda + alpha*avec^(alpha-1)
  ihvec <- lambda*avec + avec^alpha
  hvec <- cvec * gamma/tvec * bvec
  
  dihvec.dlambda <- avec
  dihvec.dalpha  <- log(avec) * avec^alpha
  dihvec.dtheta  <- -gamma*bvec*cvec/theta
  dihvec.dgamma <- cvec*log(tvec/theta)*bvec
  
  retval <- cbind(dihvec.dlambda, dihvec.dalpha, dihvec.dtheta, dihvec.dgamma)
  return(retval)
}
