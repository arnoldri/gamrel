#' Fit the MLE of the MEW model
#' 
#' @description Parameters are epar=(alpha,beta,mu,nu)
#' 
#' @export
mlefit.mew <- function(epar, datlist, fpar, ppar) {
  parvec <- log(epar)
  names(parvec) <- c("alpha","beta","mu","nu")
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
#' @description Parameterisation is log(alpha),log(beta),log(mu),log(nu)
#' 
#' @export
ofunc1.mew <- function(parvec, datlist, fpar, ppar, verbose=FALSE) {
  epar <- as.list(exp(parvec))
  names(epar) <- c("alpha","beta","mu","nu")
  
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
  names(epar) <- c("alpha","beta","mu","nu")
  alpha <- epar$alpha
  beta <- epar$beta
  mu <- epar$mu
  nu <- epar$nu

  evec <- exp( (mu*datlist$tvec)^beta )
  zvec <- evec - 1
  zdvec <- (mu*beta)*(mu*datlist$tvec)^(beta-1)*evec 
  qvec <- nu + alpha*zvec^(alpha-1)
  ihvec <- nu*zvec + zvec^alpha
  hvec <- qvec*zdvec
  
  llval <- sum(log(hvec[datlist$obs])) - sum(ihvec)
  
  if(verbose) print(c(as.vector(unlist(epar)),llval))
  return(llval)
}

#' Gradient of Log Likelihood of the MEW function for use by optim()
#' 
#' @description Parameterisation is log(alpha),log(beta),log(mu),log(nu)
#' 
#' @export
ofunc.gr.mew <- function(parvec, datlist, fpar, ppar, verbose=FALSE) {
  epar <- as.list(exp(parvec))
  names(epar) <- c("alpha","beta","mu","nu")
  alpha <- epar$alpha
  beta <- epar$beta
  mu <- epar$mu
  nu <- epar$nu

  evec <- exp( (mu*datlist$tvec)^beta )
  zvec <- evec - 1
  zdvec <- (mu*beta)*(mu*datlist$tvec)^(beta-1)*evec 
  qvec <- nu + alpha*zvec^(alpha-1)
  ihvec <- nu*zvec + zvec^alpha
  hvec <- qvec*zdvec
  llval <- sum(log(hvec[datlist$obs])) - sum(ihvec)
  
  dzvec.dbeta <- evec*log(mu*tvec)*(mu*datlist$tvec)^beta
  dzvec.dmu <- (beta/mu)*(mu*datlist$tvec)^beta * evec
  
  dzdvec.dbeta <- ( zdvec*(1/beta + log(mu*datlist$tvec))
                    + (mu*beta)*(mu*datlist$tvec)^(beta-1)*dzvec.dbeta )
  dzdvec.dmu <- (beta/mu)*zdvec*(1 + (mu*datlist$tvec)^beta)
  
  dllval.dalpha  <- ( sum( (zvec^(alpha-1)*(1+alpha*log(zvec))/qvec)[datlist$obs] )  
                     -sum( log(zvec)*zvec^alpha )
                    )
  dllval.dbeta <- ( 
    sum( ((alpha*(alpha-1)*zvec^(alpha-2))/(qvec)*dzvec.dbeta)[datlist$obs]) 
    + sum( (dzdvec.dbeta/zdvec)[datlist$obs] )
    - sum( qvec*dzvec.dbeta )
  )
  dllval.dmu <- ( 
    sum( ((alpha*(alpha-1)*zvec^(alpha-2))/(qvec)*dzvec.dmu)[datlist$obs]) 
    + sum( (dzdvec.dmu/zdvec)[datlist$obs] )
    - sum( qvec*dzvec.dmu )
  )
  dllval.dnu <- sum( 1/(qvec[datlist$obs]) ) - sum(zvec)
  
  retval <- c(alpha*dllval.dalpha, beta*dllval.beta, mu*dllval.mu, nu*dllval.nu)
  
  if(verbose) print(c(as.vector(unlist(epar)),retval))
  return(retval)
}

