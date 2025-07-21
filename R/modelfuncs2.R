####################################################
# Generic
hazf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=hazf.CON(tvec, model.list, use.Cpp),
         MEW=hazf.MEW(tvec, model.list, use.Cpp)
  )
}
chzf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=chzf.CON(tvec, model.list, use.Cpp),
         MEW=chzf.MEW(tvec, model.list, use.Cpp)
  )
}
survf <- function(tvec, model.list, use.Cpp=FALSE) {
  exp(-chzf(tvec, model.list, use.Cpp))
}
invsurvf <- function(uvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=invsurvf.CON(uvec, model.list, use.Cpp),
         MEW=invsurvf.MEW(uvec, model.list, use.Cpp)
  )
}
densf <- function(tvec, model.list, use.Cpp=FALSE, log=FALSE) {
  if(log) {
     retval <- log(hazf(tvec, model.list, use.Cpp)) - chzf(tvec, model.list, use.Cpp)
  } else {
     retval <- hazf(tvec, model.list, use.Cpp)*exp(-chzf(tvec, model.list, use.Cpp))
  }
  return(retval)
}

rgen <- function(n, model.list, use.Cpp=FALSE, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  uvec <- runif(n)
  tvec <- invsurvf(uvec, model.list=model.list, use.Cpp=use.Cpp)
  return(tvec)
}

plot.hazf <- function(model.list, xlim, use.Cpp=FALSE, n=101, add=FALSE, 
                      xlab="t", ylab=expression(lambda(t)), ...) {
  curve(hazf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

plot.chzf <- function(model.list, xlim, use.Cpp=FALSE, n=101, add=FALSE, 
                      xlab="t", ylab=expression(Lambda(t)), ...) {
  curve(chzf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

plot.survf <- function(model.list, xlim, ylim=c(0,1), use.Cpp=FALSE, n=101, add=FALSE, 
                       xlab="t", ylab=expression(bar(F)(t)), ...) {
  curve(survf(x,model.list,use.Cpp), xlim=xlim, ylim=ylim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

plot.densf <- function(model.list, xlim, use.Cpp=FALSE, n=101, add=FALSE, 
                       xlab="t", ylab=expression(f(t)), ...) {
  curve(densf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

plot.invsurvf <- function(model.list, xlim=c(0,1)+0.0001*c(1,-1), use.Cpp=FALSE, n=101, add=FALSE, 
                          xlab=expression(bar(F)(t)), ylab="t", ...) {
  curve(invsurvf(x,model.list,use.Cpp), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

####################################################
# CON
# lambda(t) = lambda0
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

####################################################
# MEW
# lambda(t) = Z'(t)[ nu + alpha*Z(t)^(alpha-1) ]
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
  #if(use.Cpp) {
  #  tvec <- invsurvf_mew_c(uvec, 
  #                         model.list$alpha,
  #                         model.list$beta,
  #                         model.list$mu,
  #                         model.list$nu)
  #} else {
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
  #}
  return(tvec)
}

####################################################
