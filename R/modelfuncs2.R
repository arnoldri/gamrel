####################################################
# Generic
hazf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=hazf.CON(tvec, model.list, use.Cpp)
  )
}
chzf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=chzf.CON(tvec, model.list, use.Cpp)
  )
}
survf <- function(tvec, model.list, use.Cpp=FALSE) {
  exp(-chzf(tvec, model.list, use.Cpp))
}
invsurvf <- function(uvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=invsurvf.CON(uvec, model.list, use.Cpp)
  )
}

rgen <- function(n, model.list, use.Cpp=FALSE) {
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
    lambdavec <- rep(model.list$lambda0, length=length(t))
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
    tvec <- invsurf_con_c(tvec, model.list$lambda0)
  } else {
    tvec <- -log(uvec)/model.list$lambda0
  }
  return(tvec)
}

####################################################
