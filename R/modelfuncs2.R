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
# Generic
hazf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=hazf.CON(tvec, model.list, use.Cpp),
         IFR=hazf.IFR(tvec, model.list, use.Cpp),
         DFR=hazf.DFR(tvec, model.list, use.Cpp),
         CIR=hazf.CIR(tvec, model.list, use.Cpp),
         CDR=hazf.CDR(tvec, model.list, use.Cpp),
         LWB=hazf.LWB(tvec, model.list, use.Cpp),
         HBT=hazf.HBT(tvec, model.list, use.Cpp),
         MEW=hazf.MEW(tvec, model.list, use.Cpp)
  )
}
chzf <- function(tvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=chzf.CON(tvec, model.list, use.Cpp),
         IFR=chzf.IFR(tvec, model.list, use.Cpp),
         DFR=chzf.DFR(tvec, model.list, use.Cpp),
         CIR=chzf.CIR(tvec, model.list, use.Cpp),
         CDR=chzf.CDR(tvec, model.list, use.Cpp),
         LWB=chzf.LWB(tvec, model.list, use.Cpp),
         HBT=chzf.HBT(tvec, model.list, use.Cpp),
         MEW=chzf.MEW(tvec, model.list, use.Cpp)
  )
}
survf <- function(tvec, model.list, use.Cpp=FALSE) {
  exp(-chzf(tvec, model.list, use.Cpp))
}
invsurvf <- function(uvec, model.list, use.Cpp=FALSE) {
  switch(model.list$model, 
         CON=invsurvf.CON(uvec, model.list, use.Cpp),
         IFR=invsurvf.IFR(uvec, model.list, use.Cpp),
         DFR=invsurvf.DFR(uvec, model.list, use.Cpp),
         CIR=invsurvf.CIR(uvec, model.list, use.Cpp),
         CDR=invsurvf.CDR(uvec, model.list, use.Cpp),
         LWB=invsurvf.LWB(uvec, model.list, use.Cpp),
         HBT=invsurvf.HBT(uvec, model.list, use.Cpp),
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

llfunc <- function(tvec, parvec, model, verbose=FALSE) {
  model.list <- as.model.list(parvec, model)
  retval <- sum(densf(tvec, model.list, use.Cpp=TRUE, log=TRUE))
  if(verbose) print(c(parvec,retval))
  return(retval)
}

as.parvec <- function(model.list) {
  if(model.list$model=="CON") {
    parvec <- log(model.list$lambda0)
  } else if(model.list$model=="MEW") {
    parvec <- log(c(model.list$alpha,model.list$beta,model.list$mu,model.list$nu))
  } else {
    stop(paste("Model",model.list$model,"not recognised"))
  }
}
as.model.list <- function(parvec,model) {
  if(model=="CON") {
    model.list <- list(model=model, lambda0=exp(parvec))
  } else if(model=="MEW") {
    model.list <- list(model=model, 
                       alpha=exp(parvec[1]), beta=exp(parvec[2]), 
                       mu=exp(parvec[3]), nu=exp(parvec[4]))
  } else {
    stop(paste("Model",model,"not recognised"))
  }
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

####################################################
# IFR
# lambda(t) = lambda0 + int_0^t G(du)
# model.list = list(model="IFR", kmax, lambda0, thetavec, wvec)
hazf.IFR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    lambdavec <- hazf_ifr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    lambdavec <- lambda0 + as.vector(outer(tvec,model.list$thetavec,">")%*%model.list$wvec)
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
    lambdavec <- lambda0 + as.vector(outer(tvec,model.list$thetavec,"<")%*%model.list$wvec)
  }
  return(lambdavec)
}
chzf.DFR <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_dfr_c(tvec, model.list$lambda0, model.list$thetavec, model.list$wvec)
  } else {
    clambdavec <- ( model.list$lambda0*tvec
                    + as.vector(outer(tvec,model.list$thetavec,pmin)%*%model.list$wvec)
    )
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
    lambdavec <- lambda0 + as.vector(outer(tvec,model.list$thetavec,function(t,theta) pmax(0,t-theta))%*%model.list$wvec)
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
    lambdavec <- lambda0 + as.vector(outer(tvec,model.list$thetavec,function(t,theta) pmax(0,theta-t))%*%model.list$wvec)
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
    lambdavec <- lambda0 + ifelse(tvec<model.list$a,
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
    lambdavec <- lambda0 + ifelse(tvec<model.list$a,
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) (t<theta & theta<model.list$a))%*%model.list$wvec),
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) (model.list$a<theta & theta<t))%*%model.list$wvec)
    )
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
    lambdavec <- lambda0 + ifelse(tvec<model.list$a,
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) (t<theta & theta<model.list$a))%*%model.list$wvec),
                                  as.vector(outer(tvec,model.list$thetavec,function(t,theta) (model.list$a<theta & theta<t))%*%model.list$wvec)
    )
  }
  return(lambdavec)
}
chzf.HCV <- function(tvec, model.list, use.Cpp=FALSE) {
  if(use.Cpp) {
    clambdavec <- chzf_hcv_c(tvec, model.list$lambda0, model.list$a, model.list$thetavec, model.list$wvec)
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
invsurvf.HCV <- function(uvec, model.list, use.Cpp=FALSE) {
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
# MEW
# lambda(t) = Z'(t)[ nu + alpha*Z(t)^(alpha-1) ]
# Z(t) = exp((mu*t)^beta)-1
# model.list = list(model="IFR", alpha, beta, mu, nu)
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
