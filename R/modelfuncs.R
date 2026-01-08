###############################################################################
# Model functions
###############################################################################
#' Initialise objects for an MCMC chain
#'
#' @export
init.objects <- function(tvec, obs=TRUE,
                         kmax=100,
                         prior.par=NULL,
                         update.par=NULL, fix.update=NULL,
                         model="IFR",
                         datscale=NULL,
                         epsilon=.Machine$double.neg.eps*100, # tiny value if needed
                         use.Cpp=TRUE,
                         seed=NULL,
                         generate="fixed") {  # generate can be "fixed" or "random"
  if(!is.null(seed)) set.seed(seed)
  
  if(is.null(datscale)) datscale <- mean(tvec,na.rm=TRUE)

  if(model%in%c("CON")) {

    if(is.null(prior.par)) {
      prior.par <- list(s1=1,s2=datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(sd.log.lambda0=0.3)
    }
    # fixed parameters
    parnames <- c("lambda0")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 s1=prior.par$s1,                  # prior for lambda0
                 s2=prior.par$s2,                    
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    # model parameter names 
    model_parnames <- parnames
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames    
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 sd.log.lambda0=update.par$sd.log.lambda0,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(lambda0=prior.par$s1/prior.par$s2)
    } else {
      epar <- list(lambda0=rgamma(1,prior.par$s1,prior.par$s2))
    }

  } else if(model%in%c("IFR","DFR","CIR","CDR")) {
      
    if(is.null(prior.par)) {
      prior.par <- list(s1=1, s2=datscale,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        g1=2, g2=2,
                        f1=2, f2=2*datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(pequal=0.5, # probability we select a support point with equal probability
                         sd.log.lambda0=0.3,
                         sd.log.gamma=0.3,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.alpha=0.3,
                         sd.log.nu=0.3)
    }
    # fixed parameters
    parnames <- c("lambda0","gamma","thetavec","vvec","alpha","beta","nu","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 s1=prior.par$s1, s2=prior.par$s2, # prior for lambda0
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 g1=prior.par$g1, g2=prior.par$g2, # prior for nu
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    # model parameter names
    model_parnames <- parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 pequal=update.par$pequal,
                 sd.log.lambda0=update.par$sd.log.lambda0,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.nu=update.par$sd.log.nu,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(lambda0=prior.par$s1/prior.par$s2,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   nu=prior.par$g1/prior.par$g2,
                   phi=prior.par$f1/prior.par$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(lambda0=rgamma(1,prior.par$s1,prior.par$s2),
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   nu=rgamma(1,prior.par$g1,prior.par$g2),
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rgamma(kmax, epar$nu, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)

  } else if(model%in%c("LWB","HBT","HCV")) {

    if(is.null(prior.par)) {
      prior.par <- list(s1=2, s2=2, 
                        c1=1, c2=2/datscale,
                        a1=3, a2=3,
                        b1=1, b2=1/datscale,
                        g1=2, g2=2, 
                        f1=2, f2=2*datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(pequal=0.5, # probability we select a support point with equal probability
                         sd.log.lambda0=0.3,
                         sd.log.gamma=0.3,
                         sd.log.a=0.1,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.alpha=0.3,
                         sd.log.nu=0.3)
    }
    # fixed parameters
    parnames <- c("lambda0","a","gamma","thetavec","vvec","alpha","beta","nu","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 s1=prior.par$c1, s2=prior.par$c2, # prior for lambda0
                 c1=prior.par$c1, c2=prior.par$c2, # prior for the cutpoint a
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 g1=prior.par$g1, g2=prior.par$g2, # prior for nu
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    # model parameter names
    model_parnames <- parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 pequal=update.par$pequal,
                 sd.log.lambda0=update.par$sd.log.lambda0,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.a=update.par$sd.log.a,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.nu=update.par$sd.log.nu,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(lambda0=prior.par$s1/prior.par$s2,
                   a=prior.par$c1/prior.par$c2,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   nu=prior.par$f1/prior.par$g2,
                   phi=prior.par$f1/prior.par$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(lambda0=rgamma(1,prior.par$s1,prior.par$s2),
                   a=rgamma(1,prior.par$c1,prior.par$c2),
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   nu=rgamma(1,prior.par$g1,prior.par$g2),
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rgamma(kmax, epar$nu, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)

  } else if(model%in%c("SBT","SCV")) {

    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        g11=1, g21=1, 
                        f11=1, f21=datscale/2, # DFR
                        g12=1, g22=1,
                        f12=2, f22=2*datscale) # IFR
    }
    if(is.null(update.par)) {
      update.par <- list(pequal=0.5, # probability we select a support point with equal probability
                         sd.log.lambda0=0.3,
                         sd.log.gamma=0.1,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.alpha=0.3,
                         sd.log.nu=0.3)
    }
    # fixed parameters
    parnames <- c("lambda01",
                  "gamma1","thetavec1","vvec1","alpha1","beta1","nu1","phi1",
                  "lambda02",
                  "gamma2","thetavec2","vvec2","alpha2","beta2","nu2","phi2")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta (lambda0/gamma2)
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha1,alpha2
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta1,beta2
                 g11=prior.par$g11, g21=prior.par$g21, # prior for nu1
                 f11=prior.par$f11, f21=prior.par$f21, # prior for phi1
                 g12=prior.par$g12, g22=prior.par$g22, # prior for nu2
                 f12=prior.par$f12, f22=prior.par$f22, # prior for phi2
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    update["lambda02"] <- FALSE # not a true parameter, lambda02 is always 0
    # model parameter names
    model_parnames <- parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 pequal=update.par$pequal,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.nu=update.par$sd.log.nu,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(lambda01=prior.par$s1/prior.par$s2,
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=prior.par$a1/prior.par$a2,
                   beta1=prior.par$b1/prior.par$b2,
                   nu1=prior.par$g11/prior.par$g21,
                   phi1=prior.par$f11/prior.par$f21,
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=prior.par$a1/prior.par$a2,
                   beta2=prior.par$b1/prior.par$b2,
                   nu2=prior.par$g12/prior.par$g22,
                   phi2=prior.par$f12/prior.par$f22)
      epar$gamma1 <- epar$alpha1/epar$beta1
      epar$gamma2 <- epar$alpha2/epar$beta2
    } else {
      epar <- list(lambda0=rgamma(1,prior.par$s1,prior.par$s2),
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=rgamma(1,prior.par$a1,prior.par$a2),
                   beta1=rgamma(1,prior.par$b1,prior.par$b2),
                   nu1=rgamma(1,prior.par$g11,prior.par$g21),
                   phi1=rgamma(1,prior.par$f11,prior.par$f21),
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=rgamma(1,prior.par$a1,prior.par$a2),
                   beta2=rgamma(1,prior.par$b1,prior.par$b2),
                   nu2=rgamma(1,prior.par$g12,prior.par$g22),
                   phi2=rgamma(1,prior.par$f12,prior.par$f22))
      epar$gamma1 <- rgamma(1,epar$alpha1,epar$beta1)
      epar$gamma2 <- rgamma(1,epar$alpha2,epar$beta2)
    }
    epar$thetavec1 <- rgamma(kmax, epar$nu1, epar$phi1)
    epar$vvec1 <- rbeta.t(kmax, 1, epar$alpha1)
    epar$thetavec2 <- rgamma(kmax, epar$nu2, epar$phi2)
    epar$vvec2 <- rbeta.t(kmax, 1, epar$alpha2)

  } else if(model%in%c("MBT")) {

    if(is.null(prior.par)) {
      prior.par <- list(s1=1, s2=datscale/100, 
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        g11=2, g21=2,
                        g12=2, g22=2,
                        f11=1, f21=datscale/2,
                        f12=2, f22=2*datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(pequal=0.5, # probability we select a support point with equal probability
                         sd.logit.pival=0.1,
                         sd.log.lambda0=0.3,
                         sd.log.gamma=0.1,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.alpha=0.3,
                         sd.log.nu=0.3)
    }
    # fixed parameters
    parnames <- c("pival",
                  "lambda01","gamma1","thetavec1","vvec1","alpha1","beta1","nu1","phi1",
                  "lambda02","gamma2","thetavec2","vvec2","alpha2","beta2","nu2","phi2")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 s1=prior.par$s1, s2=prior.par$s2, # prior for lambda01
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha1,alpha2,
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta1,beta2
                 g11=prior.par$g11, g21=prior.par$g21, # prior for nu1
                 f11=prior.par$f11, f21=prior.par$f21, # prior for phi1
                 g12=prior.par$g12, g22=prior.par$g22, # prior for nu2
                 f12=prior.par$f12, f22=prior.par$f22, # prior for phi2
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    update["lambda20"] <- FALSE # lambda20 is not a real parameter - it is always zero
    # model parameter names
    model_parnames <- parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 pequal=update.par$pequal,
                 sd.logit.pival=update.par$sd.logit.pival,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.nu=update.par$sd.log.nu,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(pival=0.5,
                   lambda01=prior.par$s11/prior.par$s21,
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=prior.par$a1/prior.par$a2,
                   beta1=prior.par$b1/prior.par$b2,
                   phi1=prior.par$f11/prior.par$f21,
                   lambda02=0,
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=prior.par$a1/prior.par$a2,
                   beta2=prior.par$b1/prior.par$b2,
                   phi2=prior.par$f12/prior.par$f22)
      epar$gamma1 <- epar$alpha1/epar$beta1
      epar$gamma2 <- epar$alpha2/epar$beta2
    } else {
      epar <- list(pival=runif(1),
                   lambda01=rgamma(1,prior.par$s11,prior.par$s21),
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=rgamma(1,prior.par$a1,prior.par$a2),
                   beta1=rgamma(1,prior.par$b1,prior.par$b2),
                   nu1=rgamma(1,prior.par$g11,prior.par$g21),
                   phi1=rgamma(1,prior.par$f11,prior.par$f21),
                   lambda02=0,
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=rgamma(1,prior.par$a1,prior.par$a2),
                   beta2=rgamma(1,prior.par$b1,prior.par$b2),
                   nu2=rgamma(1,prior.par$g12,prior.par$g22),
                   phi2=rgamma(1,prior.par$f12,prior.par$f22))
      epar$gamma1 <- rgamma(1,epar$alpha1,epar$beta1)
      epar$gamma2 <- rgamma(1,epar$alpha2,epar$beta2)
    }
    epar$thetavec1 <- rgamma(kmax, epar$nu1, epar$phi1)
    epar$vvec1 <- rbeta.t(kmax, 1, epar$alpha1)
    epar$thetavec2 <- rgamma(kmax, epar$nu2, epar$phi2)
    epar$vvec2 <- rbeta.t(kmax, 1, epar$alpha2)

  } else if(model%in%c("LCV")) {

    if(is.null(prior.par)) {
      prior.par <- list(s1=1, s2=datscale/100.,
                        sigmap.w0=100./datscale,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        g1=1, g2=1, 
                        f1=1, f2=datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(pequal=0.5, # probability we select a support point with equal probability
                         sd.w0=1./datscale,
                         sd.log.gamma=0.05,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.alpha=0.3,
                         sd.log.nu=0.3)
    }
    # fixed parameters
    parnames <- c("lambda0","w0","gamma","thetavec","vvec","alpha","beta","nu","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 s1=prior.par$s1, s2=prior.par$s2, # prior for lambda0
                 sigmap.w0=prior.par$sigmap.w0,    # prior for w0
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 g1=prior.par$f1, g2=prior.par$f2, # prior for nu
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    # model parameter names
    model_parnames <- parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 pequal=update.par$pequal,
                 sd.log.lambda0=update.par$sd.log.lambda0,
                 sd.w0=update.par$sd.w0,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.nu=update.par$sd.log.nu,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(lambda0=prior.par$s1/prior.par$s2,
                   w0=-prior.par$sigmap.w0,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   nu=prior.par$g1/prior.par$g2,
                   phi=prior.par$f1/prior.par$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(lambda0=rgamma(1,prior.par$s1,prior.par$s2),
                   w0=rnorm(1,0,prior.par$sigmap.w0),
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   nu=rgamma(1,prior.par$g1,prior.par$g2),
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rgamma(kmax, epar$nu, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)

  } else if(model%in%c("MEW")) {
    
    if(is.null(prior.par)) {
      prior.par <- list(a1=2, a2=0.5,      # alpha
                        b1=1, b2=2,        # beta
                        s1=1, s2=datscale, # mu
                        t1=1, t2=1) # nu
    }
    if(is.null(update.par)) {
      update.par <- list(sd.log.alpha=1.0,
                         sd.log.beta=0.3,
                         sd.log.mu=0.3,
                         sd.log.nu=0.3)
    }
    # fixed parameters
    parnames <- c("alpha","beta","mu","nu")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 s1=prior.par$s1, s2=prior.par$s2, # prior for mu
                 t1=prior.par$t1, t2=prior.par$t2, # prior for nu
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    if(!is.null(fix.update)) { # parameters which we want to fix update (TRUE or FALSE)
      fix.update <- fix.update[intersect(names(fix.update),names(update))]
      if(length(fix.update)>0) update[names(fix.update)] <- fix.update
    }
    # model parameter names
    model_parnames <- parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 model_parnames=model_parnames,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.beta=update.par$sd.log.beta,
                 sd.log.mu=update.par$sd.log.mu,
                 sd.log.nu=update.par$sd.log.nu,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   mu=prior.par$s1/prior.par$s2,
                   nu=prior.par$t1/prior.par$t2)
    } else {
      epar <- list(alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   mu=rgamma(1,prior.par$s1,prior.par$s2),
                   nu=rgamma(1,prior.par$t1,prior.par$t2))
    }

  } else {
    stop("Specified model has not been implemented")
  }
  # data
  datlist <- make.datlist(tvec,obs)
  # make the state
  state <- make.state(epar, datlist, fpar, ppar, model)
  # augment fpar
  fpar <- augment.fpar(state, fpar, model)
  # return all these objects
  return(list(epar=epar, datlist=datlist,
              fpar=fpar, ppar=ppar, model=model,
              state=state))
}


####################################################################################

#' Construct the state object from the parameter vectors and data
#'
#' @export
make.state <- function(epar, datlist, fpar, ppar, model) {
  if(model%in%c("CON")) {
      state <- epar

  } else if(model%in%c("IFR","DFR","CIR","CDR")) {
      state <- epar
      # complete state with useful quantities
      # derive the unscaled weights uvec
      cp <- cumprod(1-state$vvec[-fpar$kmax])
      state$uvec <- state$vvec*c(1,cp)
      state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
      # compute the scaled weights wvec
      state$wvec <- state$gamma * state$uvec

  } else if(model%in%c("LWB","HBT","HCV")) {
    state <- epar
    # complete state with useful quantities
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec[-fpar$kmax])
    state$uvec <- state$vvec*c(1,cp)
    state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec <- state$gamma * state$uvec

  } else if(model%in%c("SBT","SCV")) {
    state <- epar
    # complete state with useful quantities

    # ensure lambda02=0
    state$lambd02 <- 0
    
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec1[-fpar$kmax])
    state$uvec1 <- state$vvec1*c(1,cp)
    state$uvec1[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec1 <- state$gamma1 * state$uvec1

    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec2[-fpar$kmax])
    state$uvec2 <- state$vvec2*c(1,cp)
    state$uvec2[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec2 <- state$gamma2 * state$uvec2

  } else if(model%in%c("MBT")) {
    state <- epar
    # complete state with useful quantities

    # ensure lambda02=0
    state$lambd02 <- 0
      
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec1[-fpar$kmax])
    state$uvec1 <- state$vvec1*c(1,cp)
    state$uvec1[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec1 <- state$gamma1 * state$uvec1

    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec2[-fpar$kmax])
    state$uvec2 <- state$vvec2*c(1,cp)
    state$uvec2[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec2 <- state$gamma2 * state$uvec2

  } else if(model%in%c("LCV")) {
    state <- epar
    # complete state with useful quantities
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec[-fpar$kmax])
    state$uvec <- state$vvec*c(1,cp)
    state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec <- state$gamma * state$uvec
    
  } else if(model%in%c("MEW")) {
    state <- epar
    
  } else {
    stop("Specified model has not been implemented")
  }

  # log L of the current state
  state$llike <- llikef(state, datlist, fpar, model)
  # log prior of the current state
  state$lprior <- lpriorf(state, fpar, model)
  # was the last update accepted?
  state$accepted <- rep(0,length(ppar$update))
  names(state$accepted) <- names(ppar$update)
  # counter
  state$count <- 0
  return(state)
}


#' Set default priors
#'
#' @export
set.default.priors <- function(model, datscale) {

   if(model%in%c("CON")) {
       
      prior.par <- list(s1=1, s2=2)
      update.par <- NULL
      fix.update <- c()

   } else if(model%in%c("IFR","DFR","CIR","CDR",
                        "LWB","HBT","HCV",
                        "LCV")) {

      prior.par <- list(s1=2, s2=datscale/2.,
                        a1=5, a2=3,
                        b1=2, b2=1/datscale,
                        g1=2, g2=2,
                        f1=2, f2=2*datscale)
      update.par <- list(pequal=0.5, # probability we select a support point with equal probability
                         sd.log.lambda0=0.3,
                         sd.log.gamma=0.3,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.alpha=0.3,
                         sd.log.nu=0.3)
      fix.update <- c()

      if(model%in%c("CIR","CDR")) {
        prior.par$b1 <- 1; prior.par$b2 <- 2/(datscale^2);
        update.par$sd.log.a <- 0.3;
      }

      if(model%in%c("LWB","HBT")) {
        prior.par$c1 <- 1; prior.par$c2 <- 2/datscale;
        update.par$sd.log.a <- 0.3;
      }

      if(model%in%c("HCV")) {
        prior.par$b1 <- 1; prior.par$b2 <- 2/(datscale^2);
        prior.par$c1 <- 1; prior.par$c2 <- 2/datscale;
        update.par$sd.log.a <- 0.3;
      }

      if(model%in%c("SBT","SCV","MBT")) {
        prior.par$f1 <- NULL; prior.par$f2 <- NULL
        prior.par$g1 <- NULL; prior.par$g2 <- NULL
        prior.par$f11 <- 1; prior.par$f21 <- datscale/2;
        prior.par$f12 <- 2; prior.par$f22 <- 2*datscale;
        prior.par$g11 <- 2; prior.par$g21 <- 2;
        prior.par$g12 <- 2; prior.par$g22 <- 2;
      }

      if(model%in%c("LCV")) {
         prior.par$sigmap.w0=2/datscale;
         update.par$sd.w0=1./datscale;
      }

   } else if(model%in%c("MEW")) {
       
      prior.par <- list(a1=4, a2=0.5,      # alpha
                        b1=2, b2=4,        # beta
                        s1=1, s2=datscale, # mu
                        t1=1, t2=1)        # nu
      update.par <- list(sd.log.alpha=1.0,
                         sd.log.beta=0.3,
                         sd.log.mu=0.3,
                         sd.log.nu=0.3)
      fix.update <- c()
   }

   defaults <- list(prior.par=prior.par, update.par=update.par,
                    fix.update=fix.update)

   return(defaults)
}

################################################################################
#' Plot a representation of the current state
#'
#' @description Plot the current state
#'
#' @param state The current state
#' @param datlist The data list
#' @param fpar The fixed parameters
#' @param ppar Proposal parameters
#' @param model Type of model (IFR or DFR)
#' @param type Type of plot (ignored at the moment)
#' @param add Logical: add the state to an existing plot?
#'
#' @export
plot_state <- function(state, datlist, fpar, ppar, model,
                       type=1, main=NULL, add=FALSE, legend.loc="topright", ...) {
  # plot the state
  if(is.null(main)) {
    main <- sprintf("LP=%.3f; LL=%.3f; LPost=%.3f",
                    state$lprior,state$llike,state$llike+state$lprior)
  }
  if(model=="CON") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=c(-1,1), ylim=c(0,3/mean(datlist$tvec[datlist$obs])),
           xlab="", ylab=bquote(lambda[0]), main=main, ...)
    } 
    points(0, state$lambda0, ...)
    
  } else if(model=="IFR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)

  } else if(model=="DFR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)

  } else if(model=="LWB") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    abline(v=state$a, col="red", lwd=2)
    points(state$thetavec, state$wvec, ...)

  } else if(model=="SBT") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(c(state$thetavec1,state$thetavec2)),
           ylim=range(c(state$wvec1,state$wvec2)),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
      legend(legend.loc, pch=c(16,1), legend=c("DFR","IFR"))
    }
    points(state$thetavec1, state$wvec1, pch=16, ...)
    points(state$thetavec2, state$wvec2, pch=1, ...)

  } else if(model=="MBT") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(c(state$thetavec1,state$thetavec2)),
           ylim=range(c(state$wvec1,state$wvec2)),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
      legend(legend.loc, pch=c(16,1), legend=c("DFR","IFR"))
    }
    points(state$thetavec1, state$wvec1, pch=16, ...)
    points(state$thetavec2, state$wvec2, pch=1, ...)

  } else if(model=="LCV") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)
    
  } else if(model=="CIR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)
    
  } else if(model=="CDR") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(state$thetavec),
           ylim=range(state$wvec),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
    }
    points(state$thetavec, state$wvec, ...)
    
  } else if(model=="CVX") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(c(state$thetavec1,state$thetavec2)),
           ylim=range(c(state$wvec1,state$wvec2)),
           xlab=bquote(theta), ylab=bquote(w), main=main, ...)
      legend(legend.loc, pch=c(16,1), legend=c("DFR","IFR"))
    }
    points(state$thetavec1, state$wvec1, pch=16, ...)
    points(state$thetavec2, state$wvec2, pch=1, ...)
    
  } else if(model=="MEW") {
    if(!add) {
      # start a new plot
      plot(NA,NA, xlim=range(c(state$alpha,state$beta,state$mu,state$nu)),
           ylim=c(0,1),
           xlab="parameters", ylab="", main=main, ...)
    }
    points(c(state$alpha,state$beta,state$mu,state$nu), c(0,0,0,0), 
           pch=c("a","b","m","n"),
           ...)
    
  } else {
    stop(paste0("Model ",model," not recognised"))
  }
  invisible()
}


