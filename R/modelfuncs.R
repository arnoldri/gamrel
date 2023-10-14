###############################################################################
# Model functions
###############################################################################
#' Initialise objects for an MCMC chain
#'
#' @export
init.objects <- function(tvec, obs,
                         kmax=100,
                         prior.par=NULL,
                         update.par=NULL,
                         model="IFR",
                         datscale=1,
                         epsilon=.Machine$double.neg.eps*100, # tiny value if needed
                         use.Cpp=TRUE,
                         seed=NULL,
                         generate="fixed") {  # generate can be "fixed" or "random"
  if(!is.null(seed)) set.seed(seed)
  if(model%in%c("IFR","DFR")) {

    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        f1=2, f2=2*datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.log.eta=0.3,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.w=0.3,
                         sd.log.alpha=0.3)
    }
    # fixed parameters
    parnames <- c("eta","gamma","thetavec","vvec","alpha","beta","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta (lambda0/gamma)
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames,"wvec","thetaswap")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(eta=1/prior.par$nu,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   phi=prior.par$f1/prior.par$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(eta=rexp(1,prior.par$nu),
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rexp(kmax, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)

  } else if(model%in%c("LWB")) {

    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        c1=1, c2=2/datscale,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        f1=2, f2=2*datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.log.eta=0.3,
                         sd.log.a=0.1,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.w=0.3,
                         sd.log.alpha=0.3)
    }
    # fixed parameters
    parnames <- c("eta","a","gamma","thetavec","vvec","alpha","beta","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta (lambda0/gamma)
                 c1=prior.par$c1, c2=prior.par$c2, # prior for the cutpoint a
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames,"wvec","thetaswap")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.a=update.par$sd.log.a,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(eta=1/prior.par$nu,
                   a=prior.par$c1/prior.par$c2,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=prior.par$a1/prior.par$a2,
                   beta=prior.par$b1/prior.par$b2,
                   phi=prior.par$f1/prior.par$b2)
      epar$gamma <- epar$alpha/epar$beta
    } else {
      epar <- list(eta=rexp(1,prior.par$nu),
                   a=prior.par$c1/prior.par$c2,
                   gamma=NA,
                   thetavec=NA,
                   vvec=NA,
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   beta=rgamma(1,prior.par$b1,prior.par$b2),
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rexp(kmax, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)

  } else if(model%in%c("SBT")) {

    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        f11=1, f21=datscale/2, # DFR
                        f12=2, f22=2*datscale) # IFR
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.log.eta=0.1,
                         sd.log.gamma=0.1,
                         sd.log.theta=0.1,
                         sd.logit.v=0.1,
                         sd.log.w=0.1,
                         sd.log.alpha=0.1)
    }
    # fixed parameters
    parnames <- c("eta",
                  "gamma1","thetavec1","vvec1","alpha1","beta1","phi1",
                  "gamma2","thetavec2","vvec2","alpha2","beta2","phi2")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta (lambda0/gamma2)
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha1,alpha2
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta1,beta2
                 f11=prior.par$f11, f21=prior.par$f21, # prior for phi1
                 f12=prior.par$f12, f22=prior.par$f22, # prior for phi2
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames,"wvec1","thetaswap1","wvec2","thetaswap2")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(eta=1/prior.par$nu,
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=prior.par$a1/prior.par$a2,
                   beta1=prior.par$b1/prior.par$b2,
                   phi1=prior.par$f11/prior.par$f21,
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=prior.par$a1/prior.par$a2,
                   beta2=prior.par$b1/prior.par$b2,
                   phi2=prior.par$f12/prior.par$f22)
      epar$gamma1 <- epar$alpha1/epar$beta1
      epar$gamma2 <- epar$alpha2/epar$beta2
    } else {
      epar <- list(eta=rexp(1,prior.par$nu),
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=rgamma(1,prior.par$a1,prior.par$a2),
                   beta1=rgamma(1,prior.par$b1,prior.par$b2),
                   phi1=rgamma(1,prior.par$f11,prior.par$f21),
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=rgamma(1,prior.par$a1,prior.par$a2),
                   beta2=rgamma(1,prior.par$b1,prior.par$b2),
                   phi2=rgamma(1,prior.par$f12,prior.par$f22))
      epar$gamma1 <- rgamma(1,epar$alpha1,epar$beta1)
      epar$gamma2 <- rgamma(1,epar$alpha2,epar$beta2)
    }
    epar$thetavec1 <- rexp(kmax, epar$phi1)
    epar$vvec1 <- rbeta.t(kmax, 1, epar$alpha1)
    epar$thetavec2 <- rexp(kmax, epar$phi2)
    epar$vvec2 <- rbeta.t(kmax, 1, epar$alpha2)

  } else if(model%in%c("MBT")) {

    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        f11=1, f21=datscale/2,
                        f12=2, f22=2*datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.logit.pival=0.1,
                         sd.log.eta=0.3,
                         sd.log.gamma=0.1,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.w=0.3,
                         sd.log.alpha=0.3)
    }
    # fixed parameters
    parnames <- c("pival",
                  "eta1","gamma1","thetavec1","vvec1","alpha1","beta1","phi1",
                  "eta2","gamma2","thetavec2","vvec2","alpha2","beta2","phi2")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta1,eta2 (lambda0/gamma)
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha1,alpha2,
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta1,beta2
                 f11=prior.par$f11, f21=prior.par$f21, # prior for phi1
                 f12=prior.par$f12, f22=prior.par$f22, # prior for phi1
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames,"wvec1","thetaswap1","wvec2","thetaswap2")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    update["eta1"] <- FALSE # eta1 is not a real parameter - it is always zero
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.logit.pival=update.par$sd.logit.pival,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(pival=0.5,
                   eta1=0,
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=prior.par$a1/prior.par$a2,
                   beta1=prior.par$b1/prior.par$b2,
                   phi1=prior.par$f11/prior.par$f21,
                   eta2=2/prior.par$nu,
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
                   eta1=0,
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=rgamma(1,prior.par$a1,prior.par$a2),
                   beta1=rgamma(1,prior.par$b1,prior.par$b2),
                   phi1=rgamma(1,prior.par$f11,prior.par$f21),
                   eta2=rexp(1,prior.par$nu),
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=rgamma(1,prior.par$a1,prior.par$a2),
                   beta2=rgamma(1,prior.par$b1,prior.par$b2),
                   phi2=rgamma(1,prior.par$f12,prior.par$f22))
      epar$gamma1 <- rgamma(1,epar$alpha1,epar$beta1)
      epar$gamma2 <- rgamma(1,epar$alpha2,epar$beta2)
    }
    epar$thetavec1 <- rexp(kmax, epar$phi1)
    epar$vvec1 <- rbeta.t(kmax, 1, epar$alpha1)
    epar$thetavec2 <- rexp(kmax, epar$phi2)
    epar$vvec2 <- rbeta.t(kmax, 1, epar$alpha2)

  } else if(model%in%c("LCV")) {

    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        s1=1, s2=datscale/100.,
                        sigmap.w0=100./datscale,
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        f1=1, f2=datscale)
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.w0=1./datscale,
                         sd.log.gamma=0.05,
                         sd.log.theta=0.3,
                         sd.logit.v=0.3,
                         sd.log.w=0.05,
                         sd.log.alpha=0.3)
    }
    # fixed parameters
    parnames <- c("lambda0","w0","gamma","thetavec","vvec","alpha","beta","phi")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 s1=prior.par$s1, s2=prior.par$s2, # prior for lambda0
                 sigmap.w0=prior.par$sigmap.w0,    # prior for w0
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames,"wvec","thetaswap")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.w0=update.par$sd.w0,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
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
                   phi=rgamma(1,prior.par$f1,prior.par$b2))
      epar$gamma <- rgamma(1,epar$alpha,epar$beta)
    }
    epar$thetavec <- rexp(kmax, epar$phi)
    epar$vvec <- rbeta.t(kmax, 1, epar$alpha)

  } else if(model%in%c("CVX")) {
    
    if(is.null(prior.par)) {
      prior.par <- list(nu=1,
                        c1=1, c2=2/datscale, 
                        a1=4, a2=1,
                        b1=1, b2=1/datscale,
                        f1=2, f2=2*datscale) # IFR
    }
    if(is.null(update.par)) {
      update.par <- list(psweep=0.1, # probability of the sweep move
                         sd.log.eta=0.1,
                         sd.log.tau=0.1,
                         sd.log.gamma=0.1,
                         sd.log.theta=0.1,
                         sd.logit.v=0.1,
                         sd.log.w=0.1,
                         sd.log.alpha=0.1)
    }
    # fixed parameters
    parnames <- c("eta","tau",
                  "gamma1","thetavec1","vvec1","alpha1","beta1",
                  "gamma2","thetavec2","vvec2","alpha2","beta2","phi2")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 kmax=kmax,                        # sum truncation point
                 nu=prior.par$nu,                  # prior for eta (lambda0/gamma2)
                 c1=prior.par$c1, c2=prior.par$c2, # prior for tau
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha1,alpha2
                 b1=prior.par$b1, b2=prior.par$b2, # prior for beta1,beta2
                 f1=prior.par$f1, f2=prior.par$f2, # prior for phi
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames,"wvec1","thetaswap1","wvec2","thetaswap2")
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 ksweep=FALSE, # = all support points are updated each time
                 ksim=min(kmax,max(5,round(kmax/5))),  # only used if *not* doing a sweep update
                 kswap=min(kmax,max(5,round(kmax/5))), # number of theta values to swap
                 psweep=update.par$psweep,
                 sd.log.eta=update.par$sd.log.eta,
                 sd.log.tau=update.par$sd.log.tau,
                 sd.log.gamma=update.par$sd.log.gamma,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.logit.v=update.par$sd.logit.v,
                 sd.log.w=update.par$sd.log.w,
                 sd.log.alpha=update.par$sd.log.alpha,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(eta=1/prior.par$nu,
                   tau=prior.par$c1/prior.par$c2,
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=prior.par$a1/prior.par$a2,
                   beta1=prior.par$b1/prior.par$b2,
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=prior.par$a1/prior.par$a2,
                   beta2=prior.par$b1/prior.par$b2,
                   phi=prior.par$f1/prior.par$f2)
      epar$gamma1 <- epar$alpha1/epar$beta1
      epar$gamma2 <- epar$alpha2/epar$beta2
    } else {
      epar <- list(eta=rexp(1,prior.par$nu),
                   tau=rgamma(1,prior.par$c1,prior.par$c2),
                   gamma1=NA,
                   thetavec1=NA,
                   vvec1=NA,
                   alpha1=rgamma(1,prior.par$a1,prior.par$a2),
                   beta1=rgamma(1,prior.par$b1,prior.par$b2),
                   gamma2=NA,
                   thetavec2=NA,
                   vvec2=NA,
                   alpha2=rgamma(1,prior.par$a1,prior.par$a2),
                   beta2=rgamma(1,prior.par$b1,prior.par$b2),
                   phi=rgamma(1,prior.par$f1,prior.par$f2))
      epar$gamma1 <- rgamma(1,epar$alpha1,epar$beta1)
      epar$gamma2 <- rgamma(1,epar$alpha2,epar$beta2)
    }
    epar$thetavec1 <- runif(kmax, 0, epar$tau)
    epar$vvec1 <- rbeta.t(kmax, 1, epar$alpha1)
    epar$thetavec2 <- epar$tau + rexp(kmax, epar$phi)
    epar$vvec2 <- rbeta.t(kmax, 1, epar$alpha2)
    
  } else if(model%in%c("MEW")) {
    
    if(is.null(prior.par)) {
      prior.par <- list(s1=1, s2=1,    # lambda
                        a1=2, a2=0.05, # alpha
                        t1=1, t2=2/datscale, # theta
                        g1=1, g2=1)    # gamma
    }
    if(is.null(update.par)) {
      update.par <- list(sd.log.lambda=0.3,
                         sd.log.alpha=0.3,
                         sd.log.theta=0.3,
                         sd.log.gamma=0.3)
    }
    # fixed parameters
    parnames <- c("lambda","alpha","theta","gamma")
    fpar <- list(model=model,                      # model name
                 parnames=parnames,                # parameters
                 s1=prior.par$s1, s2=prior.par$s2, # prior for lambda
                 a1=prior.par$a1, a2=prior.par$a2, # prior for alpha
                 t1=prior.par$t1, t2=prior.par$t2, # prior for theta
                 g1=prior.par$g1, g2=prior.par$g2, # prior for gamma
                 epsilon=epsilon,
                 use.Cpp=use.Cpp)
    # parameters to update
    update_parnames <- c(parnames)
    update <- rep(TRUE, length(update_parnames))
    names(update) <- update_parnames
    # proposal parameters for updates
    ppar <- list(update_parnames=update_parnames,
                 sd.log.lambda=update.par$sd.log.lambda,
                 sd.log.alpha=update.par$sd.log.alpha,
                 sd.log.theta=update.par$sd.log.theta,
                 sd.log.gamma=update.par$sd.log.gamma,
                 update=update,
                 verbose=FALSE,
                 interactive=FALSE)
    # parameters to estimate
    if(generate=="fixed") {
      epar <- list(lambda=prior.par$s1/prior.par$s2,
                   alpha=prior.par$a1/prior.par$a2,
                   theta=prior.par$t1/prior.par$t2,
                   gamma=prior.par$g1/prior.par$g2)
    } else {
      epar <- list(lambda=rgamma(1,prior.par$s1,prior.par$s2),
                   alpha=rgamma(1,prior.par$a1,prior.par$a2),
                   theta=rgamma(1,prior.par$t1,prior.par$t2),
                   gamma=rgamma(1,prior.par$g1,prior.par$g2))
    }

  } else {
    stop("Specified model has not been implemented")
  }
  # data
  datlist <- list(n=length(tvec),      # size of data
                  n0=sum(obs),         # number of uncensored observations
                  tvec=tvec,           # failure/censoring times
                  obs=obs)             # vector of indicators: TRUE=observed, FALSE=censored
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
  if(model%in%c("IFR","DFR")) {
    state <- epar
    # complete state with useful quantities
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec[-fpar$kmax])
    state$uvec <- state$vvec*c(1,cp)
    state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec <- state$gamma * state$uvec
    # compute lambda0
    state$lambda0 <- state$gamma*state$eta

  } else if(model%in%c("LWB")) {
    state <- epar
    # complete state with useful quantities
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec[-fpar$kmax])
    state$uvec <- state$vvec*c(1,cp)
    state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec <- state$gamma * state$uvec
    # compute lambda0
    state$lambda0 <- state$gamma*state$eta

  } else if(model%in%c("SBT")) {
    state <- epar
    # complete state with useful quantities

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

    # compute lambda0 - comes from the IFR component
    state$lambda0 <- state$gamma2*state$eta

  } else if(model%in%c("MBT")) {
    state <- epar
    # complete state with useful quantities

    # ensure eta1=0
    state$eta1 <- 0
      
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec1[-fpar$kmax])
    state$uvec1 <- state$vvec1*c(1,cp)
    state$uvec1[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec1 <- state$gamma1 * state$uvec1
    # compute lambda0
    state$lambda01 <- state$gamma1*state$eta1

    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec2[-fpar$kmax])
    state$uvec2 <- state$vvec2*c(1,cp)
    state$uvec2[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec2 <- state$gamma2 * state$uvec2
    # compute lambda0
    state$lambda02 <- state$gamma2*state$eta2

  } else if(model%in%c("LCV")) {
    state <- epar
    # complete state with useful quantities
    # derive the unscaled weights uvec
    cp <- cumprod(1-state$vvec[-fpar$kmax])
    state$uvec <- state$vvec*c(1,cp)
    state$uvec[fpar$kmax] <- cp[fpar$kmax-1]
    # compute the scaled weights wvec
    state$wvec <- state$gamma * state$uvec
    
  } else if(model%in%c("CVX")) {
    state <- epar
    # complete state with useful quantities
    
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
    
    # compute lambda0 - comes from the IFR component
    state$lambda0 <- state$gamma2*state$eta
    
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

#' Log likelihood of the current state
#'
#' @export
llikef <- function(state, datlist, fpar, model) {
  # log likelihood of observations tvec
  hazfuncs <- both.lambda.funcs(datlist$tvec,
                                model.list=c(list(model=model,kmax=fpar$kmax),
                                             state),
                                use.Cpp=fpar$use.Cpp, epsilon=fpar$epsilon)
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
#' Log likelihood of the current state
#'
#' @export
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
#' Log prior of the current state
#'
#' @export
lpriorf <- function(state, fpar, model) {
  # log prior of the state
  retval <- sum(lpriorf.vector(state, fpar, model))
  return(retval)
}
#' Components of the log prior of the current state - return as a vector
#'
#' @export
lpriorf.vector <- function(state, fpar, model) {
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
  } else if(model%in%c("CVX")) {
    if(fpar$use.Cpp) {
      lprior.vec <- lprior_cvx_c(state$eta, 
                                 state$tau,
                                 state$gamma1,
                                 state$thetavec1, 
                                 state$vvec1,
                                 state$alpha1, 
                                 state$beta1,
                                 state$gamma2,
                                 state$thetavec2, 
                                 state$vvec2,
                                 state$alpha2, 
                                 state$beta2,
                                 state$phi,
                                 fpar$nu, 
                                 fpar$c1, fpar$c2,
                                 fpar$a1, fpar$a2, fpar$b1, fpar$b2, 
                                 fpar$f1, fpar$f2)
      names(lprior.vec) <- fpar$parnames
      #cat("CPP!")
    } else {
      #cat("NOT!")
      # prior for eta
      ## lprior.vec["eta"] <- dexp(state$eta, fpar$nu)
      lprior.vec["eta"] <- -state$eta*fpar$nu 
      # prior for tau
      lprior.vec["tau"] <- (fpar$c1-1)*log(state$tau) - fpar$c2*tau
      
      # prior for gamma1
      lprior.vec["gamma1"] <- dgamma(state$gamma1, state$alpha1, state$beta1, log=TRUE)
      # prior for thetavec1
      lprior.vec["thetavec1"] <- -fpar$kmax*log(state$tau)
      # prior for vvec1
      lprior.vec["vvec1"] <-  sum( (log(state$alpha1) + (state$alpha1-1)*log(1-state$vvec1[-fpar$kmax])) )
      # prior for alpha1
      #lprior.vec["alpha1"] <- dgamma(state$alpha1, fpar$a1, fpar$a2, log=TRUE)
      lprior.vec["alpha1"] <- (fpar$a1-1)*log(state$alpha1) - fpar$a2*state$alpha1
      # prior for beta1
      #lprior.vec["beta1"] <- dgamma(state$beta1, fpar$b1, fpar$b2, log=TRUE)
      lprior.vec["beta1"] <- (fpar$b1-1)*log(state$beta1) - fpar$b2*state$beta1

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
      # prior for phi
      #lprior.vec["phi"] <- dgamma(state$phi2, fpar$f12, fpar$f22, log=TRUE)
      lprior.vec["phi"] <- (fpar$f1-1)*log(state$phi) - fpar$f2*state$phi
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
#' Augment the fixed parameter object (fpar)
#'
#' @export
augment.fpar <- function(state, fpar, model) {
  # add some extra stuff to fpar
  fpar$statenames <- names(state)
  fpar$stackind <- stack(state)$ind
  fpar$statevnames <- names(unlist(sapply(sapply(state,length),
                                          function(i) 1:i)))
  return(fpar)
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
  if(model=="IFR") {
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
      plot(NA,NA, xlim=range(c(state$lambda,state$alpha,state$theta,state$gamma)),
           ylim=c(0,1),
           xlab="parameters", ylab="", main=main, ...)
    }
    points(c(state$lambda,state$alpha,state$theta,state$gamma), c(0,0,0,0), ...)
    
  } else {
    stop(paste0("Model ",model," not recognised"))
  }
  invisible()
}


