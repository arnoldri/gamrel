####################################################
# Generic
hazf <- function(t, model.list) {
  switch(model.list$model, 
         CON=hazf.CON(t, model.list)
  )
}
chzf <- function(t, model.list) {
  switch(model.list$model, 
         CON=chzf.CON(t,model.list)
  )
}
survf <- function(t, model.list) {
  exp(-chzf(t, model.list))
}
plot.hazf <- function(model.list, xlim, n=101, add=FALSE, 
                      xlab="t", ylab=expression(lambda(t)), ...) {
  curve(hazf(x,model.list), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}
plot.chzf <- function(model.list, xlim, n=101, add=FALSE, 
                      xlab="t", ylab=expression(Lambda(t)), ...) {
  curve(chzf(x,model.list), xlim=xlim, 
        n=101, add=add, xlab=xlab, ylab=ylab, ...)
  mtext(model.list$model, side=3, line=0, adj=1, cex=0.5)
}

####################################################
# CON
# lambda(t) = lambda0
hazf.CON <- function(t, model.list) rep(model.list$lambda0, length=length(t))
chzf.CON <- function(t, model.list) model.list$lambda0*t
invsurvf.CON <- function(u, model.list) -log(u)/model.list$lambda0

####################################################
