#' Log likelihood of the current state
#'
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
lpriorf <- function(state, fpar, model) {
  # log prior of the state
  retval <- sum(lpriorf.vector(state, fpar, model))
  return(retval)
}
