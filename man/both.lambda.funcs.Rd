% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hazfuncs.R
\name{both.lambda.funcs}
\alias{both.lambda.funcs}
\title{Evaluate both the hazard and integrated hazard functions at the same locations}
\usage{
both.lambda.funcs(
  tvec,
  model.list,
  use.Cpp = TRUE,
  epsilon = .Machine$double.neg.eps * 100
)
}
\arguments{
\item{tvec}{Locations at which to evaluate the functions}

\item{model.list}{Model specification}

\item{use.Cpp}{Use cpp to evaluate the functions?}
}
\description{
Hazard rate and Integrated hazard rate function generated by 
a set of discrete locations and weights and integrated
}
