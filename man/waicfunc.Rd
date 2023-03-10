% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcfuncs.R
\name{waicfunc}
\alias{waicfunc}
\title{Evaluate the WAIC value from a chain}
\usage{
waicfunc(smat, datlist, fpar, ppar, model)
}
\arguments{
\item{smat}{Matrix output of an MCMC chain}

\item{datlist}{Data object}

\item{fpar}{Fixed parameters of the model}

\item{ppar}{Parameters used in running the sampler}

\item{model}{Model name}
}
\description{
Evaluate the WAIC value from a chain
}
