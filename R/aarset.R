#' Failure times of devices from Aarset (1987)
#'
#' Data on failure times from Aarset (1987)
#'
#' @docType data
#'
#' @usage data(aarset)
#'
#' @format A vector of failure times (in hours)
#'
#' @keywords datasets
#'
#' @references Aarset ML. How to Identify a Bathtub Hazard Rate.
#' In: IEEE Transactions on Reliability, R-36(1), 106-108, 1987.
#'
#' @examples
#' data(aarset)
#' \donttest{hist(aarset, 
#'                xlim=c(0,100), breaks=10, xlab="t (hours)", ylab="Frequency",
#'                main="Failure times of devices from Aarset (1987)")}
"aarset"