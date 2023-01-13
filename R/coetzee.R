#' Failure times of rear dump trucks
#'
#' Data from Coetzee (1996)
#'
#' @docType data
#'
#' @usage data(coetzee)
#'
#' @format A vector of failure times (in hours)
#'
#' @keywords datasets
#'
#' @references Coetzee JL. Reliability degradation and the equipment replacement
#' problem. In: Proc. Int. Conf. of Maintenance Societies (ICOMS-96).
#' Melbourne, 1996, Paper 21.
#'
#' @examples
#' data(coetzee)
#' \donttest{hist(coetzee, breaks=10, xlab="t (hours)", ylab="Frequency",
#'                main="Failure times of rear dump trucks")}
"coetzee"