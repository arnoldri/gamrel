#' Failure times of three Load Haul Dump machines in a Swedish mine
#' 
#' Data on failure times on load haul dump machines from Kumar, Kelfsjo and Granholm (1989)
#'
#' @docType data
#'
#' @usage data(kumarLHD)
#'
#' @format A data frame with machine IDs (A,B,C), observation number,
#' times between failures, faiure times, and a string indicating which
#' systems failed or were maintained.  The system codes are 
#' E Engine;
#' H Hydraulics;
#' O Others including body, cabin, chassis, etc.;
#' R Transmission;
#' T Tyres and Wheels;
#' B Brakes;
#' A Air Conditioning;
#' * Maintenance
#'
#' @keywords datasets
#'
#' @references Kumar U, Klefsj\"{o} B, Granholm S. 
#' Reliability investigation for a fleet
#' of load haul dump machines in a Swedish mine.
#' Reliability Engineering and System Safety 1989;26(4):341-361.
#'
#' @examples
#' data(kumarLHD)
#' \donttest{hist(kumarLHD$t, breaks=10, xlab="t (hours)", ylab="Frequency",
#'                main="Failure times of Load Haul Dump machines")}
"kumarLHD"