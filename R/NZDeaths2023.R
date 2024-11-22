#' Deaths in New Zealand in 2023
#'
#' Data from Statistics New Zealand
#'
#' @docType data
#'
#' @usage data(NZDeaths2023)
#'
#' @format Counts of deaths in New Zeland in 2023 by age and sex.
#'
#' @keywords datasets
#'
#' @references Statistics New Zealand \url{https://www.stats.govt.nz/topics/population}
#'
#' @examples
#' data(NZDeaths2023)
#' \donttest{barplot(t(as.matrix(NZDeaths2023[,c("Male","Female")])),
#'    beside=TRUE,names=dcounts$age_group, las=2, cex.lab=1.0, 
#'    cex.names=0.8, cex.axis=0.8,
#'    legend=TRUE, args.legend=list(x="topleft"),
#'    xlab="Age group", ylab="Counts", main="Deaths in 2023")
#'           mtext("Source: Statistics New Zealand", 
#'                 side=1, line=3, adj=1, cex=0.5)}
"NZDeaths2023"
