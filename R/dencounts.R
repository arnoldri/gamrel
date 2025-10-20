#' Danish death counts
#' 
#' Counts of deaths in Denmark by age and sex 2006-2024
#'
#' @docType data
#'
#' @usage data(dencounts)
#'
#' @format A data frame with five columns: Year, Age, Male, Female, Total
#' Numbers of people dying at each age in each year.
#'
#' @keywords datasets
#'
#' @references Statistics Denmark table FOD07.
#' Available at www.statbank.dk/FOD207 (data downloaded on 20 October 2025).
#'
#' @examples
#' data(dencounts)
#' \donttest{plot(dencounts$Age[dencounts$Year==2024],
#'                dencounts$Male[dencounts$Year==2024],
#'                main="Deaths by Age in Denmark in 2024")}
"dencounts"
