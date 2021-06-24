#' Strata areas for the Scotian Shelf
#'
#'Strata areas for strata 440 - 478, 480 - 485 and 490 to 495. Each
#'  \code{AREA} is converted to tow units using the conversion \code{TOW UNITS =
#'  AREA/((35./6080.2)*1.75)} until 1981 and \code{TOW UNITS =
#'  AREA/((41./6080.2)*1.75)} after 1981. (Fanning 1985)
#'
#'@format A dataframe with 51 rows and 2 variables: \describe{
#'  \item{STRAT}{Strata number} \item{AREA}{Strata area}  }
#'@source These strata areas are from AC's file ExtractIndicators/extra
#'  info/stratweights.csv.
"strat_weights"
