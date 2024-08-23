#'@title Assigns area ID based on strata.
#'@description Adds column \code{ID} to a dataframe that already has a column
#'  \code{STRAT}, and then assigns an area ID based on the strata number in
#'  \code{STRAT}.
#'@param dat A dataframe that includes the column \code{STRAT}, which has
#'  entries corresponding to the strata on the Scotian Shelf (e.g., values from
#'  440 to 495).
#'@param area Character string indicating the spatial scale for which to assign
#'  areas IDs. Options are: \code{"strat"}, \code{"nafo"}, \code{"esswss"}, or
#'  \code{"shelf"}.
#'@return The function returns \code{dat} with an extra column: \code{ID}
#'@references Modified code from AC's ExtractIndicators/R/defineGroups.R.
#'@export

defineAreas <- function(dat, area) {
  
  if(area == "nafo") {
    ar <- data.frame(STRAT = c(440:466,470:478,480:485,490:495),
                     ID = c(rep('4VN',3), rep('4VS',10), rep('4W',14), rep('4X',21)))
  }
  
  if(area == "shelf") {
    ar <- data.frame(STRAT = c(440:466,470:478,480:485,490:495),
                     ID = rep('SHELF',48))
  }
  
  if(area == "esswss") {
    ar <- data.frame(STRAT=c(440:466,470:478,480:485,490:495),
                     ID = c(rep('ESS', 27),rep('WSS', 21)))
  }
  
  if(area == "strat") {
    ar <- data.frame(STRAT = c(440:466,470:478,480:485,490:495),
                     ID = c(440:466,470:478,480:485,490:495))
  }
  
  
  if(any(!is.na(dat$STRAT))) {
    dat <- merge(dat, ar, by = 'STRAT', all.x = T)
  } else {
    print('Data input problems. Make sure dat has a column named STRAT with values from 440 to 495, and area is "nafo", "shelf", "esswss" or "shelf".')
    dat <- dat
  }
  dat
}

