#'Key to convert from commercial species codes to research vessel species codes
#'
#'Used in \code{LANDdataframe()} to replace commercial species codes with
#'research vessel species codes. Species that have one commercial code but more
#'than one research vessel code are accounted for.
#'
#'@format A dataframe with 108 rows and 3 variables: \describe{
#'  \item{ALLCODES}{Commercial species codes} \item{SPECIES}{Research vessel
#'  codes} \item{PROPORTION_OF_LANDINGS}{The proportion contributed to the
#'  landings of each species (to account for species with one commercial code
#'  but more than one research vessel code). }  }
#'@source These species codes are from AC's file ExtractIndicators/extra
#'  info/indiseas_allcodes2res.csv.
"species_codes"
