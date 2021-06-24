#'Species included in 9 commercial groups
#'
#'A dataframe indicating which species belong in each commercial landings group.
#'\code{commercial_groups} has 11 columns: \code{ALLCODES} (commercial species
#'code), \code{ALLNAMES} (commercial species name), and a column for each
#'commercial group: \code{GROUNDFISH}, \code{CLUPEIDS}, \code{INVERTEBRATES},
#'\code{FORAGE}, \code{FLATFISH}, \code{FINFISH}, \code{LARGE_PELAGICS},
#'\code{GADOIDS}, \code{SKATES}. Species that belong to a commercial group have
#'a \code{1} in the corresponding column and \code{NA} in all other columns.
#'Applied in \code{extractLAND()}.
#'
#'@format A dataframe with 175 rows and 11 variables: \describe{
#'  \item{ALLCODES}{Commercial species code} \item{ALLNAMES}{Common species
#'  name} \item{GROUNDFISH}{For example, a 1 in this column indicates the
#'  corresponding species is included in commercial groundfish landings; an NA
#'  indicates the species is not included in this group.}}
#'@source These commercial groupings are from AC's script LandByGroup.R, u <-
#'  sqlQuery(channel,paste("select * from gomezc.indiseas_allcodes")).
"commercial_groups"
