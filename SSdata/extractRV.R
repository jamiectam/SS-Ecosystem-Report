#'@title Calls functions to extract and stratify fishery independent data
#'@description Runs \code{biomassData()} to extract survey data and
#'  \code{stratifyBiomass()} to stratify the biomass and abundance estimates.
#'@details User must define \code{channel = odbcConnect("ptran", uid = ###, pwd
#'  = ###)} in the global environment. This channel must have access to the
#'  gscat, gsdet, gsinf, and gs_lengths tables from the groundfish database and
#'  the nafo_strat table from the mfd_stomach database.
#'@inheritParams biomassData
#'@inheritParams stratifyBiomass
#'@param areas Areas (spatial scales) for which to extract and stratify the
#'  fishery independent survey data. A separate dataframe will be exported for
#'  each area. Options are "shelf", "esswss", "nafo", "strat", or any
#'  combination of those four.  Default is \code{areas = c("shelf", "esswss",
#'  "nafo", "strat")}.
#'@return Returns data extracted by \code{biomassData()} and stratified by
#'  \code{stratifyBiomass()}.
#'
#'  From \code{biomassData()}:
#'
#'  Not length-based data are stored in path/data/aggregate/. This folder
#'  includes an RData file for each year called year.RData (object name
#'  \code{dat}). \code{dat} has 9 columns: \code{MISSION}, \code{SETNO},
#'  \code{SPECIES}, \code{YEAR}, \code{STRAT}, \code{BIOMASS}, \code{ABUNDANCE},
#'  \code{QBIOMASS} and \code{QABUNDANCE}.
#'
#'  Length-based data is stored in path/data/length/.  This folder includes an
#'  RData file for each year called num_biom_at_length_year.RData (object name
#'  \code{out}). \code{out} has 10 columns: \code{MISSION}, \code{SETNO},
#'  \code{SPECIES}, \code{YEAR}, \code{STRAT}, \code{BIOMASS}, \code{ABUNDANCE},
#'  \code{QBIOMASS} and \code{QABUNDANCE}, and \code{FLEN}, where \code{FLEN} is
#'  length in 1 cm increments.
#'
#'  From \code{stratifyBiomass()}:
#'
#'  data/stratified/folder/. The name of "folder" depends on the arguments
#'  \code{lengthbased} and \code{qadjusted}:
#'
#'  if(lengthbased & qadjusted), then folder is named "lengthbased_qadj"
#'
#'  if(lengthbased & !qadjusted), then folder is named "lengthbased_notqadj"
#'
#'  if(!lengthbased & qadjusted), then folder is named "notlengthbased_qadj"
#'
#'  if(!lengthbased & qadjusted), then folder is named "notlengthbased_notqadj"
#'
#'  Inside "folder" is a folder for each entry in \code{area}. Within each area
#'  folder is a .RData file for each year from \code{s.year} to \code{e.year}
#'  called year.RData (object name \code{out}).
#'
#'  If \code{lengthbased = TRUE}, \code{out} has 5 columns : \code{SPECIES} ,
#'  \code{ID}, \code{FLEN} (cm), \code{BIOMASS} (kg), and \code{ABUNDANCE}
#'  (numbers).
#'
#'  If \code{lengthbased = FALSE}, \code{out} has 4 columns: \code{SPECIES},
#'  \code{ID}, \code{BIOMASS} (kg), and \code{ABUNDANCE} (numbers).
#'
#'@references Original code by DD.
#'@family RV functions
#'@export

extractRV <- function(path, s.year, e.year, lengthbased, qadjusted,
                           areas = c("shelf", "esswss", "nafo", "strat")){
  
  # Extract biomass & abundance data
  biomassData(path = path, s.year = s.year, e.year = e.year)
  
  # Stratify biomass & abundance data
  stratifyBiomass(path = path, s.year = s.year, e.year = e.year, 
                      lengthbased = lengthbased, qadjusted = qadjusted, areas = areas)
  
  print("Biomass data extracted & stratified")
}

