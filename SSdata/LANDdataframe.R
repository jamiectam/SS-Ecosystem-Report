#'@title Formats commercial landings data for use in the \code{marindicators}
#'  package
#'@description This function imports data from
#'  path/data/landings/landings.RData, attaches labels for the spatial scales of
#'  interest (shelf, ESS/WSS, and/or NAFO divisions) in column \code{ID}, and
#'  replaces commercial species codes with the research vessel species codes.
#'@details If \code{update_LAND = TRUE}, user must define \code{channel =
#'  odbcConnect("ptran", uid = ###, pwd = ###)} in the global environment. This
#'  channel must have access to: the nafo_summary and nafo_area_codes tables
#'  from the COMLAND database, the sub_trips_XXXX and identified_catches tables
#'  from the cl database, the marfis_catch_effort table from mfd_obfmi database,
#'  and the indiseas_marfis2allcodes table from the gomezc database.
#'
#'  Area ID's are assigned from \code{landings_groupings} (type
#'  \code{?landings_groupings} into the console for more info).
#'
#'  Commercial and research vessel species codes are matched from
#'  \code{species_codes} (type \code{?species_codes} into the console for more
#'  info). The function accounts for the proportion of landings of species with
#'  more than one RV code but one commercial code.
#'
#'@param path The filepath to the data folder created by
#'  \code{extractLandings()}.
#'@param areas Areas (spatial scales) of interest. Options are "shelf",
#'  "esswss", "nafo", or any combination of those three. A separate landings
#'  dataframe will be exported for each area. Default is \code{areas =
#'  c("shelf", "esswss", "nafo")}
#'@param csv Logical value indicating whether to export landings dataframe as a
#'  .csv file. Default is \code{csv = TRUE}.
#'@param rdata Logical value indicating whether to export landings dataframe as
#'  a .RData file. Default is \code{rdata = TRUE}.
#'@param update_LAND Logical parameter indicating whether to run
#'  \code{extractLAND()}. This may be time consuming, so if data is already
#'  extracted, use the default \code{update_LAND = FALSE}. If \code{update_LAND
#'  = TRUE}, user must define \code{channel} in the global environment (see
#'  \code{Details}).
#'@param e.year If \code{update_LAND = TRUE}, \code{e.year} is the final year
#'  for which to extract the landings data. Default is \code{e.year = NULL}.
#'@return This function creates a directory path/output/Landings and stores the
#'  landings dataframe for each area in area_land.RData (object name is
#'  \code{land}) and/or area_land.csv. The dataframe has 4 columns: \code{ID}
#'  (the area ID), \code{YEAR}, \code{SPECIES} (the research vessel species
#'  code) and \code{CATCH} (****UNITS of the corresponding species caught in the
#'  corresponding year and area). These files are formatted for the
#'  \code{marindicators} package.
#'@references Original code by DD.
#'@importFrom utils write.csv
#'@importFrom utils read.csv
#'@family LAND functions
#'@export


LANDdataframe <- function(path, areas = c("shelf", "esswss", "nafo"), update_LAND, e.year,
                          csv = TRUE, rdata = TRUE){
  
  # Extract landings data if it hasn't been extracted already
  if(update_LAND) extractLAND(path = path, e.year = e.year)
  
  # Import data
  load(paste(path, "/data/landings/landings.RData", sep=""))    # load landings data 
  names(landings)[1] <- "ALLCODES"                              # change column name to "ALLCODES"
  
  # table to match NAFO_UNIT (in landings.Rdata) to area ID codes
  grp <- landings_groupings 
  
  # table to match ALLCODES with SPECIES and account for the proportion of landings of each SPECIES
  prop.land.table <- species_codes
  
  # loop over each area
  for (j in 1:length(areas)){
    
    data.j <- landings                                  
    areas.j = areas[j]
    
    wl <- grp[, c(1, which(toupper(areas.j) == names(grp)))]               # extracts columns NAFOSUB and areas.J from grp
    names(wl) <- c('NAFO_UNIT', 'ID')                                      # name columns 
    data.j <- merge(data.j, wl, by = 'NAFO_UNIT')                          # merge landings dataframe with the area data
    
    data.j <- merge(data.j, prop.land.table, by = "ALLCODES")              # merge landings dataframe with species codes data
    data.j$CATCH <- as.numeric(data.j$CATCH)
    data.j$CATCH <- data.j$CATCH * data.j$PROPORTION_OF_LANDINGS           # account for the proportion of landings of ALLCODES of each SPECIES
    
    land <- data.j[, c("ID", "YEAR", "SPECIES", "CATCH")]                  # create dataframe with the columns of interest
    land <- aggregate(CATCH ~ ID + YEAR + SPECIES, data = land, FUN = sum) # added this line Feb 4 2020 
    # above line is needed to aggregate CATCH over the different sub-units in grp (e.g., ESS = 4W,  4VS, 4VSB,  4VSU, etc)
    
    dir.create(paste(path, "/output/Landings", sep = ""), recursive = T, showWarnings = F)
    path.output <- paste(path, "/output/Landings/", areas.j, "_land", sep = "")
    
    # save data as an Excel .csv file
    if(csv) write.csv(land, file = paste(path.output, ".csv",sep=""), row.names = FALSE)
    
    if(rdata) save(land, file = paste(path.output, ".RData", sep=""))
    
  } # end of loop over areas
  
  print("landings dataframe exported")
} # end of function