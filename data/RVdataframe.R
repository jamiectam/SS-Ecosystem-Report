#'@title Formats fishery independent survey data for use in the
#'  \code{marindicators} package
#'@description This function imports data exported by \code{stratifyBiomass()},
#'  adds column \code{YEAR}, and binds all years together.
#'@details If \code{update_RV = TRUE}, user must define \code{channel =
#'  odbcConnect("ptran", uid = ###, pwd = ###)} in the global environment. This
#'  channel must have access to the gscat, gsdet, gsinf, and gs_lengths tables
#'  from the groundfish database and the nafo_strat table from the mfd_stomach
#'  database.
#'@inheritParams biomassData
#'@inheritParams stratifyBiomass
#'@param path The filepath to the /data folder.
#'@param areas Areas (spatial scales) for which to format data. A separate
#'  dataframe will be exported for each area. Options are "shelf", "esswss",
#'  "nafo", "strat", or any combination of the four. Default is \code{areas =
#'  c("shelf", "esswss", "nafo", "strat")}.
#'@param csv Logical value indicating whether to export dataframe as a .csv
#'  file. Default is \code{csv = TRUE}.
#'@param rdata Logical value indicating whether to export dataframe as a .RData
#'  file. Default is \code{rdata = TRUE}.
#'@param update_RV Logical parameter indicating whether to run
#'  \code{extractRV()}. Extracting the abundance and biomass data for all years
#'  and areas is time consuming, so if data is already extracted, use the
#'  default \code{update_biomass = FALSE}. If \code{update_RV = TRUE}, user
#'  must define \code{channel} in the global environment (see \code{Details}).
#'@return This function creates a directory path/output/RV/area for each entry
#'  in \code{area}. In each area folder is a csv and/or RData file for the
#'  specified combination of \code{lengthbased} and \code{qadjusted} (object
#'  name \code{RVdata}). These files are formatted for the \code{marindicators}
#'  package.
#'
#'  If \code{lengthbased = TRUE}, then \code{RVdata} has 6 columns: \code{YEAR},
#'  \code{SPECIES}, \code{ID}, \code{LENGTH} (cm), \code{BIOMASS} (kg), and
#'  \code{ABUNDANCE} (numbers).
#'
#'  If \code{lengthbased = FALSE}, then \code{RVdata} has 5 columns:
#'  \code{YEAR}, \code{SPECIES}, \code{ID}, \code{BIOMASS} (kg), and
#'  \code{ABUNDANCE} (numbers).
#'
#'@references Original code by DD.
#'@importFrom utils write.csv
#'@family RV functions
#'@export

RVdataframe <- function(path, s.year, e.year, areas = c("shelf", "esswss", "nafo", "strat"), 
                        lengthbased, qadjusted, update_RV = FALSE, csv = TRUE, rdata = TRUE){
  
  # Extract and stratify biomass data if it hasn't been done already
  if(update_RV) {
    extractRV(path = path, s.year = s.year, e.year = e.year, 
              areas = areas, 
              lengthbased = lengthbased, qadjusted = qadjusted)
  }
  
  years <- c(s.year:e.year)
  
  for (j in 1:length(areas)){
    
    areas.j = areas[j]
    allData <- NULL
    
    for (i in 1:length(years)){
      
      if(lengthbased & qadjusted) u <- "lengthbased_qadj"
      if(lengthbased & !qadjusted) u <- "lengthbased_notqadj"
      if(!lengthbased & qadjusted) u <- "notlengthbased_qadj"
      if(!lengthbased & !qadjusted) u <- "notlengthbased_notqadj"
      
      path.input <- paste(path, "/data/stratified/", u, "/", areas.j,"/", sep = "")
      
      load(paste(path.input, years[i], ".RData", sep=""))		
      data.i = out
      rm(out)
      
      year.i = years[i]
      YEAR.i = rep(year.i, nrow(data.i))
      
      data.i = cbind(YEAR.i, data.i)
      
      if(lengthbased) names(data.i) = c("YEAR", "SPECIES", "ID", "LENGTH", "BIOMASS", "ABUNDANCE")
      if(!lengthbased) names(data.i) = c("YEAR", "SPECIES", "ID", "BIOMASS", "ABUNDANCE")
      
      allData <- rbind(allData, data.i)
      
    }
    
    dir.create(paste(path, "/output/RV/", areas.j, sep = ""), recursive = T, showWarnings = F)
    path.output <- paste(path, "/output/RV/", areas.j, "/", areas.j, "_", u, sep = "")
    
    RVdata <- allData
    if(csv) write.csv(RVdata, file = paste(path.output, ".csv",sep=""), row.names = FALSE)
    if(rdata) save(RVdata, file = paste(path.output, ".RData", sep=""))
    
  }
  
  print("survey dataframes exported")
  
}