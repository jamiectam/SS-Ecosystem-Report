#'@title Formats length-weight data for use in the \code{marindicators} package
#'@description This function imports data from path/data/lenwgt adds column
#'  \code{YEAR}, and attaches labels for the spatial scales of interest (shelf,
#'  ESS/WSS, NAFO divisions, and/or strata) in column \code{ID}.
#'@details If \code{update_LW = TRUE}, user must define \code{channel =
#'  odbcConnect("ptran", uid = ###, pwd = ###)} in the global environment. This
#'  channel must have access to the gsinf and gsdet tables from the groundfish database.
#'@inheritParams RVdataframe
#'@param path The filepath to the data folder created by \code{extractLW()}.
#'@param update_LW Logical parameter indicating whether to run
#'  \code{extractLW()}. This can be time consuming, so if LW data is already
#'  extracted, use the default \code{update_LW = FALSE}. If \code{update_LW =
#'  TRUE}, user must define \code{channel} in the global environment (see
#'  \code{Details}).
#'@return This function creates a directory path/output/LengthWeight and stores
#'  the length-weight dataframe for each area in area_LengthWeight.RData (object
#'  name is \code{lw}) and/or area_LengthWeight.csv. The dataframe has 5
#'  columns: \code{ID}, \code{YEAR}, \code{SPECIES}, \code{LENGTH} (cm), and
#'  \code{WEIGHT} (g). These files are formatted for the \code{marindicators}
#'  package.
#'@references Original code by DD.
#'@importFrom utils write.csv
#'@family LW functions
#'@export

LWdataframe <- function(path, s.year, e.year, areas = c("shelf", "esswss", "nafo", "strat"),
                        update_LW = FALSE, csv = TRUE, rdata = TRUE){
  
  # Extract length-weight data if it hasn't been extracted already
  if(update_LW) extractLW(path = path, s.year = s.year, e.year = e.year)
  
  years <- c(s.year:e.year)
  allData <- NULL
  
  for (i in 1:length(years)){
    
    load(paste(path, "/data/lenwgt/lw", years[i],".RData",sep=""))		
    
    data.i = wt
    rm(wt)
    
    year.i = years[i]
    YEAR.i = rep(year.i, nrow(data.i))
    
    data.i = cbind(YEAR.i, data.i)
    names(data.i) = c("YEAR", "STRAT", "SPECIES", "LENGTH", "WEIGHT")
    
    allData <- rbind(allData, data.i)
    
  }
  
  # now use the defineAreas function to identify which spatial scale each strata belongs to
  # i.e., add "ID" column to data
  for(j in 1:length(areas)){
    
    areas.j = areas[j]
    allData_ID = defineAreas(allData, area = areas[j])
    allData_ID$STRAT <- NULL
    lw <- allData_ID
    
    dir.create(paste(path, "/output/LengthWeight", sep = ""), recursive = T, showWarnings = F)
    path.output <- paste(path, "/output/LengthWeight/", areas.j, "_LengthWeight", sep = "")
    
    if(csv) write.csv(lw, file = paste(path.output, ".csv", sep=""), row.names = FALSE)
  
    if(rdata) save(lw, file = paste(path.output, ".RData", sep=""))
    
  }
  
  print("length-weight dataframes exported")
  
}
