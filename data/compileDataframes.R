#'@title Exports fishery independent and commercial landings data
#'@description Exports fishery independent and commercial landings dataframes in
#'  a format suitable for the \code{marindicators} package. These dataframes
#'  correspond to \code{X}, \code{X_length}, \code{LenWt.table}, and \code{land}
#'  in the \code{marindicators} package.
#'@details This function calls \code{biomassData()}, \code{stratifyBiomass()},
#'  \code{RVdataframe()}, \code{LWdataframe()}, and \code{LANDdataframe()} to
#'  extract and format the survey, length-weight, and commercial landings data.
#'
#'  User must define \code{channel = odbcConnect("ptran", uid = ###, pwd = ###)}
#'  in the global environment. This channel must have access to the gscat,
#'  gsdet, gsinf, and gs_lengths tables from the groundfish database and the
#'  nafo_strat table from the mfd_stomach database for the fishery-independent
#'  data and the NAFO, ZIF, and MARFIS databases for the commercial data.
#'@inheritParams biomassData
#'@param path Filepath indicating where to create folders to store the output.
#'@param s.year Year for which to begin data compilation.
#'@param e.year Year for which to end data compilation.
#'@param areas.RV Areas (spatial scales) for which to compile fishery
#'  independent data. Options are "strat", "nafo", "esswss", "shelf", or any
#'  combination of those four. Default is \code{areas.RV = c("strat", "nafo",
#'  "esswss", "shelf")}.
#'@param areas.land Areas (spatial scales) for which to compile commercial
#'  landings data. Options are "nafo", "esswss", "shelf", or any combination of
#'  those three.  Default is \code{areas.land = c("nafo", "esswss", "shelf")}.
#'@param csv Logical value indicating whether to export dataframe as a .csv
#'  file. Default is \code{csv = TRUE}.
#'@param rdata Logical value indicating whether to export dataframe as a .RData
#'  file. Default is \code{rdata = TRUE}.
#'@return The output is formatted for the \code{marindicators} package and saved
#'  in directories created by the functions:
#'
#'  \code{RVdataframe} creates a directory path/output/RV. In the RV folder is a
#'  folder for each entry in \code{areas.RV}. In each area folder is a csv
#'  and/or RData file for the specified combination of \code{lengthbased} and
#'  \code{qadjusted}.
#'
#'  \code{LWdataframe()} creates a directory path/output/LengthWeight. In the
#'  LengthWeight folder is a csv and/or RData file with length-weight data for
#'  each area in \code{areas.RV}.
#'
#'  \code{LANDdatafrmae()} creates a directory path/output/Landings: In the
#'  Landings folder is a csv and/or RData file for each area in
#'  /code{areas.land}.
#'@references Original code by DD.
#'@export


compileDataframes <- function(path, s.year, e.year, areas.RV = c("strat", "nafo", "esswss", "shelf"),
                              areas.land = c("nafo", "esswss", "shelf"), csv = TRUE, rdata = TRUE){

  # Extract fishery independent data
  biomassData(path = path, s.year = s.year, e.year = e.year, s.strat = 440, e.strat = 495,
              vessel.correction = TRUE)
  
  # Stratify fishery independent data
  stratifyBiomass(path = path, s.year = s.year, e.year = e.year, lengthbased = TRUE, qadjusted = TRUE,
                  areas = areas.RV)
  stratifyBiomass(path = path, s.year = s.year, e.year = e.year, lengthbased = TRUE, qadjusted = FALSE,
                  areas = areas.RV)
  stratifyBiomass(path = path, s.year = s.year, e.year = e.year, lengthbased = FALSE, qadjusted = TRUE,
                  areas = areas.RV)
  stratifyBiomass(path = path, s.year = s.year, e.year = e.year, lengthbased = FALSE, qadjusted = FALSE,
                  areas = areas.RV)
  
  
  # Format fishery independent data
  RVdataframe(path = path,  s.year = s.year, e.year = e.year,
              areas = areas.RV, lengthbased = TRUE, qadjusted = TRUE, csv = csv, rdata = rdata)
  RVdataframe(path = path,  s.year = s.year, e.year = e.year, 
              areas = areas.RV, lengthbased = TRUE, qadjusted = FALSE, csv = csv, rdata = rdata)
  RVdataframe(path = path,  s.year = s.year, e.year = e.year,
              areas = areas.RV, lengthbased = FALSE, qadjusted = TRUE, csv = csv, rdata = rdata)
  RVdataframe(path = path,  s.year = s.year, e.year = e.year,
              areas = areas.RV, lengthbased = FALSE, qadjusted = FALSE, csv = csv, rdata = rdata)

  # Extract and format length-weight data
  LWdataframe(path = path,  s.year = s.year, e.year = e.year, 
              areas = areas.RV, update_LW = TRUE, csv = csv, rdata = rdata)
  
  # Extract and format landings data
  LANDdataframe(path = path, areas = areas.land, update_LAND = TRUE, e.year = e.year, csv = csv, rdata = rdata)
  
  
}