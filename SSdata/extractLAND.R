#'@title Extracts commercial landings data from historic and current databases
#'@description Extracts commercial landings data from NAFO (1968 - 1985), ZIF
#'  (1986 - 2002) and MARFIS (2003 - present) databases.
#'@details User must define \code{channel = odbcConnect("ptran", uid = ###, pwd
#'  = ###)} in the global environment. This channel must have access to: the
#'  nafo_summary and nafo_area_codes tables from the COMLAND database, the
#'  sub_trips_XXXX and identified_catches tables from the cl database, the
#'  marfis_catch_effort table from mfd_obfmi database, and the
#'  indiseas_marfis2allcodes table from the gomezc database.
#'
#'  Units: tonnes 
#'
#'  This function defines and then calls three functions: \code{NAFO()},
#'  \code{ZIF()}, and \code{MARFIS()}.
#'
#'  Extracts data from 1968 to \code{e.year}.
#'
#'  Species are assigned to commercial landings groups based on
#'  \code{commercial_groups}. Type \code{?commercial_groups} for more
#'  information.
#'@inheritParams biomassData
#'@param path  Filepath indicating where to create the folder to store the
#'  extracted data.
#'@return This function creates directory path/data/landings and stores file
#'  landings.RData (object name is \code{landings}). \code{landings} has 14
#'  columns: \code{SPECIES}, \code{ALLNAMES}, \code{YEAR}, \code{NAFO_UNIT},
#'  \code{CATCH} (tonnes), and 9 landings groups. A value of \code{1} in a
#'  landings group column indicates that the corresponding \code{SPECIES}
#'  belongs to that group. Species may belong to more than one group. Species
#'  codes are the commercial landings codes.
#'@references Modified code from AC's ExtractIndicators/R/getLandings.R
#'@importFrom RODBC sqlQuery
#'@importFrom RODBC odbcConnect
#'@importFrom RODBC odbcClose
#'@family LAND functions
#'@export

extractLAND <- function(path, e.year) {
  
  print("running extractLAND()")
  
  # NAFO: 1968 - 1985
  NAFO <- function(area=paste("4VS","4VN","4X","4W", sep="','")) {
    y <- 1968:1985
    out <- list()
    
    for( i in 1:length(y)) {
      out[[i]] <- sqlQuery(channel,paste("select ",y[i]," year,description nafo_unit,species,sum(catch) catch 
                                         from comland.nafo_summary a,comland.nafo_area_codes r
                                         where a.area=r.area and upper(description) in ('",area,"') and year_of_activity = ",y[i],"
                                         group by ",y[i],",description, species",sep=""))
    }				
    
    dat <- do.call(rbind, out)
    return(dat)		
  }
  
  # ZIF: 1986 - 2002
  ZIF <- function (area=paste("4VS","4VN","4X","4W",sep="','")) {
    #match the zif data dictionary landings 
    y <- 1986:2002
    out <- list()
    for( i in 1:length(y)) {
      out[[i]]<-sqlQuery(channel,paste("select year,nafo_unit,zif2allcodes species,sum(wt) catch 
                                       from (select ",y[i]," year,b.id, nafod,unit,nafo_unit,species_code,wt 
                                       from (select catchers_recid||','||region_code||','||trip_num||','||sub_trip_num id,
                                       year_of_activity,nafo_division_code nafod ,nafo_unit_area unit,nafo_division_code||nafo_unit_area nafo_unit
                                       from cl.sub_trips_",y[i]," where nafo_division_code in ('",area,"')) a,
                                       (select catchers_recid||','||region_code||','||trip_num||','||sub_trip_num id, species_code, sum(live_wt) wt
                                       from cl.identified_catches_",y[i]," 
                                       group by catchers_recid||','||region_code||','||trip_num||','||sub_trip_num, species_code) b
                                       where b.id=a.id) a, 
                                       gomezc.indiseas_zif2allcodes b
                                       where 
                                       a.species_code=b.zif  
                                       group by year,nafo_unit, zif2allcodes",sep=""))
    }
    dat <- do.call(rbind,out)
    return(dat)		
  }
  
  # MARFIS: 2003 - e.year
  MARFIS <- function(area=paste("4V","4X","4W",sep="','"), e.year.marfis) {
    #match the vdc marfis landings
    
    data <- sqlQuery(channel, paste("select year_fished year,unit_area nafo_unit, marfis2allcodes species,sum(round(rnd_weight_kgs/1000,4)) catch
                                    from mfd_obfmi.marfis_catch_effort d, gomezc.indiseas_marfis2allcodes a
                                    where upper(nafo_div) in ('",area,"') and d.species_code=a.marfis
                                    and year_fished between '2003' and ",e.year.marfis,"
                                    and category_desc not in 'OTHER'
                                    GROUP BY year_fished ,unit_area, marfis2allcodes
                                    ;",sep=""))
    
    odbcClose(channel)
    
    return(data)
  }
  
  # call each function defined above
  ndat <- NAFO()
  zdat <- ZIF()
  mdat <- MARFIS(e.year.marfis = e.year)
  
  # format data
  nam <- commercial_groups # commercial_groups is a data object saved in the package. Type ?commercial_groups for more info
  names(nam)[1] <-'SPECIES'
  dat <- rbind(ndat, zdat, mdat)
  dat<-dat[with(dat,order(SPECIES,YEAR)),]
  dat<- dat[!dat$SPECIES %in% c(164,368,399,400,420,422,460,502,542,589,623,889,900,901,902,905,906,908,922,923,927,956,199),] #species with very sparse info and small landings
  dat$SPECIES <- ifelse(dat$SPECIES %in% c(512,514,516,518,520,525,529),529,dat$SPECIES) #combine the clams
  dat$SPECIES <- ifelse(dat$SPECIES %in% c(562,564),564,dat$SPECIES) #combined the whelks and periwinkles
  land<- merge(nam,dat)
  land$NAFO_UNIT <- toupper(land$NAFO_UNIT)
  land[land$NAFO_UNIT=='4VSA','NAFO_UNIT'] <- '4VSB'
  land[land$NAFO_UNIT=='4VSS','NAFO_UNIT'] <- '4VSU'
  land[land$NAFO_UNIT=='4VNN','NAFO_UNIT'] <- '4VN'
  
  # export data
  fp = file.path(path,'data','landings')
  dir.create(fp, recursive=T, showWarnings=F)
  landings <- land
  # landings <- landings/1000 # this will convert from tonnes to kg so that landings and RV survey are in same units
  save(landings, file = file.path(fp, "landings.RData"))
 
}   

