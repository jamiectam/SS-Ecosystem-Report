```{r, echo=FALSE}
# set up SSdata package functions
biomassData <- function(path, s.strat = 440, e.strat = 495, s.year, e.year,
                        vessel.correction = TRUE) {
  
  # display error message if s.strat or e.start are < 440 or > 495
  if(s.strat < 440 | s.strat > 495 | e.strat < 440 | e.strat > 495) {
    stop("error: s.strat or e.strat is outside of the Scotian Shelf")
  }
  
  # vector of years of interest
  yr <- s.year:e.year
  
  # start loop over years
  for(i in 1:length(yr)) {            
    
    # At Length ---------------------------------------------------------------
    
    # Table with "SPECIES", "FUNGROUP", "Q", "LENCORR" (88 observations)
    catch_coefs <- sqlQuery(channel, paste('select * from gomezc.indiseas_catchability_coeffs'))
    catch_coefs[catch_coefs == 'NULL'] <- NA
    outputs <- list()
    m = 0
    
    # Extract table with 12 columns: 
    ## YEAR, STRAT, MISSION, YYDDMMSS, XDDMMSS, SETNO, FLEN, 
    ## ABUNDANCE, QABUNDANCE, QBIOMASS, BIOMASS, SPECIES
    # QABUNDANCE and QBIOMASS are q-adjusted at the set level (BEFORE stratification). 
    for (j in 1:nrow(catch_coefs)) {
      outputs[[j]] <- qBiomass(species = catch_coefs[j, 1], fun_group = catch_coefs[j, 2], 
                               q = catch_coefs[j, 3], len_corr = catch_coefs[j, 4], year = yr[i]) 
    }
    out <- as.data.frame(do.call(rbind, outputs))
    out <- out[, -which(names(out) %in% c('YDDMMSS', 'XDDMMSS'))]  # remove lat/long columns
    
    if(vessel.correction) out <- vesselCorr(out)                   # apply vessel correction
    
    # Export length based data (biomass and abundance at 1 cm intervals)
    fna <- paste(path,"/data/length/",sep="")
    dir.create(fna, recursive = T, showWarnings = F)
    fna <- paste(fna, "num_biom_at_length", yr[i], ".RData", sep = "")
    save(out, file = fna, compress = T)
    
    # End of At Length ------------------------------------------------------
    
    # Begin aggregate ---------------------------------------------------------
    # Sum over length (aggregated biomass & abundance of species that have length-based q)
    ag.out <- aggregate(cbind(QBIOMASS, BIOMASS, QABUNDANCE, ABUNDANCE) ~ 
                          YEAR + STRAT + MISSION + SETNO + SPECIES, 
                        data = out, FUN = sum)
    
    # Extract table with 7 columns:
    ## MISSION, SETNO, STRAT, YEAR, SPEC, ABUNDANCE, BIOMASS
    dat <- sqlQuery(channel,paste("select distinct i.mission,i.setno,i.strat, to_char(sdate,'yyyy') year, spec,sum(nvl(totno,0)*1.75/i.dist) Abundance,sum(nvl(totwgt,0)*1.75/i.dist) biomass from 
							groundfish.gsinf i, groundfish.gscat c, mfd_stomach.nafo_strat sg where i.mission=c.mission and i.setno=c.setno and i.strat=sg.strat and to_char(sdate,'mm') in ('06','07','08') and
							i.strat between '",s.strat,"' and '",e.strat,"' and type=1 and spec<9000 and to_char(sdate,'yyyy')=",yr[i],"
							group by i.mission,i.setno,i.strat,slat , slong ,to_char(sdate,'yyyy'), spec;",sep=""))
    
    # Fix high biomass estimate(s) for dogfish
    dat[dat$BIOMASS > 8000 & dat$SPEC != 220,'BIOMASS'] <- 0 # added in November 08, 2013 
    
    # Change NA for abundance or biomass based on mean size of animals captured (for invs only)
    dat[is.na(dat$ABUNDANCE), 'ABUNDANCE'] <- 0
    dat[is.na(dat$BIOMASS), 'BIOMASS'] <- 0
    
    # Handle the missing biomass or abudance data where the other is present
    # by using mean weight = biomass/abundance
    if(any(dat$ABUNDANCE == 0 | dat$BIOMASS==0)) {
      
      wt <- sqlQuery(channel, paste("select * from mean_wts_for_fill_in;",sep="")) # table of mean fish weight ("SPEC" and "MEAN_WT_FISH", 640 observations)
      dat <- merge(dat, wt, by = c('SPEC'), all.x = T)                             # merge dat and wt
      dat[is.na(dat$MEAN_WT_FISH), 'MEAN_WT_FISH'] <- 0.5                          # if mean fish weight is NA, assign mean weight of 0.5
      
      if(any(dat$ABUNDANCE == 0)) {
        dat$ABUNDANCE[dat$ABUNDANCE == 0] <- dat$BIOMASS[dat$ABUNDANCE==0]/dat$MEAN_WT_FISH[dat$ABUNDANCE==0] # abundance = biomass/mean weight  			
      } # DD added this bracket so that biomass correction is OUTSIDE of the if statement for abundance
      if(any(dat$BIOMASS == 0)) {
        dat$BIOMASS[dat$BIOMASS == 0] <- dat$ABUNDANCE[dat$BIOMASS==0]*dat$MEAN_WT_FISH[dat$BIOMASS==0]					
      }
    }
    
    # Change ABUNDANCE = inf to ABUNDANCE = 1 
    if(unique(dat$YEAR)==2012 & any(!is.finite(dat$ABUNDANCE))) dat[which(!is.finite(dat$ABUNDANCE)),'ABUNDANCE'] <- 1
    
    # If any biomass or abundance is still NA, stop code
    if(any(is.na(dat[,c('BIOMASS','ABUNDANCE')]))) browser()
    
    # Replace 0's with small values
    dat[dat$BIOMASS == 0 ,'BIOMASS'] <- 0.01	
    dat[dat$ABUNDANCE == 0 ,'ABUNDANCE'] <- 1
    
    dat$QABUNDANCE <- dat$ABUNDANCE 
    dat$QBIOMASS <- dat$BIOMASS
    dat$SPECIES <- dat$SPEC       # need this line? Drop SPECIES column at line 145
    if(vessel.correction) dat <- vesselCorr(dat) # Apply vessel correction
    dat <- dat[, -which(names(dat) == 'SPECIES')]
    
    # Add in zero sets. 
    # These are needed to calculate the stratified means
    extra.sets <- sqlQuery(channel,paste("select distinct i.mission,i.setno,i.strat,to_char(sdate,'yyyy') year from groundfish.gsinf i, nafo_strat sg where 
			i.strat=sg.strat and to_char(sdate,'yyyy') =",yr[i]," and to_char(sdate,'mm') in ('06','07','08') and i.strat between '",s.strat,"' and '",e.strat,"' and type=1;",sep=""))
    s <- unique(dat$SPEC)	
    m <- matrix(0, nrow = dim(extra.sets)[1], ncol=length(s),
                dimnames = list(c(1:nrow(extra.sets)), s))	
    m <- cbind(extra.sets, m)
    h <- melt(m, id = c("MISSION", "SETNO", "STRAT", "YEAR"))
    names(h)[which(names(h) %in% c('variable','value'))] <- c('SPECIES','BIOMASS')
    h$QBIOMASS <- h$QABUNDANCE <- h$MEAN_WT_FISH <- h$ABUNDANCE <- 0
    names(dat)[1] <- 'SPECIES'
    
    l <- rbind(dat, h, all = T) # not sure why need all = T. This adds an extra
    # row where MISSION = TRUE and SPECIES = TRUE (other columns = 1).
    # I remove this row below
    
    dat <- l[!duplicated(l[,c('MISSION','SETNO','SPECIES')], fromLast = F),]
    dat <- dat[-which(dat$MISSION == TRUE), ] # remove row where MISSION == TRUE
    
    # Add in the qadjusted data
    dat1 <- dat
    w <- merge(dat, ag.out, by = c('MISSION','SETNO','STRAT','SPECIES','YEAR'), all.x = T)	
    # Add columns for ABUNDANCE, QABUNDANCE, BIOMASS, and QBIOMASS
    w$ABUNDANCE <- 0		
    w$QABUNDANCE <- 0		
    w$BIOMASS <- 0		
    w$QBIOMASS <- 0		
    
    # Fill in columns for ABUNDANCE, QABUNDANCE, BIOMASS, and QBIOMASS with appropriate data
    # .x is from dat (not q adj); .y is from ag.out (qadj)
    w$ABUNDANCE <- ifelse(is.na(w$ABUNDANCE.y), w$ABUNDANCE.x, w$ABUNDANCE.y)
    w$BIOMASS <- ifelse(is.na(w$BIOMASS.y), w$BIOMASS.x, w$BIOMASS.y)
    w$QABUNDANCE <- ifelse(is.na(w$QABUNDANCE.y), w$QABUNDANCE.x, w$QABUNDANCE.y)
    w$QBIOMASS <- ifelse(is.na(w$QBIOMASS.y), w$QBIOMASS.x, w$QBIOMASS.y)
    
    # Export aggregated data
    dat <- w[, c('MISSION','SETNO','SPECIES','YEAR','STRAT','ABUNDANCE','BIOMASS','QABUNDANCE','QBIOMASS')]		
    fna <- paste(path,"/data/aggregate/",sep="")
    dir.create(fna, recursive = T, showWarnings = F)
    save(dat, file = paste(fna, "num_biom", yr[i], ".RData", sep=""))
    
    rm(dat, w)
    
    print(paste("biomassData() function: finished year", yr[i]))
  } # end of loop over years
  
} # end of function

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

extractLW <- function(path, s.year, e.year) {
  
  # Create file path to store data
  fna <- file.path(path, "data", "lenwgt")  
  dir.create(fna, recursive = T, showWarnings = F)
  
  yr <- s.year:e.year
  
  for(i in 1:length(yr)){
    
    # Extract data from SQL
    wt <- sqlQuery(channel,paste("select distinct strat,spec species,flen,fwt from groundfish.gsinf i, groundfish.gsdet d where i.mission=d.mission and i.setno=d.setno and to_char(sdate,'yyyy') = ",yr[i]," and to_char(sdate,'mm') in ('06','07','08') and strat between '440' and '495' and fwt is not null and flen is not null;",sep=""))
    
    # Save data
    save(wt, file = paste(fna, "/lw", yr[i], ".Rdata",sep=""), compress = T)
  }
  
}

extractRV <- function(path, s.year, e.year, lengthbased, qadjusted,
                      areas = c("shelf", "esswss", "nafo", "strat")){
  
  # Extract biomass & abundance data
  biomassData(path = path, s.year = s.year, e.year = e.year)
  
  # Stratify biomass & abundance data
  stratifyBiomass(path = path, s.year = s.year, e.year = e.year, 
                  lengthbased = lengthbased, qadjusted = qadjusted, areas = areas)
  
  print("Biomass data extracted & stratified")
}

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

plot_LW <- function(species, year){
  
  area='4VWX'
  
  bb<-sqlQuery(channel,paste("select ",species,",flen,fwt from groundfish.gsdet d, groundfish.gsinf i where i.mission=d.mission and i.setno=d.setno and i.strat in (select strat from 
                           mflib.gsmgt where unit in '",area,"') and to_char(sdate,'mm') in ('06','07','08') AND to_char(sdate,'yyyy') in ('",year,"') and spec in (",species,") and fwt is not null;",sep=""))
  
  if(nrow(bb)<100) {
    bb<-sqlQuery(channel,paste("select ",species,",flen,fwt from groundfish.gsdet d, groundfish.gsinf i where i.mission=d.mission and i.setno=d.setno and i.strat in (select strat from 
                                             mflib.gsmgt where unit in '",area,"') and to_char(sdate,'mm') in ('06','07','08') AND to_char(sdate,'yyyy') between ('",year,"'-5) and ('",year,"'+5) and spec in (",species,") and fwt is not null;",sep=""))
  }
  if(nrow(bb)<100){
    bb<-sqlQuery(channel,paste("select ",species,",flen,fwt from groundfish.gsdet d, groundfish.gsinf i where i.mission=d.mission and i.setno=d.setno and i.strat in (select strat from 
                                            mflib.gsmgt where unit in '",area,"') and to_char(sdate,'mm') in ('06','07','08') and spec in (",species,") and fwt is not null;",sep=""))
  }
  
  out <- with(bb, lm(log(FWT)~log(FLEN), bb))
  a.est <- coef(out)[1]
  b.est <- coef(out)[2]
  
  
  with(bb, plot(FLEN, FWT))
  curve(exp(a.est+log(x)*b.est),add=T,col='red')
  title(paste(species,"-",year,"-",area),
        sub = (paste("a.est =", signif(a.est, 3), "& b.est =", signif(b.est, 3))))
  
}

qBiomass <- function(species, year, fun_group = NA, q = 0, len_corr = 1) {
  
  area = "4VWX"
  
  # make a dummy table of lengths (this is to deal with the herring switch from cm to mm)
  # it makes the code WAY slower. Leaving for now, but there may be a way to speed up!
  suppressWarnings(sqlQuery(channel,paste("DROP table gs_len;")))
  flens = data.frame(CLASS=1, FLEN=1:500)
  sqlSave(channel,dat=flens, tablename='GS_LEN',rownames=F)
  
  # Check to see if there is enough data to bother with this species
  initial <- sqlQuery(channel,paste("select count(distinct c.mission||','||c.setno||','||",species,"||','||fshno) from groundfish.gsdet c, groundfish.gsinf i where i.mission=c.mission and i.setno=c.setno and to_char(sdate,'mm') in ('06','07','08') and strat in (select distinct strat from mflib.gsmgt where unit in ('",area,"')) and spec=",species,";",sep=""))
  go <- ifelse(initial>30,TRUE,FALSE)
  if(go==TRUE) {
    
    # Obtain the numbers at length for 1 cm length interval
    aa<-sqlQuery(channel, paste("select distinct year,strat,mission,YDDMMSS,XDDMMSS,setno,flen,SUM(clen * DECODE(NVL(totwgt,0),0,1,DECODE(NVL(sampwgt,0),0,1,totwgt/sampwgt))) catch
	  			from 
	 				(SELECT distinct year, icl.strat, icl.mission, icl.setno, icl.flen FLEN,sampwgt,totwgt,nvl(clen,0)*1.75/dist clen, YDDMMSS,XDDMMSS
						 FROM
	 						(SELECT mission,setno, flen, SUM(clen) clen, AVG(fwt) avg_fwt
	    						FROM   groundfish.gsdet
	    						WHERE flen IS NOT NULL AND spec=",species,"
	    						GROUP BY mission,setno, FLEN) d,
	 						(SELECT year, mission, setno, strat, dist, totwgt, sampwgt, flen,YDDMMSS,XDDMMSS
	     						FROM
	       							(SELECT flen
	        							FROM GS_LEN
	        							WHERE class=1
	          							AND flen <=(SELECT max(flen) + 1
	              										fROM groundfish.gsdet
	              										WHERE spec=",species," AND flen IS NOT NULL
	                									AND (mission, setno) IN
	                   										(SELECT DISTINCT i.mission, i.setno
	                     										FROM groundfish.gsinf i
	                     										WHERE to_char(sdate,'mm') in ('06','07','08') 
	                      										AND to_char(sdate,'yyyy') in (",year,") AND i.type=1))) l,
	       							(SELECT year, i.mission, i.setno, strat, dist, totwgt,YDDMMSS,XDDMMSS, sampwgt
	        							FROM
	        								(SELECT mission, setno, totwgt, sampwgt
	                							FROM groundfish.gscat WHERE spec=",species,") c,
	          								(SELECT to_char(sdate,'yyyy') year, i.mission, i.setno,slat YDDMMSS,slong*-1 XDDMMSS, i.strat, dist
	                							FROM groundfish.gsinf i
	           									WHERE  to_char(sdate,'mm') in ('06','07','08')
	                      						AND to_char(sdate,'yyyy') in (",year,") AND i.type=1) i
	        					WHERE i.mission=c.mission(+)
	          					AND i.setno=c.setno(+)) ic) icl
	 			WHERE icl.mission=d.mission(+) AND icl.setno=d.setno(+) AND icl.flen=d.flen(+) and icl.strat in (select distinct strat from mflib.gsmgt where unit in ('",area,"')))
	 			 group by year,strat,mission,setno,flen, XDDMMSS,YDDMMSS;",sep=""	))
    
    run <- ifelse(nrow(aa) < 1, FALSE, TRUE)
    
    if(run==TRUE) {
      
      #order the data by fork length
      aa <- aa[order(aa$FLEN),]
      
      # Not length-based q for some species (e.g., herring,)
      if(is.na(fun_group) && q > 0) {
        aa$q <- q
        aa$corCatch <- aa$CATCH/aa$q
        aa$model<-paste('single_q-', q)
      } 
      
      # adjust the catch at lengths by appropriate q and length correction
      
      # cod length-based q
      if(fun_group == 'cod' && q == 0) {
        a1 = -5.14
        b1 = 0.141
        g1 = 0.870
        aa$q <- g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch <- aa$CATCH/aa$q
        aa$model <- 'cod'
      }
      
      # haddock length-based q
      if(fun_group == 'haddock'&& q == 0) {	
        a1=-2.80
        b1 = 0.0661
        g1 = 1.5
        aa$q <- g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch <- aa$CATCH/aa$q
        aa$model <- 'haddock'
      }
      
      # pelagic gadoids length-based q
      if(fun_group == 'pg'&& q == 0 ) {
        a1 = -4.61
        b1 = 0.0789
        g1 = 0.58
        aa$q <-g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch<-aa$CATCH/aa$q
        aa$model<-'pg'
      }
      
      # demersal gadoids length-based q
      if(fun_group == 'dg'&& q == 0) {
        a1=-3.50
        b1=0.0925
        g1=0.968	
        aa$q<-g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch<-aa$CATCH/aa$q
        aa$model<-'dg'
      }
      
      # flatfish length-based q
      if(fun_group == 'ff'&& q == 0) {
        a1=-4.35
        b1=0.111
        g1=0.831	
        aa$q<-g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch<-aa$CATCH/aa$q
        aa$model<-'ff'
      }
      
      # ling length-based q	 	
      if(fun_group == 'ling' && q == 0) {
        a1=-13.9
        b1=0.193
        g1=1.66	
        aa$q<-g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch<-aa$CATCH/aa$q
        aa$model<-'ling'
      }
      
      # small pelagics length-based q from H.Benoit
      if(fun_group == 'sp' && q == 0) {
        a1=-17.7165
        b1=0.175
        g1=37710.7	
        aa$q<-g1*(exp(a1+b1*(aa$FLEN*len_corr)))/(1+exp(a1+b1*(aa$FLEN*len_corr)))
        aa$corCatch<-aa$CATCH/aa$q
        aa$model<-'sp'
      }
      
      
      # calculate appropriate length weight parameters (a's and b's)
      
      bb<-sqlQuery(channel,paste("select ",species,",flen,fwt from groundfish.gsdet d, groundfish.gsinf i where i.mission=d.mission and i.setno=d.setno and i.strat in (select strat from 
	 	mflib.gsmgt where unit in '",area,"') and to_char(sdate,'mm') in ('06','07','08') AND to_char(sdate,'yyyy') in ('",year,"') and spec in (",species,") and fwt is not null;",sep=""))
      
      if(nrow(bb)<100) {bb<-sqlQuery(channel,paste("select ",species,",flen,fwt from groundfish.gsdet d, groundfish.gsinf i where i.mission=d.mission and i.setno=d.setno and i.strat in (select strat from 
	 	mflib.gsmgt where unit in '",area,"') and to_char(sdate,'mm') in ('06','07','08') AND to_char(sdate,'yyyy') between ('",year,"'-5) and ('",year,"'+5) and spec in (",species,") and fwt is not null;",sep=""))
      }
      if(nrow(bb)<100){bb<-sqlQuery(channel,paste("select ",species,",flen,fwt from groundfish.gsdet d, groundfish.gsinf i where i.mission=d.mission and i.setno=d.setno and i.strat in (select strat from 
	 	mflib.gsmgt where unit in '",area,"') and to_char(sdate,'mm') in ('06','07','08') and spec in (",species,") and fwt is not null;",sep=""))
      }
      
      out <- with(bb, lm(log(FWT)~log(FLEN), bb))
      a.est <- coef(out)[1]
      b.est <- coef(out)[2]
      
      #q-adj biomass estimate
      aa$est.wt <- exp(a.est+log(aa$FLEN)*b.est)
      aa$qBiomass <- aa$corCatch*aa$est.wt
      aa$rvBiomass <- aa$CATCH*aa$est.wt
      names(aa)[8:14] <- c('ABUNDANCE','Q','QABUNDANCE','MODEL','ESTWT','QBIOMASS','BIOMASS')
      
      aa <- aa[,c(1:8,10,13,14)]
      aa$SPECIES <- species
      aa$QBIOMASS <- aa$QBIOMASS/1000 #for kg
      aa$BIOMASS <- aa$BIOMASS/1000  #for kg
      return(aa)
    }
  }
}

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

stratifyBiomass <- function(path, s.year, e.year, lengthbased, qadjusted, 
                            areas = c('esswss','nafo','shelf','strat')	) {
  
  print("running stratifyBiomass()")
  
  # strata weights
  st.weights <- strat_weights
  
  out.fp <- file.path(path, "data", "stratified") # where to export data: path/data/stratified
  yr <- s.year:e.year                             # years to stratify  
  
  for(i in 1:length(yr)) {
    
    #trawlable units differnce (fanning 1985)
    if(yr[i]<=1981) st.weights$TUNITS 	<- st.weights$AREA/((35./6080.2)*1.75)
    if(yr[i]>1981) st.weights$TUNITS 	<- st.weights$AREA/((41./6080.2)*1.75)
    
    if(!lengthbased) {
      
      u <- 'notlengthbased'
      
      # Load data (exported from biomassData())
      fna 	<- paste(path,"/data/aggregate/num_biom",yr[i],".RData",sep="") # path to biomass data 
      # (path/data/aggregate/sumq/num_biom in AC's code)
      load(fna)                                                              # load the data into an object called dat
      dat$SPECIES <- as.numeric(dat$SPECIES)                                 # make sure species codes are numeric
      
      #Begin Stratified Estimates
      if(any(dat$YEAR<1970)) dat <- dat[-which(dat$YEAR<1970),]              # remove years before 1970
      
      if(qadjusted) {
        dat <- dat[,-which(names(dat) %in% c('BIOMASS','ABUNDANCE'))]
        names(dat)[which(names(dat) %in% c('QBIOMASS'))] <- c('BIOMASS')
        names(dat)[which(names(dat) %in% c('QABUNDANCE'))] <- c('ABUNDANCE')
        u <- paste(u,'_qadj',sep="") 
      }
      if(!qadjusted) {
        dat <- dat[,-which(names(dat) %in% c('QBIOMASS','QABUNDANCE'))]
        u <- paste(u,'_notqadj',sep="") 
      }
      
      #strata means	
      dat <- aggregate(cbind(BIOMASS,ABUNDANCE) ~ YEAR + STRAT + SPECIES, 
                       data = dat, FUN = mean)	
      
      #add in the strata weightings for the appropriate groups
      dat <- merge(dat, st.weights, by='STRAT', all.x = T)	
      out <- dat
      
      #strata totals				
      out$BIOMASS <- out$BIOMASS*out$TUNITS
      out$ABUNDANCE <- out$ABUNDANCE*out$TUNITS
      out$SPECIES <- as.numeric(out$SPECIES)
      ou <- out[out$BIOMASS!=0,]	                # remove the zeros (all species were assumed across all areas)	
      
      # assign area ID and export each year
      for(j in 1:length(areas)) {
        out <- defineAreas(ou, area = areas[j])
        out <- aggregate(cbind(BIOMASS,ABUNDANCE)~SPECIES+ID,data=out,FUN=sum)
        
        fp1 <- file.path(out.fp,u,areas[j])
        dir.create(fp1,recursive=T,showWarnings=F)
        save(out, file = file.path(fp1, paste(yr[i],".RData",sep="")))
      }
    } # end of not length based
    
    if(lengthbased) {
      
      u <- 'lengthbased'
      
      # load length based data  (exported from biomassData())
      fna <- paste(path,"/data/length/num_biom_at_length", yr[i], ".RData", sep="")  # path to length based data
      load(fna)                                                                       # load data into object called out
      
      # load not length based data to add in species without length info
      fna <- paste(path,"/data/aggregate/num_biom",yr[i],".RData",sep="")		        # this is the same as /data/aggregate/sumq/num_bio in AC's code
      load(fna)
      
      if(any(dat$YEAR<1970)) dat <- dat[-which(dat$YEAR<1970),]
      spo <- unique(out$SPECIES)
      dat <- dat[!dat$SPECIES %in% spo,]
      dat$FLEN <- -99
      
      out <- rbind(out,dat)
      
      dat <- out
      rm(out)
      
      #Make Stratified Estimates
      if(any(dat$YEAR<1970)) { dat <- dat[-which(dat$YEAR<1970),] }
      
      if(qadjusted) {
        dat <- dat[,-which(names(dat) %in% c('BIOMASS','ABUNDANCE'))]
        names(dat)[which(names(dat) %in% c('QBIOMASS'))] <- c('BIOMASS')
        names(dat)[which(names(dat) %in% c('QABUNDANCE'))] <- c('ABUNDANCE')
        u <- paste(u,'_qadj',sep="") 
      }
      
      if(!qadjusted) {
        dat <- dat[,-which(names(dat) %in% c('QBIOMASS','QABUNDANCE'))]
        u <- paste(u,'_notqadj',sep="") 
      }
      
      
      dat <- aggregate(cbind(BIOMASS,ABUNDANCE)~YEAR+STRAT+SPECIES+FLEN,data=dat,FUN=mean)	
      
      #add in the strata weightings for the appropriate groups
      dat <- merge(dat,st.weights,by='STRAT',all.x=T)	#combine to ensure only sampled strat are included
      
      dat$BIOMASS <- dat$TUNITS*dat$BIOMASS			
      dat$ABUNDANCE <- dat$TUNITS*dat$ABUNDANCE	
      dat$SPECIES <- as.numeric(dat$SPECIES)
      dat <- dat[dat$BIOMASS>0,]
      ou <- dat
      
      for(j in 1:length(areas)) {
        out <- defineAreas(ou, area = areas[j])
        out <- aggregate(cbind(BIOMASS,ABUNDANCE)~SPECIES+ID+FLEN,data=out,FUN=sum)
        fp1 <- file.path(out.fp,u,areas[j])
        dir.create(fp1,recursive=T,showWarnings=F)
        save(out,file=file.path(fp1,paste(yr[i],".RData",sep="")))
      }
      
    } # end of if(lengthbased)
    
  } # end of loop over years
} # end of function

vesselCorr <- function(x) {
  
  cf = NULL
  x$CFVESSEL = 1  # initialise
  
  # vessel change correction factors apply to these years:
  HAM=1   #  Lady Hammond (1979 - 1981)
  ATC=2   #  A.T. Cameron (1982 - 1983)
  
  # species codes used by the database
  cod=10
  haddock=11
  whitehake=12
  silverhake=14
  plaicelarge=40
  plaicesmall=40
  witch=41
  yellowtail=42
  winterflounder=43
  cf$cod[HAM]         = 0.8
  cf$haddock[HAM]     = 1.0
  cf$whitehake[HAM]   = 1.0
  cf$silverhake[HAM]  = 1.0
  cf$plaicesmall[HAM] = 1   # <=28cm
  cf$plaicelarge[HAM] = 1   # > 28cm
  cf$witch[HAM]       = 0.8
  cf$yellowtail[HAM]  = 0.8
  cf$winterflounder[HAM] = 1.0
  cf$cod[ATC]         = 0.8
  cf$haddock[ATC]     = 1.2
  cf$whitehake[ATC]   = 1.0
  cf$silverhake[ATC]  = 1.0
  cf$plaicesmall[ATC] = 0.7
  cf$plaicelarge[ATC] = 1.0
  cf$witch[ATC]       = 0.8
  cf$yellowtail[ATC]  = 0.8
  cf$winterflounder[ATC] = 1.0
  attach (x)
  x$CFVESSEL[ which( (substring(MISSION,1,3)=="HAM" & SPECIES==cod)) ] = cf$cod[HAM]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="HAM" & SPECIES==witch)) ] = cf$witch[HAM]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="HAM" & SPECIES==yellowtail)) ] = cf$yellowtail[HAM]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="ATC" & SPECIES==cod)) ] = cf$cod[ATC]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="ATC" & SPECIES==haddock)) ] = cf$haddock[ATC]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="ATC" & SPECIES==plaicesmall && FLEN<=28)) ] = cf$plaicesmall[ATC]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="ATC" & SPECIES==witch)) ] = cf$witch[ATC]
  x$CFVESSEL[ which((substring(MISSION,1,3)=="ATC" & SPECIES==yellowtail)) ] = cf$yellowtail[ATC]
  detach (x)
  x$BIOMASS = x$BIOMASS * x$CFVESSEL #similar coefficients in manning 1985      
  x$ABUNDANCE = x$ABUNDANCE * x$CFVESSEL
  if(any(names(x) %in% c('QBIOMASS','QABUNDANCE'))) {
    x$QBIOMASS = x$QBIOMASS * x$CFVESSEL #similar coefficients in manning 1985      
    x$QABUNDANCE = x$QABUNDANCE * x$CFVESSEL
  }
  x <- x[,-which(names(x)=='CFVESSEL')]
  return (x)
}


```

