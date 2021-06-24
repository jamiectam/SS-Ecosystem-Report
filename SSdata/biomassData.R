#'@title Extracts and exports fishery independent biomass and abundance data
#'@description Extracts and exports annual q-adjusted and not q-adjusted,
#'  length-based (separated into 1 cm length classes) and not length-based
#'  (aggregated over length) biomass and abundance data from the Scotian Shelf
#'  summer research vessel surveys. q-adjustments are applied at the set level
#'  (before stratification).
#'@details User must define \code{channel = odbcConnect("ptran", uid = ###, pwd
#'  = ###)} in the global environment. This channel must have access to the
#'  gscat, gsdet, gsinf, and gs_lengths tables from the groundfish database and
#'  the nafo_strat table from the mfd_stomach database.
#'
#'  Units: biomass - kg; abundance - numbers; length - cm.
#'
#'  If \code{BIOMASS} or \code{ABUNDANCE} data is missing when the other is
#'  present, the missing field is filled in based on the relationship
#'  \eqn{Weight_{Avg} = \code{BIOMASS}/\code{ABUNDANCE}}. \eqn{Weight_{Avg}} is
#'  from a sql call to "mean_wts_for_fill_in." If \eqn{Weight_{Avg}} is not
#'  available for a species, it is assumed that \eqn{Weight_{Avg} = 0.5}.
#'
#'  Abnormally high dogfish biomass estimates (> 8000 kg at the set level) are
#'  replaced with \code{BIOMASS = 0} kg.
#'
#'  Abundance estimates of \code{Inf} in 2012 are replaced with \code{ABUNDANCE
#'  = 1}.
#'
#'  Code will stop if any \code{BIOMASS} or \code{ABUNDANCE} estimates are
#'  \code{NA}.
#'
#'  Estimates of zero are replaced with small values: \code{BIOMASS = 0} is
#'  replaced with \code{BIOMASS = 0.01}; \code{ABUNDANCE = 0} is replaced with
#'  \code{ABUNDANCE = 1}.
#'@param path Filepath indicating where to create folders to store the extracted
#'  data.
#'@param s.strat Stratum for which to begin data extraction. Default is
#'  \code{s.strat = 440}. Code will stop with an error message if \code{s.strat
#'  < 440} or \code{s.strat > 495}.
#'@param e.strat Stratum for which to end data extraction. Default is
#'  \code{e.strat = 495}.  Code will stop with an error message if \code{e.strat
#'  < 440} or \code{e.strat > 495}.
#'@param s.year Year for which to begin data extraction.
#'@param e.year Year for which to end data extraction.
#'@param vessel.correction Logical value indicating whether to apply vessel
#'  correction to \code{BIOMASS} and \code{ABUNDANCE}. Default is
#'  \code{vessel.correction} = TRUE.
#'@return Creates directories to store extracted data.
#'
#'  Not length-based data are stored in path/data/aggregate/. This folder
#'  includes an RData file for each year called year.RData (object name
#'  \code{dat}). \code{dat} has 9 columns: \code{MISSION}, \code{SETNO},
#'  \code{SPECIES}, \code{YEAR}, \code{STRAT}, \code{BIOMASS}, \code{ABUNDANCE},
#'  \code{QBIOMASS} and \code{QABUNDANCE}.
#'
#'  Length-based data is stored in path/data/length/. This folder includes an
#'  RData file for each year called num_biom_at_length_year.RData (object name
#'  \code{out}). \code{out} has 10 columns: \code{MISSION}, \code{SETNO},
#'  \code{SPECIES}, \code{YEAR}, \code{STRAT}, \code{BIOMASS}, \code{ABUNDANCE},
#'  \code{QBIOMASS} and \code{QABUNDANCE}, and \code{FLEN}, where \code{FLEN} is
#'  length in 1 cm increments.
#'
#'@references Modified code from AC's ExtractIndicators/R/biomassData.R
#'@importFrom stats aggregate
#'@importFrom reshape melt
#'@importFrom RODBC sqlQuery
#'@family RV functions
#'@export

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


