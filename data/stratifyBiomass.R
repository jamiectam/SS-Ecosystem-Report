#'@title Stratifies biomass and abundance estimates
#'@description Stratifies the biomass and abundance estimates based on the
#'  weights in \code{strat_weights} (type \code{?strat_weights} into the console
#'  for more information). User can choose whether to return length-based and/or
#'  q-corrected values.
#'@inheritParams biomassData
#'@param path Filepath indicating where to create folder to store the stratified
#'  data.
#'@param areas Areas (spatial scales) for which to stratify the fishery
#'  independent survey data. A separate dataframe will be exported for each
#'  area. Options are "shelf", "esswss", "nafo", "strat", or any combination of
#'  the four. Default is \code{areas = c("shelf", "esswss", "nafo", "strat")}.
#'@param lengthbased Logical parameter indicating whether to return stratified
#'  length-based data (\code{lengthBased = TRUE}) or not length-based data
#'  (\code{lengthBased = FALSE}).
#'@param qadjusted Logical parameter indicating whether to return q-adjusted
#'  biomass and abundance data (\code{qadjusted = TRUE}) or not q-adjusted
#'  biomass and abundance data (\code{qadjusted = FALSE}). q-adjustments are
#'  applied at the set level.
#'@return This function creates a directory to store the stratified data:
#'  data/stratified/folder/. The name of "folder" depends on the arguments
#'  \code{lengthbased} and \code{qadjusted}:
#'
#'  \code{if(lengthBased & qadjusted)}, then folder is named "lengthbased_qadj"
#'
#'  \code{if(lengthBased & !qadjusted)}, then folder is named
#'  "lengthbased_notqadj"
#'
#'  \code{if(!lengthBased & qadjusted)}, then folder is named
#'  "notlengthbased_qadj"
#'
#'  \code{if(!lengthBased & qadjusted)}, then folder is named
#'  "notlengthbased_notqadj"
#'
#'  Inside "folder" is a folder for each entry in \code{area}. Within each area
#'  folder is a .RData file for each year from s.year to e.year called
#'  year.RData (object name \code{out}).
#'
#'  \code{out} has 5 columns if \code{lengthbased = TRUE}: \code{SPECIES} ,
#'  \code{ID}, \code{FLEN}, \code{BIOMASS}, and \code{ABUNDANCE}.
#'
#'  \code{out} has 4 columns if \code{lengthbased = FALSE}: \code{SPECIES},
#'  \code{ID}, \code{BIOMASS}, and \code{ABUNDANCE}.
#'@references Modified code from AC's ExtractIndicators/R/stratifyBiomass.R
#'@family RV functions
#'@export


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