#'@title Plot length-weight data
#'@description Plots length-weight data for given species and year. Displays the
#'  estimated values of a and b.
#'@param year Year for which to extract data. Must be one of the 88 species in
#'  the \code{catch_coefs} table in \code{biomassData()}.
#'@param species Species code for which to extract data.
#'@return A plot of weight (g) vs. length (cm).
#'@references Modified from AC's code in script qBiomass.
#'@importFrom utils write.csv
#'@family LW functions
#'@export

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
