#'@title Applies q-adjustment to biomass and abundance data at the set level
#'@description Applies q-adjustments to biomass and abundance data at the set
#'  level and returns to the results to \code{biomassData()}. Length-based
#'  adjustments are applied for several functional groups; not length-based
#'  adjustments are applied for some species; no adjustment is applied for
#'  remaining species.
#'
#'@details \code{qBiomass()} is called by the \code{biomassData()} function.
#'  Arguments are supplied to \code{qBiomass()} from the \code{catch_coefs}
#'  table (SQL call to gomezc.indiseas_catchability_coeffs), which includes
#'  columns \code{FUNGROUP}, \code{Q}, and \code{LENCORR} for 88 \code{SPECIES}.
#'  \code{qBiomass()} loops over each row of this table for each year supplied
#'  to \code{biomassData()}.
#'
#'  \code{qBiomass()} obtains numbers (\code{ABUNDANCE}) at length for the given
#'  \code{species} for 1 cm intervals from the groundfish database.
#'
#'  Biomass estimates are based on the length-weight relationship:
#'  \deqn{ln(Weight) = a + ln(Length) * b} The function obtains length and
#'  weight data for \code{species} in \code{year} from the groundfish database
#'  and estimates \eqn{a} and \eqn{b}. The average weight of fish at each length
#'  is then estimated using: \deqn{Weight_{Avg} = exp(a + ln(Length)*b)} and the
#'  biomass estimate is: \deqn{\code{BIOMASS} = \code{ABUNDANCE} *
#'  Weight_{Avg}}
#'
#'  If there are less than 100 length-weight observations in the database for
#'  the given \code{species} in \code{year}, the function will include
#'  observations from the 5 years before and after \code{year}. If there are
#'  still less than 100 observations, data from all years will be included. This
#'  can result is some discrepancies when biomass data is extracted in different
#'  years.
#'
#'  \bold{q-correction}
#'
#'  q is calculated based on the arguments \code{fun_group} and \code{q}.
#'
#'  If \code{q} > \code{0} and \code{fun_group} = \code{NA}, then a single q is applied for
#'  all length classes. Otherwise, \code{q} is calculated using the equation:
#'  \deqn{q =
#'  g1*(exp(a1+b1*(\code{LENGTH}*\code{lencorr})))/(1+exp(a1+b1*(\code{LENGTH}*\code{lencorr})))}
#'   where \eqn{g1}, \eqn{a1}, and \eqn{b1} are hard-coded for each
#'  \code{fun_group}, \code{len_corr} is an argument provided by
#'  \code{biomassData()}, and \code{LENGTH} is the length interval (cm).
#'
#'  q-corrected abundance: \eqn{\code{QABUNDANCE} = \code{ABUNDANCE}/q}.
#'
#'  q-corrected biomass: \eqn{\code{QBIOMASS} = \code{QABUNDANCE} *
#'  Weight_{Avg}}
#'
#'@param species Species code from the \code{catch_coeffs} table. This argument
#'  is supplied by \code{biomassData()}.
#'@param year Year to calculate the data, supplied by \code{biomasData()}
#'@param fun_group Functional group from the \code{catch_coeffs} table. This
#'  argument is supplied by \code{biomassData()}. If \code{fun_group} = \code{NA}, then
#'  this \code{q} values is applied to all length classes.
#'@param q Initial q-correction from the \code{catch_coeffs} table. This
#'  argument is supplied by \code{biomassData()}. If \code{q} > \code{0} and
#'  \code{fun_group} = \code{NA}, then this \code{q} value is applied to all length
#'  classes.
#'@param len_corr Initial length correction from the \code{catch_coeffs} table.
#'  This argument is supplied by \code{biomassData()}.
#'@return Returns the not q-adjusted and q-adjusted length-based biomass and
#'  abundance to \code{biomassData()} for species that have q-adjustments.
#'@references Modified code from AC's ExtractIndicators/R/qBiomass (function
#'  \code{biomass_q_adj()})
#'@importFrom stats coef
#'@importFrom graphics curve
#'@importFrom graphics title
#'@importFrom RODBC sqlQuery
#'@family RV functions
#'@export

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

