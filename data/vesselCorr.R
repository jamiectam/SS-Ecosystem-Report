#'@title Applies vessel correction factor to key species
#'@description Called by \code{biomassData()} to apply vessel correction factor
#'  to key species.
#'  
#'@param x Fishery indepdendent survey data as extracted by \code{biomassData()}.
#'@return This function returns the corrected biomass and abundance data to
#'  function \code{biomassData()}.
#'@family RV functions
#'@references Modified code from AC's ExtractIndicators/R/vesselCorrection.R.


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
