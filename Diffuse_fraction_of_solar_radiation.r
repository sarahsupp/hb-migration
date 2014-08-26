# Diffuse_fraction_of_solar_radiation.R
#
# Pieter Beck (psabeck@gmail.com)
# 09-12-2013
# Modified by S.R. Supp 2014 (sarah@weecology.org)
# 
# functions to calculate what fraction of incoming solar radiation reaches the surface as diffuse radiation
#
# Changelog:
#   VERS  | DATE  		 | CHANGES				| BY
#   ------|------------|----------------|---- 
#		1.0.0	|	09-12-2013 | Wrote script		| PB

#Id is diffuse radiation
#I is 'global radiation' (DSWRF) (for a horizontal plane, global irradiance is the sum of the diffuse and the direct component)
#Id/I (IdoverI) is the fraction of diffuse radiation in global radiation 
#IdoverI is estimated as a function of kt (Liu and Jordan 1960) (and in some cases solar zenith angle)

#kt is the clearness index
#The clearness index measures the proportion of horizontal extraterrestrial radiation (Io)
#reaching the surface.
#it is defined as kt = I / (Io cos solarzen).
#DSWRF = Downward short wave radiation flux 

IdoverI.calc <- function(kt,solarzen){
  #calculate the fraction of diffuse radiation in global radiation
  #using the Reindl* method (Helbig 2009)
  # Args: 
  #  kt: The clearness index (sensu Liu and Jordan 1960)
  #  solarzen: solar zenith angle in radians
  # Returns:
  #   the fraction of diffuse radiation (Id) in global radiation (I)
  
  if (kt <= 0.3) {IdoverI <- 0.1020 - 0.248*kt}else{
  if(kt < 0.78){
    solar.elev <- pi/2 - solarzen
    IdoverI <- 1.4 - 1.749*kt + 0.177*sin(solar.elev)}else{
    IdoverI <- 0.147}
  }
  return(IdoverI)
}

SpSd.calc <- function(Rsurface,R_extra_terr,solarzen){
  #calculate direct (Sp) and diffuse (Sd) radiation from global radiation
  #and horizontal extraterrestrial radiation 
  # Args: 
  #  Rsurface: incoming shortwave radiation at the surface
  #  R_extra_terr: horizontal extraterrestrial radiation 
  #  solarzen: solar zenith angle in radians
  # Returns:
  #  A two-column matrix giving incoming direct [,1], and diffuse [,2] radiation at the surface
  
  if (solarzen > 2*pi){cat("STOP STOP STOP provide solarzenith in radiance to SpSd.calc\n");browser()}
  kt <- Rsurface/R_extra_terr
  #TODO: Added this check (e.g. if Rsurface and R_extra_terr == 0, kt == NaN which is a problem)
  if (is.na(kt)){
    kt = 0
  }
  #the formula in Lanini 2010 p1 is
  #kt <- Rsurface/(Io*cos(solarzen))
  #but this is for Io*cos(solarzen) being the horizontal extraterrestrial radiation
  #and thus is adjusted for latitude/solarzen before comparison with Rsurface
  #I believe sirad:extrat provides latitude-corrected (ie horizontal) extraterrestrial radiation
  IdoverI <- IdoverI.calc(kt,solarzen)
  IdoverI[IdoverI < 0] <- 0 ; IdoverI[IdoverI > 1] <- 1
  SpSd <- Rsurface * cbind(1-IdoverI,IdoverI)
  return(SpSd)
}

solarzen.calc <- function(coords,datetimePOSIXct){
  #calculate solar zenith based on lat, lon & time of day
  # Args: 
  #  coords: coordinates
  #  datetimePOSIXct: a POSIXct object giving date and time
  # Returns:
  #  solar zenith angle

  require(maptools)
  solarelv <- solarpos(crds=coords,dateTime=datetimePOSIXct)
  solarelv <- solarelv[,2]
  #convert to radians
  solarelv <- pi*solarelv/180
  #convert to solar zenith (0 when sun is overhead) in radians
  solarzen <- pi/2 - solarelv
  #set below-horizon zeniths to pi/2
  solarzen[solarzen > pi/2] <- pi/2
  cat("in zolar zenith values, pi/2 (90degs) represents horizon/sub-horizon\n
      while 0 represents directly over-head\n")
  return(solarzen)  
}

R_extra_terr.calc<-function(thisdate,lat.in.deg){
  #calculate hourly extraterrestrial irradiance in W/m2 using sirad package
  # Args: 
  #  thisdate: string object giving date, e.g. "2011-12-31"
  #  lat.in.deg: lattitude in degrees
  # Returns:
  #  hourly extraterrestrial radiation
  # Example:
  #  R_extra_terr.calc(thisdate=c("2012-01-19","2012-01-19"),lat.in.deg=c(-69))
  
  require(sirad)
  JulianDay <- sirad::dayOfYear(thisdate)#dayOfYear("2011-01-01")
  lat.in.rad <- lat.in.deg*pi/180
  R_extra_terr <- extrat(i=JulianDay,lat.in.rad)#[[2]]
  #extrat[[1]] uses MJ/(day m^2) units
  #extrat[[2]] uses MJ/(hr m^2) units
  #convert output to W/m2
  conv.fac.daily <- (1000000/ 86400)
  conv.fac.hourly <- (1000000/ 86400)*24
  R_extra_terr[[1]] <- conv.fac.daily *  R_extra_terr[[1]]
  R_extra_terr[[2]] <- conv.fac.hourly *  R_extra_terr[[2]]
  #set negative hourly extrat irrad to 0
  #R_exta_terr <- R_exta_terr[[2]]
  R_extra_terr[[2]][R_extra_terr[[2]]<0] <- 0
  R_extra_terr <- R_extra_terr[[2]]
  R_extra_terr <- matrix(R_extra_terr,ncol=24)
  #cat("calculated 24hrs of extraterrestrial irradiance for\n", nrow(R_extra_terr)," latitude-date combinations\n")
  #output is hourly starting at solar midning
  return(R_extra_terr)
}

shift.vec <- function(vec,n,wrap=TRUE,pad=FALSE){
  # Shift a vector over by n spots
  #Args: 
  #  vec: vector to be shifted
  #  n: number of spots to shift the vector over
  #  Wrap adds the entry at the beginning to the end
  #  pad does nothing unless wrap is false, in which case it specifies whether to pad with NAs
  # Returns:
  #  the vector vec, shifted over n sots
  # Source:
  #  http://stackoverflow.com/questions/6828937/what-to-do-with-imperfect-but-useful-functions
  
  if(length(vec)<abs(n)) { 
    #stop("Length of vector must be greater than the magnitude of n \n") 
  }
  if(n==0) { 
    return(vec) 
  } else if(length(vec)==n) { 
    # return empty
    length(vec) <- 0
    return(vec)
  } else if(n>0) {
    returnvec <- vec[seq(n+1,length(vec) )]
    if(wrap) {
      returnvec <- c(returnvec,vec[seq(n)])
    } else if(pad) {
      returnvec <- c(returnvec,rep(NA,n))
    }
  } else if(n<0) {
    returnvec <- vec[seq(1,length(vec)-abs(n))]
    if(wrap) {
      returnvec <- c( vec[seq(length(vec)-abs(n)+1,length(vec))], returnvec )
    } else if(pad) {
      returnvec <- c( rep(NA,abs(n)), returnvec )
    }
    
  }
  return(returnvec)
}

shift_to_UTC <- function(hourly_tser,original.zone){
  #shift an hourly vector starting at midnight to start at midnight UTC
  # Args: 
  #  hourly_tser: an hourly vector starting at midnight local time
  #  original.zone: the local time zone
  # Returns:
  #  an hourly vector, starting at midnight UTC
  # Example:
  #  shift_to_UTC(hourly_tser=1:24,original.zone="America/Mexico_City")
  
  require(timeDate)
  tt1<-timeDate("2010-01-01 00:00:00",zone=original.zone)
  tt2<-timeDate("2010-01-01 00:00:00",zone="UTC")
  shift.by <- as.numeric(tt2-tt1)
  cat("the 24 hour series will be shifted earlier hrs by ",shift.by," hours to match UTC\n")
  if (is.matrix(hourly_tser)){
    shifted_tser<-t(apply(hourly_tser,1,shift.vec,n=shift.by))
    }else{shifted_tser <- shift.vec(hourly_tser,n=shift.by)}
  return(shifted_tser)
}

R_extra_for_site.vec <- function(thisdate,lat.in.deg,original.zone){
  #given a date, latitude, and time zone, calculate extraterrestrial radiation
  #for that latitude
  #Args:
  #  thisdate: date
  #  lat.in.deg: latitude in degrees
  #  original.zone: timezone
  # Returns:
  #  Extraterrestrial radiation for the given latitude and day, in 6 hour steps
  
  if(length(thisdate)!=length(lat.in.deg)){
    cat("Please provide thisdate and lat.in.deg of equal length\n");browser()}
  #get the 24 R_extra_terr values for this lat & date
  R_extra_terr_24 <- R_extra_terr.calc(thisdate=thisdate,lat.in.deg=lat.in.deg)
  #determine by how many hours the series need to be shifted
  R_extra_terr_24UTC <- shift_to_UTC(R_extra_terr_24,original.zone=original.zone)
  #convert the R_extra_terr values from 24 hour to epoch
  epochmn <- function(x){tapply(x,rep(1:4,each=6),mean)}
  R_extra_terr_epochUTC <- t(apply(R_extra_terr_24UTC,1,epochmn))
  cat("all 24 hour series of extra-terrestrial R converted to 4 6-hour UTC epochs\n")
  rm(R_extra_terr_24,R_extra_terr_24UTC)
  return(R_extra_terr_epochUTC)
}
