# Add Physiologically informed variables, based on Don's and Pieter's equations.
library(sp)
library(ggplot2)
library(maptools)
library(plyr)

# define pathnames
file.dir <- "C:/Users/sarah/Dropbox/hb_migration_data/ebird_annotated_fil/"
function.dir <- "C:/Users/sarah/Documents/github/hb-migration/"

# source function script
source(paste0(function.dir, "hb_RS_functions.R"))
source(paste0(function.dir, "Standard_Operative_Temparture_from_meteorology.R"))
source(paste0(function.dir, "Diffuse_fraction_of_solar_radiation.r"))

# species codes
spcodes <- c("rthu", "bchu", "bthu", "cahu", "ruhu")

for (spp in unique(spcodes)){
  
  #read in annotaed data
  ann <- read.csv(paste0(file.dir, spp, "/", spp, "_lag0_allYears_fil.csv"), as.is=T)
  tstamp <- as.POSIXct(ann$timestamp, format='%Y-%m-%d %H:%M:%S') #format time variables TODO: extract time variables
  ann$time <- factor(format(tstamp, "%H:%M:%S"), ordered=T) #TODO: Check that these times are OK and make sense - filter data with "bad"/unlikely times
  
  #calculate daylength from the dates and locations
  daylength <- apply(ann, 1, function(x){
    coords <- matrix(c(as.numeric(x["location.long"]), as.numeric(x["location.lat"])), nrow=1)
    #coords <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84"))
    sunrise <- sunriset(coords, as.POSIXct(x["timestamp"]), direction="sunrise", POSIXct.out=TRUE)
    sunset <- sunriset(coords, as.POSIXct(x["timestamp"]), direction="sunset", POSIXct.out=TRUE)
    sunset$time - sunrise$time
  })
  
  #calculate windspeed (m/s) from the east-west(u10m) and north-south(v10m) wind components
  windspeed <- apply(ann, 1, function(x){
    speed = sqrt( as.numeric(x["u10m"])^2 + as.numeric(x["v10m"])^2 ) 
  })
  
  #calculate wind direction from the east-west(U) and north-south(V) wind components 
  # where  with 0° or 360° indicating a wind blowing to the north, 90° indicating a wind blowing to the east, 
  # 180° indicating a wind blowing to the south and 270° indicating a wind blowing to the west.
  windir <- function(U, V){
    dir = (270 - atan2(V,U) * 180/pi)%%360
    return(dir)
  }
  
  winddir <- apply(ann, 1, function(x){
    dir = windir(as.numeric(x["u10m"]), as.numeric(x["v10m"])) 
  })
  
  #calculate solar zenith (sun angle) and horizontal extraterrestrial radiation
    #coords need to be a spatial points or matrix object
    # apply across all rows of the annotated data
  solarzen <- apply(ann, 1, function(x){
    coords <- matrix(c(as.numeric(x["location.long"]), as.numeric(x["location.lat"])), nrow=1)
    solarzen.calc(coords, as.POSIXct(x["timestamp"]))
  })
  
  R_extra_terr <- apply(ann, 1, function(x){
    date <- as.Date(x[1], format="%Y-%m-%d %H:%M:%S.000")
    Ret_hrs <- R_extra_terr.calc(date, as.numeric(x["location.lat"])) 
    #print(Ret_hrs)
    tt <- strptime(x[1], format="%Y-%m-%d %H:%M:%S.000")
    roundtime <- as.numeric(format(round(tt, units="hours"), format="%H"))
    Ret <- as.numeric(Ret_hrs[,roundtime])
    return(Ret)
  })
  
  #Append new variables to dataframe
  ann$daylength <- daylength
  ann$windspeed <- windspeed
  ann$winddir <- winddir
  ann$solarzen <- solarzen
  ann$R_extra_terr <- as.vector(R_extra_terr)
  
  #calculate standard operative temperature from surface temperature (Te)
  #Use temperature and wind at 10m for consistency. Hb are typically observed by bird watchers at relatively low heights.
  #TODO: Commented out flag for Ta < 263 (but should never get a presence point for such data) - check this is OK - maybe do an additional filter on presence data where extremely cold 
  Tes <- apply (ann, 1, function(x){
      Tes.calc.compl.incl.rad(as.numeric(x["t10m"]), abs(as.numeric(x["windspeed"])), as.numeric(x["lwrf"]), 
                              as.numeric(x["swrf"]), as.numeric(x["R_extra_terr"]), as.numeric(x["solarzen"]))
    })
  
  #Append Tes in K and C to dataframe
  ann$TesK <- Tes  
  ann$TesC <- ann$TesK - 273.15 #convert Tes (K) to Celsius
    
  #calculate physiological demand (Joules) based on Tes (c)
    #   Assume S. rufus BMR = 3.3 mL O2 g-1h-1 (Lasiewski 1963) 
    # Thermoregulation (Joules) = Thermoregulation (mL O2 g-1h-1) * 20.1
    #TODO: get BMR for other species, or assume all similar to S. rufus and C. costae (~3 g)?
  PdJ <- apply (ann, 1, function(x){
    Tes <- as.numeric(x["TesC"])
    if (Tes <= 35){ thermoreg <- 17.8 - (0.396 * Tes) - 3.3 }
    else if (Tes > 35){ thermoreg <- (0.214 * Tes) - 7.49 }
    Joules = thermoreg * 20.1
    return(Joules)
  })
  
  # append physiological demand in Joules to the dataframe
  ann$Pdj <- PdJ
  
  #write new dataframe to file
  write.table(ann, file=paste(file.dir, spp, "/", spp, "_lag0_allYears_fil_phys.csv", sep=""), sep=",", row.names=FALSE)
  
  
#----------------------------- Plot the physiological data
pres = ann[ann$presence==1,]
abs = ann[ann$presence==0,]

ggplot(ann, aes(t10m, TesK)) + geom_point(alpha=0.05) + theme_classic()
ggplot(ann, aes(windspeed, TesK)) + geom_point(alpha=0.05) + theme_classic()

ggplot(abs, aes(TesC, Pdj)) + geom_point(alpha=0.05, col="red") + theme_classic() +
  ylab ("Physiological demand (Joules)") + xlab("Operative Temperature (C) at 10 m") + 
  geom_point(data=pres, aes(TesC, Pdj), alpha=0.05)

ggplot(abs, aes(t10m - 273.15, Pdj)) + geom_point(col="red", alpha=0.05) + theme_classic() + 
  ylab ("Physiological demand (Joules)") + xlab("Ambient Temperature (C) at 10 m") + 
  geom_point(data=pres, aes(Temp_sfc - 273.15, Pdj), alpha=0.05)

ggplot(abs, aes(windspeed, Pdj)) + geom_point(col="red", alpha=0.05) + theme_classic() + 
  ylab ("Physiological demand (Joules)") + xlab("wind speed at 10 m (m/s)") + 
  geom_point(data=pres, aes(abs(windspeed), Pdj), alpha=0.05)
  
}

