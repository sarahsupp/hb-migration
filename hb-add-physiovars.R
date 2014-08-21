# Add Physiologically informed variables, based on Don's and Pieter's equations.
library(sp)
library(ggplot2)

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
  
  #calculate daylength from the dates and locations
  daylength <- apply(ann, 1, function(x){
    coords <- matrix(c(as.numeric(x["location.long"]), as.numeric(x["location.lat"])), nrow=1)
    #coords <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84"))
    sunrise <- sunriset(coords, as.POSIXct(x["timestamp"]), direction="sunrise", POSIXct.out=TRUE)
    sunset <- sunriset(coords, as.POSIXct(x["timestamp"]), direction="sunset", POSIXct.out=TRUE)
    sunset$time - sunrise$time
  })
  
  # TODO: Test to make sure windspeed and winddir work
  #calculate windspeed (m/s) from the east-west(u10m) and north-south(v10m) wind components
  windspeed <- apply(ann, 1, function(x){
    speed = sqrt( x["u10m"]^2 + x["v10m"]^2 ) 
    print(speed)
  })
  
  #calculate wind direction from the east-west(u10m) and north-south(v10m) wind components 
  winddir <- apply(ann, 1, function(x){
    dir = atan( x["v10m"]/x["u10m"] ) 
    print(dir)
  })
  
  #calculate solar zenith (sun angle) and horizontal extraterrestrial radiation
    #coords need to be a spatial points or matrix object
    # apply across all rows of the annotated data
  solarzen <- apply(ann, 1, function(x){
    coords <- matrix(c(as.numeric(x["location.long"]), as.numeric(x["location.lat"])), nrow=1)
    solarzen.calc(coords, as.POSIXct(x["timestamp"]))
  })
  
    R_extra_terr <- apply(ann, 1, function(x){
      date <- as.Date(x[1], format='%Y-%m-%d %H:%M:%S.000')
      mean(R_extra_terr.calc(date, as.numeric(x["location.lat"]))) #TODO: output is Ret for 24 hours. Use mean? max? mode? noon?
    })
  
  #Append new variables to dataframe
  ann$solarzen <- solarzen
  ann$R_extra_terr <- R_extra_terr
  ann$daylength <- daylength
  ann$windspeed <- windspeed
  ann$winddir <- winddir
  
  #calculate standard operative temperature from surface temperature (Te)
  #Use temperature and wind at 10m for consistency. Hb are typically observed by bird watchers at relatively low heights.
  #TODO: wind is in m/s (should we be using absolute value for wind?)
  #TODO: u10m should be combined with v10m for actual wind speed (new windspeed variable)? @TinaCormier
  #TODO: Commented out flag for Ta < 263 (but should never get a presence point for such data) - check this is OK - maybe do an additional filter on presence data where extremely cold
    Tes <- apply (ann, 1, function(x){
      Tes.calc.compl.incl.rad(as.numeric(x["t10m"]), abs(as.numeric(x["windspeed"])), as.numeric(x["lwrf"]), 
                              as.numeric(x["swrf"]), as.numeric(x["R_extra_terr"]), as.numeric(x["solarzen"]))
    })
  
  #Append Tes in K and C to dataframe
  ann$TesK <- Tes  
  ann$TesC <- ann$Tes - 273.15 #convert Tes (K) to Celsius
    
  #calculate physiological demand (Joules) based on Tes (c)
    #   Assume S. rufus BMR = 3.3 mL O2 g-1h-1 (Lasiewski 1963) 
    # Thermoregulation (Joules) = Thermoregulation (mL O2 g-1h-1) * 20.1
    #TODO: get BMR for other species, or assume all similar to S. rufus and C. costae (~3 g)?
  PdJ <- apply (ann, 1, function(x){
    Tes <- as.numeric(x["TesC"])
    if (Tes <= 35){ thermoreg <- 17.8 - (0.396 * Tes) - 3.3 }
    else if (Tes > 35){ thermoreg <- (0.214 * Tes) - 7.49 }
    thermoreg * 20.1
  })
  
  # append physiological demand in Joules to the dataframe
  ann$Pdj <- PdJ
  
  #write new dataframe to file
  write.csv(ann, file=paste0(file.dir, spp, "/", spp, "_lag0_allYears_fil_phys.csv"), sep=",", row.names=FALSE)
  
  
#----------------------------- Plot the physiological data
pres = ann[ann$presence==1,]
abs = ann[ann$presence==0,]

ggplot(pres, aes(t10m, TesK)) + geom_point(alpha=0.05) + theme_classic()
ggplot(pres, aes(abs(windspeed), TesK)) + geom_point(alpha=0.05) + theme_classic()

ggplot(abs, aes(t10m - 273.15, Pdj)) + geom_point(col="red", alpha=0.05) + theme_classic() + 
  ylab ("Physiological demand (Joules)") + xlab("Temperature (C) at 10 m") + 
  geom_point(data=pres, aes(Temp_sfc - 273.15, Pdj), alpha=0.05)

ggplot(abs, aes(abs(windspeed), Pdj)) + geom_point(col="red", alpha=0.05) + theme_classic() + 
  ylab ("Physiological demand (Joules)") + xlab("wind speed at 10 m (m/s)") + 
  geom_point(data=pres, aes(abs(windspeed), Pdj), alpha=0.05)
  
}

