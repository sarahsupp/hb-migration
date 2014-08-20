# Add Physiologically informed variables, based on Don's and Pieter's equations.

# define pathnames
agan.dir <- "C:/Users/sarah/Dropbox/ebird_annotated_raw/"
function.dir <- "C:/Users/sarah/Documents/github/hb-migration/"

# source function script
source(paste0(function.dir, "hb_RS_functions.R"))
source(paste0(function.dir, "Standard_Operative_Temparture_from_meteorology.R"))
source(paste0(function.dir, "Diffuse_fraction_of_solar_radiation.r"))

# species codes
spcodes <- c("rthu", "bchu", "bthu", "cahu", "ruhu")

for (spp in unique(spcodes)){
  
  #read in annotaed data
  ann <- read.csv(paste0(agan.dir, spp, "/", spp, "_lag0_allYears.csv"), as.is=T)
  
  #use only complete data
  ann <- ann[complete.cases(ann[,-4]),]
  
  #calculate solar zenith (sun angle) and horizontal extraterrestrial radiation
    #coords need to be a spatial points or matrix object
    # apply across all rows of the annotated data
  solarzen <- apply(ann, 1, function(x){
    coords <- matrix(c(as.numeric(x["location.long"]), as.numeric(x["location.lat"])), nrow=1)
    solarzen.calc(coords, as.POSIXct(x["timestamp"]))
  })
  
    R_extra_terr <- apply(ann, 1, function(x){
      date <- as.Date(x[1], format='%Y-%m-%d %H:%M:%S.000')
      mean(R_extra_terr.calc(date, as.numeric(x["location.lat"]))) #TODO: output is Ret for 24 hours. Use mean? max? mode?
    })
    
  #Append solar zen and R_extra_terr to dataframe
  ann$solarzen <- solarzen
  ann$R_extra_terr <- R_extra_terr
  
  
  #calculate standard operative temperature from surface temperature (Te)
  #TODO: wind is in m/s (should we be using absolute value for wind?)
  #TODO: Commented out flag for Ta < 263 (but should never get a presence point for such data) - check this is OK
    Tes <- apply (ann, 1, function(x){
      Tes.calc.compl.incl.rad(as.numeric(x["Temp_sfc"]), abs(as.numeric(x["u10m"])), as.numeric(x["lwrf"]), 
                              as.numeric(x["swrf"]), as.numeric(x["R_extra_terr"]), as.numeric(x["solarzen"]))
    })
  
    #Append Tes to data frame and convert Tes (K) to Celsius
  ann$TesK <- Tes  
  ann$TesC <- ann$Tes - 273.15
    
  
  #calculate physiological demand (Joules) based on Tes (c)
    #   Assume S. rufus BMR = 3.3 mL O2 g-1h-1 (Lasiewski 1963) 
    # Thermoregulation (Joules) = Thermoregulation (mL O2 g-1h-1) * 20.1
    #TODO: get BMR for other species, or assume all similar to S. rufus (~3 g)?
  PdJ <- apply (ann, 1, function(x){
    Tes <- as.numeric(x["TesC"])
    if (Tes <= 35){ thermoreg <- 17.8 - (0.396 * Tes) - 3.3 }
    else if (Tes > 35){ thermoreg <- (0.214 * Tes) - 7.49 }
    thermoreg * 20.1
  })
  
  # append physiological demand in Joules to the dataframe
  ann$Pdj <- PdJ
  
  
#----------------------------- Plot the physiological data
pres = ann[ann$presence==1,]
abs = ann[ann$presence==0,]

ggplot(pres, aes(Temp_sfc, TesK)) + geom_point(alpha=0.05) + theme_classic()
ggplot(pres, aes(abs(u10m), TesK)) + geom_point(alpha=0.05) + theme_classic()

ggplot(abs, aes(Temp_sfc - 273.15, Pdj)) + geom_point(col="red", alpha=0.05) + theme_classic() + 
  ylab ("Physiological demand (Joules)") + xlab("Surface Temperature (C)") + 
  geom_point(data=pres, aes(Temp_sfc - 273.15, Pdj), alpha=0.05)

ggplot(abs, aes(abs(u10m), Pdj)) + geom_point(col="red", alpha=0.05) + theme_classic() + 
  ylab ("Physiological demand (Joules)") + xlab("wind speed at 10 m (m/s)") + 
  geom_point(data=pres, aes(abs(u10m), Pdj), alpha=0.05)
  
}
  