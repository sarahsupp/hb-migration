# Add Physiologically informed variables, based on Don's and Pieter's equations.

# define pathnames
agan.dir <- "C:/Users/sarah/Dropbox/ebird_annotated_raw/"
function.dir <- "C:/Users/sarah/Documents/github/hb-migration/"

# source function script
source(paste0(function.dir, "hb_RS_functions.R"))
source(paste0(function.dir, "Standard_Operative_Temperature_from_meteorology.R"))
source(paste0(function.dir, "Diffuse_fraction_of_solar_radiation.r"))

# species codes
spcodes <- c("rthu", "bchu", "bthu", "cahu", "ruhu")

for (spp in unique(spcodes)){
  
  #read in annotaed data
  ann <- read.csv(paste0(agan.dir, spp, "/", spp, "_lag0_allYears.csv"), as.is=T)
  
  
  #calculate solar zenith (sun angle) and horizontal extraterrestrial radiation
  solarzen <- solarzen.calc(coords, datetimePOSIXct)
  R_extraterr <- R_extra_terr.calc(thisdate, lat.in.deg)
  
  #calculate standard operative temperature from surface temperature (Te)
  Tes <- Tes.calc.compl.incl.rad(Ta, u, Li, Rsurface, R_extra_terr,solarzen){
    #   Ta: ambient T (Kelvin)
    #   u: wind speed (m/s)
    #   Li: incoming longwave radition (W m^(-2))
    #   RSurface: incoming shortwave radiation at the surface
    #   R_extra_terr: horizontal extraterrestrial radiation 
    #   solarzen:  solar zenith angle in radians
  
    #Convert Tes (K) to Celsius
    ann$Tes <- ann$Tes - 273.15
    
  #calculate physiological demand based on Tes
    #   Assume S. rufus BMR = 3.3 mL O2 g-1h-1 (Lasiewski 1963) #TODO: get BMR for other species, or assume all the same?
    if (Tes <= 35){
      Pd <- 17.8 - (0.396 * Tes) - 3.3
    }
    else if (Tes > 35){
      Pd <- (0.214 * Tes) - 7.49
    }
  
  #convert Pd into energetic cost (Joules)
    # Thermoregulation (Joules) = Thermoregulation (mL O2 g-1h-1) * 20.1
    PdJ <- Pd * 20.1
  
  #calculate transport costs (cost of flight) at different Tes
    #   Hovering metabolic rate (HMR) in S. rufus is: HMR = 58.9 mL O2 g-1h-1 OR 17.8 x BMR 
    #   Forward flight costs (FLMR) in S. rufus is: FLMRRufous = 0.49 x HMR
    #   transportTime is in hours
    Tc <- transportTime * 0.49 * 17.8 * BMR
  
  
  
  