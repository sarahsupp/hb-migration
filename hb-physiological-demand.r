# Calculates Operative Temperature (Te), Thermoregulation Cost (Ct), and Transport/Flight Cost (Cf)
# equations from published materials (Bakken 1976, Bakken 1980, Walsberg 1992) and Don Powers (pers. comm.)
# Model functions written orginally by P.A. Beck 2013 (unpublished code) and Don Powers
# (c) 2014 Sarah R. Supp

setwd("C:/Users/sarah/Documents/GitHub/hb-migration/")
source("Diffuse_fraction_of_solar_radiation.r")
source("Standard_Operative_Temparture_from_meteorology.r")

# Calculate standard operative temperature (need ambient temp (K),  wind speed (u), need swrf [direct (Sp) and diffuse (Sd)] 
# and lwrf (Li), and body temperature (Tb) in K, which is assumed to be ~42 C or 315 K
Tes = Tes.calc.comp(Ta, u, Sp, Sd, Li, 315)

# Based on Tes:
#   Calculate Thermoregulation cost (Ct)


#   Calculate Transport (Flight) cost (Cf)