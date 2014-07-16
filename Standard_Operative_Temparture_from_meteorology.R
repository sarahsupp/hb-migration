# Diffuse_fraction_of_solar_radiation.R
#
# Pieter Beck (psabeck@gmail.com)
# 09-12-2013
# 
#Calculate Standard Operative T for hummingbirds using meteo input data 
#based on Don Powers's email "Estimating Operative Temperature"
#
#
# Changelog:
#   VERS  | DATE    	 | CHANGES				| BY
#   ------|------------|----------------|---- 
#		1.0.0	|	09-12-2013 | Wrote script		| PB


ra.calc <- function(u){
  #caclulate aerodynamic resistance to convective heat transfer
  #using an equation below was estimated for a small sparrow
  # Args:
  #  u: wind speed in m/s
  #
  # Returns:
  #  aerodynamic resistance to convective heat transfer (s m-1)
  # 
  ra <- 37.76*u^(-0.945)
  return(ra)
}  

Rabs.calc <- function(Sp, ApOVERA = 0.25, Sd, aS = 0.89, aL = 0.9, Li){
  #Calculate long + shortwave radiation absorbed
  #
  # Args:
  #  Sp: direct shortwave radiation (W m^(-2))
  #  ApOVERA: ratio of projected shadow area to total surface area. 
  #           This is basically used to measure how much of the surface area is 
  #           directly absorbing shortwave radiation. 0.25 is a conservative starting estimate
  #   Sd: diffuse shortwave radiation. This is usually measured with a shadow radiometer
  #   aS: mean shortwave absorptivity. 0.89 was estimated for 2 bird species
  #   aL: mean longwave absorptivity. Don will find value 0.9 is place holder !!!!!!
  #   Li: incoming longwave radiation
  
  Rabs.shortwave <- ((Sp*ApOVERA) + Sd) * aS
  Rabs.longwave <-  Li * aL
  Rabs <- Rabs.shortwave + Rabs.longwave
  return(Rabs)
}

Te.calc <- function(Ta, Rabs, sigma = 5.67e-8, epsilon = 0.95, rhoCp = 1200, ra){
  #Calculate operative temperature
  #
  # Args:
  #  Ta: ambient T
  #  Rabs: long-wave+shortwave radiation absorbed (W m-2)
  #  sigma: 5.67e-8, i.e. Stefan Boltman constant W m^(-2) K^(-4)
  #  epsilon: 3missivity of the bird's surface
  #  rhoCp: product of the density and specific heat capacity of air (J m-3 K^(-1))
  #  ra: aerodynamic resistance to convective heat transfer (s m-1)
  #
  # Returns:
  #  Operative temperature
  #
  if (Ta < 263){cat("STOP !!! Ta needs to be provided in K !!!\n");
                cat("provided Ta range is: ",range(Ta,na.rm=T),"\n");browser()}
  Te <- Ta + (Rabs - sigma*epsilon*(Ta^4))/(rhoCp / 4*sigma*(Ta^3) + ra)
  return(Te)
}

Tes.calc <- function(Te, Tb, u){
  #Calculate standard operative temperature
  #
  # Args
  #   Te: Operative temperature
  #   Tb: 42+273, #42 C, or 42 + 273 K
  #   u: wind speed in m/s)
  #
  # Returns
  #   Standard Operative Temperature
  Tes <- Tb - (1 + (0.26 * sqrt(u))*(Tb-Te))
  return(Tes)
}

Tes.calc.compl <- function(Ta, u, Sp, Sd, Li){
  #A comprensive function to estimate standard operative temperature from basic meteo. variables,
  #but processed radiation data
  #
  # Args:
  #   Ta: ambient T (Kelvin)
  #   u: wind speed (m/s)
  #   Sp: direct shortwave radiation (W m^(-2)) 
  #   Sd: diffuse shortwave radiation (W m^(-2)) 
  #   Li: incoming longwave radition (W m^(-2))
  #
  # Returns:
  #   Standard Operative Temperature
  #
  #calc aerodynamic resistance to convective heat transfer
  ra <- ra.calc(u=u)
  #calc longwave+shortwave radiation absorbed
  Rabs <- Rabs.calc(Sp=Sp, Sd=Sd, Li=Li) 
  Te <- Te.calc(Ta=Ta,ra=ra,Rabs=Rabs)
  Tes <- Tes.calc(Te,u=u)
  return(Tes)
}


Tes.calc.compl.incl.rad <- function(Ta,
                           u,
                           Li,
                           Rsurface,
                           R_extra_terr,
                           solarzen                                  
                  ){
#A comprensive function to estimate standard operative temperature from basic meteo. variables,
#and raw radiation data (drawing on functions in Diffuse_fraction_of_solar_radiation.r)
#
# Args:
#   Ta: ambient T (Kelvin)
#   u: wind speed (m/s)
#   Li: incoming longwave radition (W m^(-2))
#   RSurface: incoming shortwave radiation at the surface
#   R_extra_terr: horizontal extraterrestrial radiation 
#   solarzen:  solar zenith angle in radians
#
# Returns:
#   Standard operative temperature
#
  #load the function SpSd.calc
  source("C:\\share\\pbeck\\Hummer_NASA\Code_copy\\Diffuse_fraction_of_solar_radiation.R")
  #
  SpSd <- SpSd.calc(Rsurface=Rsurface,R_extra_terr=R_extra_terr,solarzen=solarzen)
  #win.graph();par(mfrow=c(2,2))
  #hist(Rsurface);hist(R_extra_terr);hist(SpSd[,1]);hist(SpSd[,2])
  Tes <- Tes.calc.compl(Ta=Ta,u=u,Sp=SpSd[,1],Sd=SpSd[,2],Li=Li)
  return(Tes)
}


