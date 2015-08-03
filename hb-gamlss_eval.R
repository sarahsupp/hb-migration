# this code is to run gamlss models for the eBird hummingbird data with environmental correlates
# (c) 2015 by Sarah R. Supp, Laura Graham, and Frank La Sorte

# import libraries
library(nlme)
library(lme4)
library(lubridate)
library(arm)
library(mgcv)
library(gamm4)
library(gamlss)

# define pathnames (FAL will need to change pathnames)
function.dir <- "/home/sarah/Documents/GitHub/hb-migration/hb_RS_functions.R"
agan.dir <- "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_annotated_raw/combined/"
migtime.dir <- "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"
fig.dir <- "/home/sarah/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/"
stat.dir <- "/home/sarah/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/stat-tests/"

source(function.dir)

#assign based on the time frame (number of days) used to compute the alpha hulls
alpha_window = 5

#list of species codes
species <- c("bchu", "bthu", "cahu", "ruhu", "rthu")

###########################################################################################
#             FUNCTION TO EVALUATE GAMLSS MODELS
###########################################################################################

eval_gammlss_models = function(dat=dat, sp=sp, lag=FALSE, means=FALSE) {
  # evaluates gamlss models for each dataset based on the hypothesis testing framework
  # models are run on centered and scaled (such that mean=0, sd=1) predictor variables
  # models input year and window (alpha hull) as random effects
  # if lag==FALSE, models will be run with the set of hypotheses for within-region signal of presence
  # if lag==TRUE, models will be run with the set of hypotheses for migration trajectory with timelag
  # if means==FALSE, the data represents all the records from the alpha hulls
  # if means==TRUE, the data has been averaged by window (each window represents a unique alpha hull)
  # outputs model results, AIC comparison, R2, and plots the data

  if(!means){
    cat("data represents all records from the alpha hulls \n")
    
    #-----For within-region signal of presence (presence vs. absence points)
    if(!lag){
      cat("models will compare within region presence vs. absence \n")
      
      # H0a. Remotely sensed variables will not be predictive of hummingbird presence, 
      #     because they do not adequately capture the processes (resources, demand) 
      #     that determine where birds are found along their migration route. 
      MNull = gamlss(pres ~ 1 + (1|yr) + (1|window), data=dat)
      
      # H1. In the spring, birds should be more strongly tied to EVI, assuming that greenness 
      #     is a better proxy for food resources in the spring than in the fall.
      # H2. In the fall birds should be constrained to temperatures and Elevation*EVI interaction (resources), 
      #     particularly in the west where there is more extreme topography and habitat heterogeneity.
      M1 = gamlss(pres ~ scale(EVI) * scale(t10m) * scale(STRM_elev) + (1|yr) + (1|window),data=dat)
      M2 = gamlss(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|yr) + (1|window), data=dat)
      M3 = gamlss(pres ~ scale(EVI) + scale(t10m) + (1|yr) + (1|window), data=dat)
      M4 = gamlss(pres ~ scale(EVI) + (1|yr) + (1|window), data=dat)
      M5 = gamlss(pres ~ scale(t10m) + (1|yr) + (1|window), data=dat)
      
      GAIC(MNull, M1, M2, M3, M4, M5) #returns matrix with AIC and df 
      
      # return Cox Snell R2, R-squared = 1-(L(0)/L(fitted))^(2/n)
      # from: http://rpackages.ianhowson.com/cran/gamlss/man/Rsq.html
      MNull.rsq= Rsq(MNull, type="Cox Snell")
      M1.rsq = Rsq(M1, type="Cox Snell")
      M2.rsq = Rsq(M2, type="Cox Snell")
      M3.rsq = Rsq(M3, type="Cox Snell")
      M4.rsq = Rsq(M4, type="Cox Snell")
      M5.rsq = Rsq(M5, type="Cox Snell")
    }
    
    #-----For migratory trajectory (presence vs +/- 15 day points)
    if(lag){
      cat("models will compare presence with migratory trajectory lags (+/- 15 days \n")
      
      # H0b. Remotely sensed variables will not be predictive of hummingbird presence, 
      #     because they do not adequately capture the processes (resources, demand) 
      #     that determine where birds are found along their migration route.
      MNull = gammlss(pres ~ 1 + (1|yr) + (1|window), data=pmin.spring)
      
      # H3. In the spring, birds will be more strongly tied to EVI “greening up” 
      #     as they move northward to breeding grounds than in the fall.
      # H4. In the fall, birds will be more limited by temperature, as late summer/early fall 
      #     temperatures tend to be higher (e.g. closer to lethal limit) than in the spring.
      # H5. In the fall, birds should be associated with higher elevations than in the spring, 
      #     because flower phenology (blooming time) is later at higher altitudes 
      #     (so it is associated with more resources) and cooler temperatures are found at higher elevations.
      M1 = gammlss(pres ~ scale(EVI) * scale(t10m) * scale(SRTM_elev) + (1|yr) + (1|window), data=dat)
      M2 = gammlss(pres ~ scale(temp) * scale(elev) + scale(EVI) * scale(SRTM_elev) + (1|yr) + (1|window), data=dat)
      M3 = gammlss(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|yr) + (1|window), data=dat)
      M4 = gammlss(pres ~ scale(EVI) + scale(t10m) + (1|yr) + (1|window), data=dat)
      M5 = gammlss(pres ~ scale(EVI) + (1|yr) + (1|window), data=dat)
      M6 = gammlss(pres ~ scale(t10m) + (1|yr) + (1|window), data=dat)
      
      
      GAIC(MNull, M1, M2, M3, M4, M5, M6) #returns matrix with AIC and df 
      
      # return Cox Snell R2, R-squared = 1-(L(0)/L(fitted))^(2/n);  http://rpackages.ianhowson.com/cran/gamlss/man/Rsq.html
      MNull.rsq= Rsq(MNull, type="Cox Snell")
      M1.rsq = Rsq(M1, type="Cox Snell")
      M2.rsq = Rsq(M2, type="Cox Snell")
      M3.rsq = Rsq(M3, type="Cox Snell")
      M4.rsq = Rsq(M4, type="Cox Snell")
      M5.rsq = Rsq(M5, type="Cox Snell")
      M6.rsq = Rsq(M6, type="Cox Snell")
    }
  }
  if(means){
    cat("data has been averaged within unique alpha hulls (window) \n")
    cat("effort will be accounted for as number of checklists where species was present \n")
    
    #-----For within-region signal of presence (presence vs. absence points)
    if(!lag){
      cat("models will compare within region presence vs. absence \n")
      
      #models use the same basic hypotheses as above, for no lag, but include effort
      MNull = gamlss(pres ~ 1 + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M1 = gamlss(pres ~ scale(mean.EVI) * scale(mean.t10m) * scale(mean.elev) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M2 = gamlss(pres ~ scale(mean.EVI) + scale(mean.t10m) + scale(mean.elev) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M3 = gamlss(pres ~ scale(mean.EVI) + scale(mean.t10m) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M4 = gamlss(pres ~ scale(mean.EVI) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M5 = gamlss(pres ~ scale(mean.t10m) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      
      GAIC(MNull, M1, M2, M3, M4, M5, M6) #returns matrix with AIC and df 
      
      # return Cox Snell R2, R-squared = 1-(L(0)/L(fitted))^(2/n);  http://rpackages.ianhowson.com/cran/gamlss/man/Rsq.html
      MNull.rsq= Rsq(MNull, type="Cox Snell")
      M1.rsq = Rsq(M1, type="Cox Snell")
      M2.rsq = Rsq(M2, type="Cox Snell")
      M3.rsq = Rsq(M3, type="Cox Snell")
      M4.rsq = Rsq(M4, type="Cox Snell")
      M5.rsq = Rsq(M5, type="Cox Snell")
    }
    
    #-----For migratory trajectory (presence vs +/- 15 day points)
    if(lag){
      cat("models will compare presence with migratory trajectory lags (+/- 15 days \n")
      
      #models use the same basic hypotheses as above, for +/- 15 day lag, but include effort
      MNull = gammlss(pres ~ 1 + (1|yr), sigma.formula=pb(N), data=pmin.spring)
      M1 = gammlss(pres ~ scale(mean.EVI) * scale(mean.t10m) * scale(mean.elev) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M2 = gammlss(pres ~ scale(mean.t10m) * scale(mean.elev) + scale(mean.EVI) * scale(mean.elev) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M3 = gammlss(pres ~ scale(mean.EVI) + scale(mean.t10m) + scale(mean.elev) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M4 = gammlss(pres ~ scale(mean.EVI) + scale(mean.t10m) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M5 = gammlss(pres ~ scale(mean.EVI) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      M6 = gammlss(pres ~ scale(mean.t10m) + (1|yr) + (1|window), sigma.formula=pb(N), data=dat)
      
      GAIC(MNull, M1, M2, M3, M4, M5, M6) #returns matrix with AIC and df 
      
      # return Cox Snell R2, R-squared = 1-(L(0)/L(fitted))^(2/n);  http://rpackages.ianhowson.com/cran/gamlss/man/Rsq.html
      MNull.rsq= Rsq(MNull, type="Cox Snell")
      M1.rsq = Rsq(M1, type="Cox Snell")
      M2.rsq = Rsq(M2, type="Cox Snell")
      M3.rsq = Rsq(M3, type="Cox Snell")
      M4.rsq = Rsq(M4, type="Cox Snell")
      M5.rsq = Rsq(M5, type="Cox Snell")
      M6.rsq = Rsq(M6, type="Cox Snell")
    } 
  }
  
  #Whatever set of models was run, find a way to write out the model eval step and the R2 values into a big table
  #save the modesl to the stats folder
  #print the table to the screen, and save to a txt file to see later.
  mod.eval = cbind(R2, GAIC)
  return(mod.eval)
}


###########################################################################################
#             LOOP THROUGH THE SPECIES AND DATASETS
###########################################################################################

#list of species codes
species <- c("bchu", "bthu", "cahu", "ruhu", "rthu")

for (sp in species){
  spfiles = list.files(path = agan.dir, pattern = glob2rx(paste0(sp,"*.RData")), recursive=FALSE, full.names=TRUE)
  mfiles = list.files(path = migtime.dir, pattern = glob2rx(paste0(sp,"*migration*")), recursive=FALSE, full.names=TRUE)
  
  #import migration dates from previous analysis
  if(sp == "rthu"){ migdates = read.table(mfiles[1], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")}
  if(sp %in% c("bchu", "bthu", "cahu", "ruhu")){ migdates = read.table(mfiles[2], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")}
  migdates=cleanColNames(migdates)
  
  #load all the files into a list
  abs = importANDformat(spfiles[1], 0, migdates, alpha_window)
  pres = importANDformat(spfiles[2], 1, migdates, alpha_window)
  min.05 = importANDformat(spfiles[5], -5, migdates, alpha_window)
  min.10 = importANDformat(spfiles[3], -10, migdates, alpha_window)
  min.15 = importANDformat(spfiles[4], -15, migdates, alpha_window)
  pls.05 = importANDformat(spfiles[8], 5, migdates, alpha_window)
  pls.10 = importANDformat(spfiles[6], 10, migdates, alpha_window)
  pls.15 = importANDformat(spfiles[7], 15, migdates, alpha_window)
  
  #subset data for models - note that models will be run on +/- 15 day comparisons
  pa = subset(rbind(pres, abs), season %in% c("spring", "breeding", "fall"))
  pmin = subset(rbind(pres, min.15), season %in% c("spring", "breeding", "fall"))
  ppls = subset(rbind(pres, pls.15), season %in% c("spring", "breeding", "fall"))
  
  pa.spring = subset(rbind(pres, abs), season=="spring")
  pmin.spring = subset(rbind(pres, min.15), season=="spring")
  ppls.spring = subset(rbind(pres, pls.15), season=="spring")
  
  pa.fall = subset(rbind(pres, abs), season=="fall")
  pmin.fall = subset(rbind(pres, min.15), season=="fall")
  ppls.fall = subset(rbind(pres, pls.15), season=="fall")
  
  #subset data for models - note that models will be run on +/- 15 day comparisons of MEANS
  pa.spring.mean = RS_means(pa.spring)
  pa.fall.mean = RS_means(pa.fall)
  
  pmin.spring.mean = RS_means(pmin.spring)
  pmin.fall.mean = RS_means(pmin.fall)
  
  ppls.spring.mean = RS_means(ppls.spring)
  ppls.fall.mean = RS_means(ppls.fall)
  
  #-----------------------Run hypothesis models and compare for best fit
  gamlss.list <- list()
  print(sp)
  
  #----- all data points
  eval_gamlss_models(pa.spring, sp=sp, lag=FALSE)
  eval_gamlss_models(pa.fall, sp=sp, lag=FALSE)
    
  eval_gamlss_models(pmin.spring, sp=sp, lag=TRUE)
  eval_gamlss_models(pmin.fall, sp=sp, lag=TRUE)
  
  eval_gamlss_models(ppls.spring, sp=sp, lag=TRUE)
  eval_gamlss_models(ppls.fall, sp=sp, lag=TRUE)
  
  #----- means
  eval_gamlss_models(pa.spring.mean, sp=sp, lag=FALSE, means=TRUE)
  eval_gamlss_models(pa.fall.mean, sp=sp, lag=FALSE, means=TRUE)
  
  eval_gamlss_models(pmin.spring.mean, sp=sp, lag=FALSE, means=TRUE)
  eval_gamlss_models(pmin.fall.mean, sp=sp, lag=FALSE, means=TRUE)
  
  eval_gamlss_models(ppls.spring.mean, sp=sp, lag=FALSE, means=TRUE)
  eval_gamlss_models(ppls.fall.mean, sp=sp, lag=FALSE, means=TRUE)
  
  
  #Identify best model using AICc. Assess pvalues and R2 for best model(s).
  #Look at tables in mod.eval output. 
  #Look at list of model output saved as gamlss.list files in the stats folder.
  
  
  