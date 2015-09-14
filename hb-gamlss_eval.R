# this code is to run gamlss models for the eBird hummingbird data with environmental correlates
# (c) 2015 by Sarah R. Supp, Laura Graham, and Frank La Sorte

# import libraries
library(ggplot2)
library(nlme)
library(lme4)
library(lubridate)
library(arm)
library(mgcv)
library(gamm4)
library(gamlss)

# dropbox pathname
dropbox <- "/home/lorra/Dropbox/"
#dropbox <- "/home/sarah/Dropbox/Hummingbirds/"

# define pathnames (FAL will need to change pathnames)
function.dir <- "hb_RS_functions.R"
agan.dir <- paste0(dropbox, "hb_migration_data/ebird_annotated_raw/combined/")
migtime.dir <- paste0(dropbox, "hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/")
fig.dir <- paste0(dropbox, "NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/")
stat.dir <- paste0(dropbox, "NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/stat-tests/")

source(function.dir) 

#assign based on the time frame (number of days) used to compute the alpha hulls
alpha_window = 5

#list of species codes
species <- c("bchu", "bthu", "cahu", "ruhu", "rthu")

###########################################################################################
#             FUNCTION TO EVALUATE GAMLSS MODELS
###########################################################################################

eval_gamlss_models = function(dat=dat, sp=sp, season="season", lag=FALSE, means=FALSE, output.dir=stat.dir) {
  # evaluates gamlss models for each dataset based on the hypothesis testing framework
  # models are run on centered and scaled (such that mean=0, sd=1) predictor variables
  # models input year and window (alpha hull) as random effects
  # if lag==FALSE, models will be run with the set of hypotheses for within-region signal of presence
  # if lag==TRUE, models will be run with the set of hypotheses for migration trajectory with timelag
  # if means==FALSE, the data represents all the records from the alpha hulls
  # if means==TRUE, the data has been averaged by window (each window represents a unique alpha hull)
  # outputs model results, AIC comparison, R2, and plots the data
  # top model(s) are written to output.dir
  
  #convert year to numeric for the models
  #dat$year = as.numeric(dat$year)
  
  if(!means){
    cat("data represents all records from the alpha hulls \n")
    
    #only include columns for the means (NAs in rows with N<2 give errors in gamlss)
    dat = dat[,names(dat) %in% c("pres", "year", "window", "location.lat", "location.long", "EVI", "t10m", "SRTM_elev")]
    dat = na.omit(dat) #remove any lingering rows with NA
    
    #scale data for the models
    dat$EVI <- as.vector(scale(dat$EVI))
    dat$t10m <- as.vector(scale(dat$t10m))
    dat$SRTM_elev <- as.vector(scale(dat$SRTM_elev))
    dat$year <- as.factor(dat$year)
    dat$window <- as.factor(dat$window)
    
    #drop the levels for windows with no data
    dat$window <- droplevels(dat$window)
    
    #-----For within-region signal of presence (presence vs. absence points)
    if(!lag){
      cat("models will compare within region presence vs. absence \n")
      
      # H0a. Remotely sensed variables will not be predictive of hummingbird presence, 
      #     because they do not adequately capture the processes (resources, demand) 
      #     that determine where birds are found along their migration route. 
      MNull = gamlss(pres ~ 1, family="BI", data=dat)
      
      # H1. In the spring, birds should be more strongly tied to EVI, assuming that greenness 
      #     is a better proxy for food resources in the spring than in the fall.
      # H2. In the fall birds should be constrained to temperatures and Elevation*EVI interaction (resources), 
      #     particularly in the west where there is more extreme topography and habitat heterogeneity.
      M.ETV = gamlss(pres ~ pb(EVI, df=5) * pb(t10m, df=5) * pb(SRTM_elev, df=5) + random(year) + random(window), family="BI", data=dat)
      M.E_T_V = gamlss(pres ~ pb(EVI, df=5) + pb(t10m, df=5) + pb(SRTM_elev, df=5) + random(year) + random(window), family="BI", data=dat)
      M.E_T = gamlss(pres ~ pb(EVI, df=5) + pb(t10m, df=5) + random(year) + random(window), family="BI", data=dat)
      M.E = gamlss(pres ~ pb(EVI, df=5) + random(year) + random(window), family="BI", data=dat)
      M.T = gamlss(pres ~ pb(t10m, df=5) + random(year) + random(window), family="BI", data=dat)
      
      #AIC comparison of models, returns matrix with AIC and df, ordered with best model first 
      AIC.df = GAIC(MNull, M.ETV, M.E_T_V, M.E_T, M.E, M.T) 
      # Calculate AIC weights (Hobbs and Hilborn 2006)
      AIC.df$w = round(Akaike.weight(AIC.df$AIC),2)
      
      # save model summaries for those with weight > 0
      gamlss.list = list()
      sig.models = row.names(AIC.df[AIC.df$w>0,])
      for (mod in 1:length(sig.models)) { gamlss.list[[sig.models[mod]]] = get(sig.models[mod]) }
      save(gamlss.list, file=paste0(output.dir, "/", sp, ".", season, ".gamlss.list.rda"))
      rm(gamlss.list)
    }
    
    #-----For migratory trajectory (presence vs +/- 15 day points)
    if(lag){
      cat("models will compare presence with migratory trajectory lags (+/- 15 days \n")
      
      # H0b. Remotely sensed variables will not be predictive of hummingbird presence, 
      #     because they do not adequately capture the processes (resources, demand) 
      #     that determine where birds are found along their migration route.
      MNull = gamlss(pres ~ 1, family="BI", data=dat)
      
      # H3. In the spring, birds will be more strongly tied to EVI “greening up” 
      #     as they move northward to breeding grounds than in the fall.
      # H4. In the fall, birds will be more limited by temperature, as late summer/early fall 
      #     temperatures tend to be higher (e.g. closer to lethal limit) than in the spring.
      # H5. In the fall, birds should be associated with higher elevations than in the spring, 
      #     because flower phenology (blooming time) is later at higher altitudes 
      #     (so it is associated with more resources) and cooler temperatures are found at higher elevations.
      M.ETV = gamlss(pres ~ pb(EVI, df=5) * pb(t10m, df=5) * pb(SRTM_elev, df=5) + random(year) + random(window), family="BI", data=dat)
      M.EV_TV = gamlss(pres ~ pb(t10m, df=5) * pb(SRTM_elev, df=5) + pb(EVI, df=5) * pb(SRTM_elev, df=5) + random(year) + random(window), family="BI", data=dat)
      M.E_T_V = gamlss(pres ~ pb(EVI, df=5) + pb(t10m, df=5) + pb(SRTM_elev) + random(year) + random(window), family="BI", data=dat)
      M.E_T = gamlss(pres ~ pb(EVI, df=5) + pb(t10m, df=5) + random(year) + random(window), family="BI", data=dat)
      M.E = gamlss(pres ~ pb(EVI, df=5) + random(year) + random(window), family="BI", data=dat)
      M.T = gamlss(pres ~ pb(t10m, df=5) + random(year) + random(window), family="BI", data=dat)
      
      #AIC comparison of models, returns matrix with AIC and df, ordered with best model first 
      AIC.df = GAIC(MNull, M.ETV, M.EV_TV, M.E_T_V, M.E_T, M.E, M.T) 
      # Calculate AIC weights (Hobbs and Hilborn 2006)
      AIC.df$w = round(Akaike.weight(AIC.df$AIC),2)
      
      # save model summaries for those with weight > 0
      gamlss.list = list()
      sig.models = row.names(AIC.df[AIC.df$w>0,])
      for (mod in 1:length(sig.models)) { gamlss.list[[sig.models[mod]]] = get(sig.models[mod]) }
      save(gamlss.list, file=paste0(output.dir, "/", sp, ".", season, ".gamlss.list.rda"))
      rm(gamlss.list)
    }
  }
  
  if(means){
    cat("data has been averaged within unique alpha hulls (window) \n")
    cat("effort will be accounted for as number of checklists where species was present \n")
    
    #only include columns for the means (NAs in rows with N<2 give errors in gamlss)
    dat = dat[,names(dat) %in% c("pres", "yday1", "year", "window", "N", "mean.lat", "mean.lon", "mean.EVI", "mean.t10m", "mean.elev")]
    dat = na.omit(dat) #remove any lingering rows with NA
    
    #scale data for the models
    dat$mean.EVI <- as.vector(scale(dat$mean.EVI))
    dat$mean.t10m <- as.vector(scale(dat$mean.t10m))
    dat$mean.SRTM_elev <- as.vector(scale(dat$mean.elev))
    dat$year <- as.factor(dat$year)
    dat$window <- as.factor(dat$window)
    
    #drop the levels for windows with no data
    dat$window <- droplevels(dat$window)
    
    #-----For within-region signal of presence (presence vs. absence points)
    if(!lag){
      cat("models will compare within region presence vs. absence \n")
      
      #models use the same basic hypotheses as above, for no lag, but include effort 
      #models are labeled E=EVI, T=Temp, V=Elev
      MNull = gamlss(pres ~ 1 + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat) # null model doesn't seem to work when include random factors
      M.ETV1 = gamlss(pres ~ pb(mean.EVI, df=5) * pb(mean.t10m, df=5) * pb(mean.elev, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.E_T_V = gamlss(pres ~ pb(mean.EVI, df=5) + pb(mean.t10m, df=5) + pb(mean.elev, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.E_T = gamlss(pres ~ pb(mean.EVI, df=5) + pb(mean.t10m, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.E = gamlss(pres ~ pb(mean.EVI, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.T = gamlss(pres ~ pb(mean.t10m, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      
      #AIC comparison of models, returns matrix with AIC and df, ordered with best model first 
      AIC.df = GAIC(MNull, M.ETV, M.E_T_V, M.E_T, M.E, M.T) 
      # Calculate AIC weights (Hobbs and Hilborn 2006)
      AIC.df$w = round(Akaike.weight(AIC.df$AIC),2)
      
      # save model summaries for those with weight > 0
      gamlss.list = list()
      sig.models = row.names(AIC.df[AIC.df$w>0,])
      for (mod in 1:length(sig.models)) { gamlss.list[[sig.models[mod]]] = get(sig.models[mod]) }
      save(gamlss.list, file=paste0(output.dir, "/", sp, ".", season, ".gamlss.means.list.rda"))
      rm(gamlss.list)
    }
    
    #-----For migratory trajectory (presence vs +/- 15 day points)
    if(lag){
      cat("models will compare presence with migratory trajectory lags (+/- 15 days \n")
      
      #models use the same basic hypotheses as above, for +/- 15 day lag, but include effort
      #models are labeled E=EVI, T=Temp, V=Elev, _ indicates additive
      MNull = gamlss(pres ~ 1, sigma.fix=TRUE, family="BI", data=dat) #null model doesn't seem to work if include random factors
      M.ETV = gamlss(pres ~ pb(mean.EVI, df=5) * pb(mean.t10m, df=5) * pb(mean.elev, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.EV_TV = gamlss(pres ~ pb(mean.t10m, df=5) * pb(mean.elev, df=5) + pb(mean.EVI, df=5) * pb(mean.elev, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.E_T_V = gamlss(pres ~ pb(mean.EVI, df=5) + pb(mean.t10m, df=5) + pb(mean.elev) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.E_T = gamlss(pres ~ pb(mean.EVI, df=5) + pb(mean.t10m, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.E = gamlss(pres ~ pb(mean.EVI, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      M.T = gamlss(pres ~ pb(mean.t10m, df=5) + random(year) + random(window), sigma.formula=pb(N), family="BI", data=dat)
      #AIC comparison of models, returns matrix with AIC and df, ordered with best model first 
      AIC.df = GAIC(MNull, M.ETV, M.EV_TV, M.E_T_V, M.E_T, M.E, M.T) 
      # Calculate AIC weights (Hobbs and Hilborn 2006)
      AIC.df$w = round(Akaike.weight(AIC.df$AIC),2)
      
      # save model summaries for those with weight > 0
      gamlss.list = list()
      sig.models = row.names(AIC.df[AIC.df$w>0,])
      for (mod in 1:length(sig.models)) { gamlss.list[[sig.models[mod]]] = get(sig.models[mod]) }
      save(gamlss.list, file=paste0(output.dir, "/", sp, ".", season, ".gamlss.means.list.rda"))
      rm(gamlss.list)
    }
  }
  #Whatever set of models was run, find a way to write out the model eval step and the R2 values into a big table
  #save the modesl to the stats folder
  #print the table to the screen, and save to a txt file to see later.

  # return Cox Snell R2, R-squared = 1-(L(0)/L(fitted))^(2/n);  http://rpackages.ianhowson.com/cran/gamlss/man/Rsq.html
  AIC.df$Rsq = rep(NA,nrow(AIC.df))
  for(r in 1:length(row.names(AIC.df))) { r2=round(Rsq(get(row.names(AIC.df)[r]), type="Cox Snell"),4); AIC.df$Rsq[r]=r2 }
  AIC.df$nobs = rep(NA,nrow(AIC.df))
  for(r in 1:length(row.names(AIC.df))) { nobs=get(row.names(AIC.df)[r])$N; AIC.df$nobs[r]=nobs }
  
  return(AIC.df)
}


###########################################################################################
#             LOOP THROUGH THE SPECIES AND DATASETS
###########################################################################################

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
  
  print (paste0("Imported data for species: ", sp))
  
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
  #NOTE: The labels seem to be backwards (pmin when matched on window or yday1 represents where hb will be, ppls represents where hb were)
  pa.spring.mean = RS_means(pa.spring)
  pa.fall.mean = RS_means(pa.fall)
  
  pmin.spring.mean = RS_means(pmin.spring)
  pmin.fall.mean = RS_means(pmin.fall)
  
  ppls.spring.mean = RS_means(ppls.spring)
  ppls.fall.mean = RS_means(ppls.fall)
  
  ############################################################################
  #-----------------------Run hypothesis models and compare for best fit
  ############################################################################
  mod.list <- list()
  print(sp)
  
  #----- all data points
  mod.list[["pa.spring"]] = eval_gamlss_models(pa.spring, sp=sp, season="spring", lag=FALSE, means=FALSE, output.dir=stat.dir)
  mod.list[["pa.fall"]] = eval_gamlss_models(pa.fall, sp=sp, season="fall", lag=FALSE, means=FALSE, output.dir=stat.dir)
    
  mod.list[["pmin.spring"]] = eval_gamlss_models(pmin.spring, sp=sp, season="springmin", lag=TRUE, means=FALSE, output.dir=stat.dir)
  mod.list[["pmin.fall"]] = eval_gamlss_models(pmin.fall, sp=sp, season="fallmin", lag=TRUE, means=FALSE, output.dir=stat.dir)
  
  mod.list[["ppls.spring"]] = eval_gamlss_models(ppls.spring, sp=sp, season="springpls", lag=TRUE, means=FALSE, output.dir=stat.dir)
  mod.list[["ppls.fall"]] = eval_gamlss_models(ppls.fall, sp=sp, lag=TRUE, season="fallpls", means=FALSE, output.dir=stat.dir)
  
  #----- means
  mod.list[["pa.mean.spring"]] = eval_gamlss_models(pa.spring.mean, season="spring", sp=sp, lag=FALSE, means=TRUE, output.dir=stat.dir)
  mod.list[["pa.mean.fall"]] = eval_gamlss_models(pa.fall.mean, sp=sp, season="fall", lag=FALSE, means=TRUE, output.dir=stat.dir)
  
  mod.list[["pmin.mean.spring"]] = eval_gamlss_models(pmin.spring.mean, season="springmin", sp=sp, lag=TRUE, means=TRUE, output.dir=stat.dir)
  mod.list[["pmin.mean.fall"]] = eval_gamlss_models(pmin.fall.mean, sp=sp, season="fallmin", lag=TRUE, means=TRUE, output.dir=stat.dir)
  
  mod.list[["ppls.mean.spring"]] = eval_gamlss_models(ppls.spring.mean, season="springpls", sp=sp, lag=TRUE, means=TRUE, output.dir=stat.dir)
  mod.list[["ppls.mean.fall"]] = eval_gamlss_models(ppls.fall.mean, sp=sp, season="fallpls", lag=TRUE, means=TRUE, output.dir=stat.dir)
 
  save(mod.list, file=paste0(stat.dir, "/", sp, ".mod.list.rda"))
  
  print (paste0("Finished with models for ", sp))
}

  
  
  