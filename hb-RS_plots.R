# this code is to run make figures for the eBird hummingbird data with environmental correlates
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
#dropbox <- "/home/lorra/Dropbox/"
dropbox <- "/home/sarah/Dropbox/Hummingbirds/"

# github pathname
github <- "/home/sarah/Documents/GitHub/hb-migration/"

# define pathnames (FAL will need to change pathnames)
function.dir <- paste0(github, "hb_RS_functions.R")
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
  
  spring.mean = rbind(pmin.spring.mean, ppls.spring.mean[ppls.spring.mean$pres==15,])
  fall.mean = rbind(pmin.fall.mean, ppls.fall.mean[ppls.fall.mean$pres==15,])
  
  means = rbind(spring.mean, fall.mean)
  
  ############################################################################
  #-----------------------Visualize data
  ############################################################################
  print(paste0("making plots for ", sp))
  
  ggplot(means, aes(yday1, mean.elev)) + geom_line(aes(col=factor(pres,labels=c("present","+15 days", "15 days ago")))) + 
    geom_vline(data=migdates, aes(xintercept=peak_lat)) + facet_wrap(~year) + 
    xlab("Julian day") + ylab("mean Elevation (meters)") + 
    scale_colour_manual(values=c("black", "indianred", "deepskyblue3")) +
    labs(col="migration route") + ggtitle(sp)
  ggsave(file=paste0(fig.dir, sp, "/mean_elev.png"))
  
  ggplot(means, aes(yday1, mean.EVI)) + geom_line(aes(col=factor(pres,labels=c("present","+15 days", "15 days ago")))) + 
    geom_vline(data=migdates, aes(xintercept=peak_lat)) + facet_wrap(~year) + 
    xlab("Julian day") + ylab("mean EVI") + 
    scale_colour_manual(values=c("black", "indianred", "deepskyblue3")) +
    labs(col="migration route") + ggtitle(sp)
  ggsave(file=paste0(fig.dir, sp, "/mean_EVI.png"))
  
  ggplot(means, aes(yday1, mean.t10m-273.15)) + geom_line(aes(col=factor(pres,labels=c("present","+15 days", "15 days ago")))) + 
    geom_vline(data=migdates, aes(xintercept=peak_lat)) + facet_wrap(~year) + 
    xlab("Julian day") + ylab("mean temperature (Celsius)") + 
    scale_colour_manual(values=c("black", "indianred", "deepskyblue3")) +
    labs(col="migration route") + ggtitle(sp)
  ggsave(file=paste0(fig.dir, sp, "/mean_temp.png"))
  
}

  
  