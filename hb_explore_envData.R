# this code is to analyze and plot the eBird hummingbird presence data with environmental correlates
# (c) 2014-2015 by Sarah R. Supp, Tina Cormier, Laura Graham

# import libraries
library(randomForest)
library(sm)
library(ggplot2)
library(ltm)
library(Rarity)
library(nlme)
library(lme4)
library(lubridate)
library(arm)
library(ggmap)
library(GGally)
library(mgcv)
library(gamm4)
library(gamlss)

# define pathnames
function.dir <- "/home/sarah/Documents/GitHub/hb-migration/hb_RS_functions.R"
agan.dir <- "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_annotated_raw/combined/"
migtime.dir <- "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"
fig.dir <- "/home/sarah/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/"
stat.dir <- "/home/sarah/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/stat-tests/"

# pathnames on SB PC
#function.dir <- "hb_RS_functions.R"
#agan.dir <- "D:/hb-migration-files/combined"
#migtime.dir <- "D:/hb-migration-files/aggregate_by_species/"
#fig.dir <- "D:/hb-migration-files/figures/"
#stat.dir <- "D:/hb-migration-files/stat-tests/"

# Laura pathnames
#function.dir <- "hb_RS_functions.R"
#agan.dir <- "~/Dropbox/hb_migration_data/ebird_annotated_raw/combined/"
#migtime.dir <- "~/Dropbox/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"
#fig.dir <- "~/Dropbox/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/"
#stat.dir <- "~/Dropbox/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/stat-tests/"

#function.dir <- "/Users/tcormier/Documents/scripts/git_repos/hb-migration/hb_RS_functions.R"
#agan.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/downloaded_annotations/"
#migtime.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/supp_migration/"
#fig.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/figures/ruhu/"

source(function.dir)

#assign based on the time frame (number of days) used to compute the alpha hulls
alpha_window = 5

###########################################################################################
#SRS started working in this section 5/30/15
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
  #datlist <- sapply(spfiles, function(x) get(load(x)), simplify = FALSE) 
  abs = importANDformat(spfiles[1], 0, migdates, alpha_window)
  pres = importANDformat(spfiles[2], 1, migdates, alpha_window)
  min.05 = importANDformat(spfiles[5], -5, migdates, alpha_window)
  min.10 = importANDformat(spfiles[3], -10, migdates, alpha_window)
  min.15 = importANDformat(spfiles[4], -15, migdates, alpha_window)
  pls.05 = importANDformat(spfiles[8], 5, migdates, alpha_window)
  pls.10 = importANDformat(spfiles[6], 10, migdates, alpha_window)
  pls.15 = importANDformat(spfiles[7], 15, migdates, alpha_window)
  
  #subset dat for plots
  pminpls = subset(rbind(pres, min.05, min.10, min.15, pls.05, pls.10, pls.15), season %in% c("spring", "fall", "breeding"))
  all_ssn = subset(rbind(abs, pres, min.05, min.10, min.15, pls.05, pls.10, pls.15), season %in% c("spring", "fall", "breeding"))
  pres_ssn = subset(pres, season %in% c("spring", "fall", "breeding"))
  abs_ssn = subset(abs, season %in% c("spring", "fall", "breeding"))
  min_ssn = subset(min.15, season %in% c("spring", "fall", "breeding"))
  pls_ssn = subset(pls.15, season %in% c("spring", "fall", "breeding"))
  
  #subset data for glmm / gamm
  pa = subset(rbind(pres, abs), season %in% c("spring", "breeding", "fall"))
  pmin = subset(rbind(pres, min.15), season %in% c("spring", "breeding", "fall"))
  ppls = subset(rbind(pres, pls.15), season %in% c("spring", "breeding", "fall"))
  
  pa.spring = subset(rbind(pres, abs), season=="spring")
  pmin.spring = subset(rbind(pres, min.15), season=="spring")
  ppls.spring = subset(rbind(pres, pls.15), season=="spring")
  
  pa.fall = subset(rbind(pres, abs), season=="fall")
  pmin.fall = subset(rbind(pres, min.15), season=="fall")
  ppls.fall = subset(rbind(pres, pls.15), season=="fall")
  
  #data summaries
  #NOTE: yday1 does not map directly to averages that should be compared, but compare.win does
  all.sum=summarySE(all_ssn, measurevar="EVI", groupvars=c("pres", "yday1", "year", "compare.win"), na.rm=TRUE, conf.interval=0.95)
  all.sum2=summarySE(all_ssn, measurevar="EVI", groupvars=c("pres", "window", "compare.win"), na.rm=TRUE, conf.int=0.95)
  #----------------------- PLOT THE DATA
  # plot summary data 
  
  # Use 95% confidence interval or SE
  noabs = subset(all.sum, pres %in% c(-15, -10, -5, 1, 5, 10, 15))
  noabs$pres = ordered(noabs$pres, levels=c(-15, -10, -5, 1, 5, 10, 15))
  ggplot(data=noabs, aes(x=yday1, y=EVI)) +
      geom_errorbar(aes(ymin=EVI-ci, ymax=EVI+ci), width=0.1, col=gray) + #position=pd #dodge position 
      geom_errorbar(aes(ymin=EVI-se, ymax=EVI+se), width=0.1, col=gray) +
      geom_line(aes(group=compare.win, col=pres, size=N)) +
      scale_colour_manual(values=c("navy", "royalblue3", "turquoise2", "springgreen4", "violetred1", "tomato3", "red4")) +
    facet_wrap(~year)
  
  noabs2 = subset(all.sum2, pres %in% c(-15, -10, -5, 1, 5, 10, 15))
  noabs2$pres = ordered(noabs$pres, levels=c(-15, -10, -5, 1, 5, 10, 15))
  breedwindow=count(unique(pres[pres$season=="breeding", c(18,1)])$window)
  ggplot(data=noabs2, aes(x=window, y=EVI))  +
    geom_vline(data=breedwindow, aes(xintercept=x, size=freq), col="red") +
    geom_errorbar(aes(ymin=EVI-1.96*sd, ymax=EVI+1.96*sd), width=0.1, col=gray, data=noabs2[noabs2$pres==1,]) + #position=pd #dodge position 
    #geom_errorbar(aes(ymin=EVI-se, ymax=EVI+se), width=0.1, col=gray) +
    geom_line(aes(group=compare.win, col=pres, size=N)) + 
    scale_colour_manual(values=c("navy", "royalblue3", "turquoise2", "springgreen4", "violetred1", "tomato3", "red4"))
  
  
  
  
#   # The errorbars overlapped, so use position_dodge to move them horizontally
#   pd <- position_dodge(0.1) # move them .05 to the left and right
#   ggplot(data=all.sum, aes(x=yday1, y=EVI, size=N, col=pres)) + 
#     geom_errorbar(aes(ymin=EVI-ci, ymax=EVI+ci), width=.1, position=pd, col="gray") +
#     geom_line(position=pd,alpha=0.5, aes(group=compare.win)) + geom_point(position=pd,alpha=0.5) + 
#     facet_wrap(~year)
  
  #plot environmental patterns for locations that birds were seen at (7 connected dots for location trajectories)
  png(file.path(path=paste0(fig.dir,sp,"/"), filename=paste0(sp,"_EVI.png")), height=7.5, width=10, units="in", res=300)
  ggplot(pminpls, aes(yday,EVI)) + geom_point(data=abs, aes(yday, EVI)) + geom_point(data=pres_ssn, aes(col=Npp_1km), alpha=0.75) + 
    geom_line(aes(group=id, col=EVI), alpha=0.25) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") +
    scale_color_gradient(low="moccasin", high="chartreuse4") +
    xlab("Julian Day of the Year") + ylab("EVI") + theme_bw() +  
    geom_vline(data=migdates, aes(xintercept=peak_lat)) +
    facet_wrap(~year)
  dev.off()
  
  png(file.path(path=paste0(fig.dir,sp,"/"), filename=paste0(sp,"_Temp.png")), height=7.5, width=10, units="in", res=300)  
  ggplot(pminpls, aes(yday, t10m-273.15)) + geom_point(data=abs, aes(yday, t10m-273.15)) + geom_point(data=pres_ssn, aes(col=t10m-273.15), alpha=0.75) + 
    geom_line(aes(group=id,col=t10m-273.15), alpha=0.25) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") +  
    scale_color_gradient(low="blue", high="firebrick") +
    xlab("Julian Day of the Year") + ylab("Temperature (C)") + theme_bw() + 
    geom_vline(data=migdates, aes(xintercept=peak_lat)) +
    facet_wrap(~year)
  dev.off()
  
  png(file.path(path=paste0(fig.dir,sp,"/"), filename=paste0(sp,"_SRTMelev.png")), height=7.5, width=7.5, units="in", res=300)
  ggplot(pminpls, aes(yday, SRTM_elev)) + geom_point(data=abs, aes(yday, SRTM_elev)) + geom_point(data=pres, aes(col=EVI),alpha=0.5) + 
    #geom_line(aes(group=id, col=EVI), alpha=0.25) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") +
    scale_color_gradient(low="moccasin", high="chartreuse4") + 
    xlab("Julian Day of the Year") + ylab("SRTM Elevation") + theme_bw() + 
    geom_vline(data=migdates, aes(xintercept=peak_lat)) +
    facet_wrap(~year)
  dev.off()
  
  #pairs plots, colored by season
  evitemp = ggplot(pres_ssn, aes(t10m-273.15, EVI, col=season)) + 
    geom_point(alpha=0.25, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  elevtemp = ggplot(pres_ssn, aes(SRTM_elev, t10m-273.15, col=season)) + 
    geom_point(alpha=0.25, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  latevi = ggplot(pres_ssn, aes(location.lat, EVI, col=season)) + 
    geom_point(alpha=0.25, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  lonevi = ggplot(pres_ssn, aes(location.long, EVI, col=season)) + 
    geom_point(alpha=0.25, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  # make a North America base map
  noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")
  spmap = ggmap(noam) + geom_point(data=pres, aes(location.long, location.lat, col=season), alpha=0.25) + 
    scale_color_manual(values=c("black", "orange", "cadetblue", "purple")) + ggtitle(sp)
  
  png(file.path(path=paste0(fig.dir,sp,"/"), filename=paste0(sp,"_multiplot.png")), height=7.5, width=7.5, units="in", res=300)
  multiplot(spmap, evitemp, elevtemp, NA, latevi, lonevi,  cols=2)
  dev.off()
  
  #pairs plot of environmental variables
  png(file.path(path=paste0(fig.dir,sp,"/"), filename=paste0(sp,"_pairsplot.png")), height=15, width=20, units="in", res=300)
  ggpairs(all_ssn[,c(3,2,5:7,10,12:14,21,17)])
  dev.off()
  
  #plot the distribution of the presence and absence data across years
  png(file.path(path=paste0(fig.dir,sp,"/"), filename=paste0(sp,"_numobs.png")), height=15, width=20, units="in", res=300)
  ggplot(pa, aes(yday)) + geom_histogram(aes(fill=pres), alpha=0.5, binwidth=14) + facet_wrap(~year) + theme_bw() + 
    geom_vline(data=migdates,aes(xintercept=peak_lat))
  dev.off()
  
  
  # # Code to make  triangle plots - should plot with means?
  # require(ggtern)
  # ggtern(pres_ssn, aes(x=scale(EVI), y=scale(t10m), z=scale(SRTM_elev), group=window, color=season)) +
  #   geom_confidence(aes(col=season)) +
  #   geom_point(alpha=0.5) +
  #   geom_polygon(alpha=0.5,size=2) +
  #   geom_path(aes(color=season),linetype=2,size=0.5) +
  #   tern_limits(T=1,L=1,R=1) +
  #   theme_bw(base_size = 16) + 
  #   scale_color_manual(values=c("black", "orange", "cadetblue")) + facet_wrap(~year)
  
  #--------------------------------------- STATISTICAL TESTS
  glmm.list <- list()
  print(sp)
  #Is there an environmental signal in the population's immediate region for presence?
  glmm.list[["pa.spring"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=pa.spring)
  glmm.list[["pa.fall"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=pa.fall)
  glmm.list[["pa.all"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=pa)
  
  #Is there an environmental signal for where species are in their migration pathway?
  glmm.list[["pmin.spring"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=pmin.spring)
  glmm.list[["pmin.fall"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=pmin.fall)
  glmm.list[["pmin.all"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=pmin)
  
  glmm.list[["ppls.spring"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=ppls.spring)
  glmm.list[["ppls.fall"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=ppls.fall)
  glmm.list[["ppls.all"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year), family="binomial", data=ppls)
  save(glmm.list, file=paste0(stat.dir, "/", sp, ".glmm.list.rda"))
  
  # there might be an effort difference between windows - include as random factor (models the same otherwise)  
  glmm.window.list <- list()
  #Is there an environmental signal in the population's immediate region for presence?
  glmm.window.list[["pa.spring"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=pa.spring)
  glmm.window.list[["pa.fall"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=pa.fall)
  glmm.window.list[["pa.all"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=pa)
  
  #Is there an environmental signal for where species are in their migration pathway?
  glmm.window.list[["pmin.spring"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=pmin.spring)
  glmm.window.list[["pmin.fall"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=pmin.fall)
  glmm.window.list[["pmin.all"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=pmin)
  
  glmm.window.list[["ppls.spring"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=ppls.spring)
  glmm.window.list[["ppls.fall"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=ppls.fall)
  glmm.window.list[["ppls.all"]] <- glmer(pres ~ scale(EVI) + scale(t10m) + scale(SRTM_elev) + (1|year) + (1|window), family="binomial", data=ppls)
  save(glmm.window.list, file=paste0(stat.dir, "/", sp, ".glmm.window.list.rda"))
  
  
  # glmm vs. gamm? 
  # the same model structures are being repeated but as GAMMs. All fixed effects
  # are currently being smoothed - will need to check outputs to see if this is
  # necessary/realistic
  gamm.list <- list()
  
  try(gamm.list[["pa.spring"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = pa.spring))
  try(gamm.list[["pa.fall"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = pa.fall))
  try(gamm.list[["pa.all"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = pa))
  
  try(gamm.list[["pmin.spring"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = pmin.spring))
  try(gamm.list[["pmin.fall"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = pmin.fall))
  try(gamm.list[["pmin.all"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = pmin))
  
  try(gamm.list[["ppls.spring"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = ppls.spring))
  try(gamm.list[["ppls.fall"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = ppls.fall))
  try(gamm.list[["ppls.all"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year), data = ppls))
  
  save(gamm.list, file=paste0(stat.dir, "/", sp, ".gamm.list.rda"))
  
  gamm.window.list <- list()
  
  try(gamm.window.list[["pa.spring"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = pa.spring))
  try(gamm.window.list[["pa.fall"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = pa.fall))
  try(gamm.window.list[["pa.all"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = pa))
  
  try(gamm.window.list[["pmin.spring"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = pmin.spring))
  try(gamm.window.list[["pmin.fall"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = pmin.fall))
  try(gamm.window.list[["pmin.all"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = pmin))
  
  try(gamm.window.list[["ppls.spring"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = ppls.spring))
  try(gamm.window.list[["ppls.fall"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = ppls.fall))
  try(gamm.window.list[["ppls.all"]] <- gamm4(pres ~ s(EVI) + s(t10m) + s(SRTM_elev), family="binomial", random=~(1|year) + (1|window), data = ppls))
  
  save(gamm.window.list, file=paste0(stat.dir, "/", sp, ".gamm.window.list.rda"))
  
  #compare years? - current feeling is that this might be difficult due to
  #varying effort across years. Need to look into
  
}

#-------------------- FIXME: SRS stopped editing here (July 31, 2015).

###########################################################################################

###########################################################################################
