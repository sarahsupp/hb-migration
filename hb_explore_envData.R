# this code is to analyze and plot the eBird hummingbird presence data with environmental correlates

# import libraries
library(randomForest)
#library(party)
library(sm)
library(ggplot2)
library(ltm)
library(Rarity)
library(lme4)
library(lubridate)
library(arm)

# define pathnames
function.dir <- "/home/sarah/Documents/GitHub/hb-migration/hb_RS_functions.R"
agan.dir <- "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_annotated_raw/combined/"
migtime.dir <- "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"
fig.dir <- "/home/sarah/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/"

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
  
  #load all the files into a list
  #datlist <- sapply(spfiles, function(x) get(load(x)), simplify = FALSE) 
            #TODO: Update this so that it can import ALL needed files - Files [3] and [6] correspond to +/- 10 day lag
  abs = get(load(spfiles[1]))
  pres = get(load(spfiles[2]))
  min = get(load(spfiles[3]))
  pls = get(load(spfiles[6]))

  if(sp == "rthu"){ migdates = read.table(mfiles[1], header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")}
  else{ migdates = read.table(mfiles[2], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")}
  
  #strip extra "" and "\\" from fields
  names(migdates) = gsub('.{1}$', '', substring(names(migdates),3))   #make the column names look nicer
  migdates$spr_begin = as.numeric(gsub("\"", "", migdates$spr_begin, fixed=TRUE))
  migdates$peak_lat = as.numeric(gsub("\"", "", migdates$peak_lat, fixed=TRUE))
  migdates$fal_end = as.numeric(gsub("\"", "", migdates$fal_end, fixed=TRUE))
  migdates$spr_spd = as.numeric(gsub("\"", "", migdates$spr_spd, fixed=TRUE))
  migdates$fal_spd = as.numeric(gsub("\"", "", migdates$fal_spd, fixed=TRUE))
  migdates$species = gsub("\"", "", migdates$species, fixed=TRUE)
  migdates$year = as.factor(gsub("\"", "", migdates$year, fixed=TRUE))

  #add code for present, absent, and present, but forward (2) or backward (-1) in time (# days in time window)
  abs$pres = as.factor(0)
  pres$pres = as.factor(1)
  min$pres = as.factor(-1) #note the time window for alpha hulls is -5, so this represents -5, -10 or -15 days from matching pres points, with the same locations
  pls$pres = as.factor(2)  #note the time window for alpha hulls is +5 from matching pres points, so this represents +5, +10, or +15 days, with the same locations
  
  #add rownames as a column for later matching
  abs$id = as.factor(row.names(abs))
  pres$id = as.factor(row.names(pres))
  min$id = as.factor(row.names(min))
  pls$id = as.factor(row.names(pls))
  
  all = rbind(pres, abs, min, pls)
  
  #lubridate to pull month and year from filename
  all$month = as.factor(month(as.Date(all$timestamp)))
  all$year =  as.factor(year(as.Date(all$timestamp)))
  all$yday = as.numeric(yday(as.Date(all$timestamp)))
  
  pres = subset(all, pres==1)
  abs = subset(all, pres==0)
  min = subset(all, pres==-1)
  pls = subset(all, pres==2)
  
  #assign time windows based on presence and absence data, and label by season
  pres = assign_window_season(migdates, pres, alpha_window)
  abs = assign_window_season(migdates, abs, alpha_window)
  min = assign_window_season(migdates, min, alpha_window)
  pls = assign_window_season(migdates, pls, alpha_window)
  
  pres$season = as.factor(pres$season)
  abs$season = as.factor(abs$season)
  min$season = as.factor(min$season) #forward in time
  pls$season = as.factor(pls$season) #backward in time

  pminpls = rbind(pres, min, pls)
  pminpls = subset(pminpls, season %in% c("spring", "fall", "breeding"))
  pres_ssn = subset(yrwindows, season %in% c("spring", "fall", "breeding"))
  
  #plot environmental patterns for locations that birds were seen at (3 connected dots for location trajectories)
  ggplot(pminpls, aes(yday, EVI)) + geom_point(data=abs, aes(yday, EVI)) + geom_point(aes(col=EVI), alpha=0.5) + 
    geom_line(aes(group=id, col=EVI), alpha=0.5) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") +
    scale_color_gradient(low="moccasin", high="chartreuse4") +
    xlab("Julian Day of the Year") + ylab("EVI") + theme_bw() + 
    facet_wrap(~year)
  
  ggplot(pminpls, aes(yday, t10m-273.15)) + geom_point(data=abs, aes(yday, t10m-273.15)) + geom_point(aes(col=t10m-273.15), alpha=0.5) + 
    geom_line(aes(group=id,col=t10m-273.15), alpha=0.5) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") +  
    scale_color_gradient(low="blue", high="firebrick") +
    xlab("Julian Day of the Year") + ylab("Temperature (C)") + theme_bw() +
    facet_wrap(~year)
  
  ggplot(pminpls, aes(yday, SRTM_elev)) + geom_point(data=abs, aes(yday, SRTM_elev)) + geom_point(aes(col=rh2m),alpha=0.5) + geom_line(aes(group=id, col=rh2m)) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") +
    scale_color_gradient(low="moccasin", high="navy") + 
    xlab("Julian Day of the Year") + ylab("SRTM Elevation") + theme_bw() + 
    facet_wrap(~year)
  
  ggplot(pminpls, aes(yday, swrf)) + geom_point(data=abs, aes(yday, swrf)) + geom_point(aes(col=rh2m),alpha=0.5) + geom_line(aes(group=id, col=rh2m)) +
    stat_smooth(data=pres_ssn, method="gam", formula = y~s(x, k=10), col="black") + facet_wrap(~year)+  
    scale_color_gradient(low="moccasin", high="navy") + theme_bw()
  
  pa = subset(rbind(pres, abs), season %in% c("spring", "breeding", "fall"))
  pmin = subset(rbind(pres, min), season %in% c("spring", "breeding", "fall"))
  ppls = subset(rbind(pres, pls), season %in% c("spring", "breeding", "fall"))
  pres2 = subset(pres, season %in% c("spring", "breeding", "fall"))
  
  #pairs plots, colored by season
  evitemp = ggplot(pres2, aes(t10m-273.15, EVI, col=season)) + 
    geom_point(alpha=0.5, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  elevtemp = ggplot(pres2, aes(SRTM_elev, t10m-273.15, col=season)) + 
    geom_point(alpha=0.5, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  latevi = ggplot(pres2, aes(location.lat, EVI, col=season)) + 
    geom_point(alpha=0.5, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  lonevi = ggplot(pres2, aes(location.long, EVI, col=season)) + 
    geom_point(alpha=0.5, size=3) + #shape=pres, 
    scale_color_manual(values=c("black", "orange", "cadetblue")) + theme_bw() + stat_smooth(method=lm)
  
  # make a North America base map
  noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")
  spmap = ggmap(noam) + geom_point(data=pres, aes(location.long, location.lat, col=season), alpha=0.5) + 
    scale_color_manual(values=c("black", "orange", "cadetblue", "purple")) + ggtitle(sp)
  
  multiplot(spmap, evitemp, elevtemp, NA, latevi, lonevi,  cols=2)
    
  #statistical model
  #Is there an environmental signal in the population's immediate region for presence
  glm1 = glmer(pres ~ scale(EVI) + scale(t10m) + scale(swrf) + scale(SRTM_elev) + (1|year) + (1|season), family="binomial", data=pa)
  glm2 = glmer(pres ~ scale(EVI) + scale(t10m) + scale(swrf) + scale(SRTM_elev) + season + (1|year), family="binomial", data=pa)
  
  glm1 = glmer(pres ~ scale(EVI) + scale(t10m) + scale(swrf) + scale(SRTM_elev) + (1|year) + (1|season), family="binomial", data=pmin)
  glm2 = glmer(pres ~ scale(EVI) + scale(t10m) + scale(swrf) + scale(SRTM_elev) + season + (1|year), family="binomial", data=pmin)
  
  #FIXME: Stopped editing here.
  
  ggplot(pres, aes(as.Date(timestamp), EVI)) + geom_point(aes(col=t10m - 273.15)) + 
    scale_color_gradient(low="blue", high="firebrick")
  
  ggplot(pmin, aes(window, EVI, col=pres)) + geom_point(alpha=0.5) + 
    stat_smooth(method="gam", formula = y ~ s(x, k=20)) + facet_wrap(~year) +
    scale_color_manual(values = c("indianred", "gray")) + theme_bw()
  
  ggplot(ppls, aes(window, EVI, col=pres)) + geom_point(alpha=0.5) + 
    stat_smooth(method="loess") + facet_wrap(~year) +
    scale_color_manual(values = c("indianred", "gray")) + theme_bw()
  
  ggplot(pa, aes(yday, EVI, col=pres)) + geom_point(alpha=0.5) + 
    stat_smooth(method="gam", formula = y ~ s(x, k=20)) + facet_wrap(~year)  +
    scale_color_manual(values = c("indianred", "gray"))
  
  plot(pres$EVI, min$EVI, pch=19)
  abline(a=0, b=1, col="red")
  
  plot(pres$EVI, pls$EVI, pch=19)
  abline(a=0, b=1, col="red")
  
  plot(pres$EVI, abs$EVI, pch=19)
  abline(a=0, b=1, col="blue", lwd=2)
  
  plot(as.numeric(yday(as.Date(pres$timestamp))), as.numeric(yday(as.Date(pls$timestamp))), pch=19)
  
}




###########################################################################################
#spp abbreviation (same one used in file paths and names)
spp <- "ruhu"

#what lag do you want to explore?
lag <- 0


###########################################################################################

#ann file
agan.file <- paste0(agan.dir, spp, "/", spp, "_lag",lag,"_allYears.csv")

ann <- read.csv(agan.file, as.is=T)
ann$timestamp <- as.Date(ann$timestamp, format='%Y-%m-%d %H:%M:%S.000')
ann$presence <- as.factor(ann$presence)
ann$GlobCover <- as.factor(ann$GlobCover)
ann$month <- as.numeric(format(ann$timestamp, "%m"))
ann$year <- as.numeric(format(ann$timestamp, "%Y"))
ann$doy <- as.numeric(format(ann$timestamp, "%j"))


#glob cover labels - legend for map comes from http://dup.esrin.esa.it/files/p68/GLOBCOVER2009_Product_Description_Manual_1.0.pdf
#glob.label <- unique(ann$GlobCover)
#oops, missing some here - go back and figure out which classes are missing.
#glob.names <- c("closed Broadleaved","Closed to Open Grassland","Mosaic Vegetation Cropland","Closed to Open Broadleaved Evergreen and/or Semi-Deciduous Forest",
#                 "Closed Needle-Leaved Evergreen Forest","Urban","Water","Closed to Open Mixed Broadleaved and Needleleaved Forest","Mosaic Forest or Shrubland and Grassland",
#                 "Closed to Open Shrubland","Mosaic Cropland/Vegetation", "Closed Broadleaved Forest Regularly Flooded, Fresh Water",
#                 "Closed Broadleaved Semi-Deciduous and/or Evergreen Forest Regularly Flooded, Saline Water", "Sparse Vegetation",
#                 "Close to Open Grassland or Shurbland or Woody Vegetation on Regularly Flooded or Waterlogged soil, fresh brakish, or saline water",
#                 "Open Broadleaved Deciduous Forest/Woodland","Bare","Permanent Snow and Ice","No Data")


#separate pres from abs
 pres <- ann[ann$presence == 1,]
 abs <- ann[ann$presence == 0,]

#plot doy representation to check that it is similar for pres and abs
ggplot(data=pres, aes(doy)) + geom_histogram(fill="red", alpha=0.5) + 
  geom_histogram(data=abs, aes(doy), fill="blue", alpha=0.5) + 
  theme_classic() + ggtitle("blue = absent, red = present")


#Separate by season
#first read in the migration timing table and format
migtime.file <- paste0(migtime.dir, spp, "/", "west_migration_", spp, ".txt")
migtime <- read.table(migtime.file, header=T, sep=" ")
years <- unique(ann$year)

#set up dataframes that will hold data during spring and fall migrations 
# Note: data that falls outside of seasons is not included in the new dataframe (sf)
sf <- as.data.frame(matrix(data=NA, nrow=0, ncol=ncol(ann)+1))
names(sf) <- c(names(ann), "season")

for (yr in years) {
  seasonal <- seasonalSub(ann, migtime, yr)
  #append this year's data to sf table
  sf <- rbind(sf, seasonal)
}#end yr loop


###########################################################################################
# habitat utilization graphs (i.e. variable histograms)

#convert K to celsius
sf$t10m <- sf$t10m - 273.15
sf$Temp_sfc <- sf$Temp_sfc - 273.15

sf$occurrence[sf$presence==1] <- "present"
sf$occurrence[sf$presence==0] <- "absent"

#subset by presence, absence
pres <- sf[sf$occurrence=="present",]
abs <- sf[sf$occurrence=="absent",]

#subset by spring, fall
spr <- sf[sf$season == "spring",]
fall <- sf[sf$season == "fall",]

# Kernel density plots for spring and fall migrations
#Set up variables and titles for looping
seasons <- c("spr", "fall")
season.titles <- c("Spring", "Fall")

vars <- c("t10m", "Temp_sfc", "Total_precipitation_sfc", "EVI", "swrf", "lwrf", "u10m", "uplift", "SRTM_elev")
var.titles <- c("Temperature at 10 m", "Surface Temperature", "Surface Precipitation", "EVI", "Downward Shortwave Radiation Flux", 
                "Downward longwave Radiation Flux", "East-West Wind at 10 m", "Uplift", "Elevation")
xlab.titles <- c(expression("Temperature at 10 m" ~ (degree~C)), expression("Surface Temperature" ~ (degree~C)), "Precipitation (mm)",
                 "EVI", expression("Downward Shortwave Radiation Flux" ~ (W ~ m^-2)), expression("Downward Longwave Radiation Flux" ~ (W ~ m^-2)), 
                  "East-West Wind at 10 m (m/s)", "Uplift", "Elevation (m)")

#for each variable, let's look at presence vs. absence in spring and fall
for (i in 1:length(vars)) {
  #print(paste(vars[i]))
  for (s in 1:length(seasons)) {
    #print(paste(vars[i], seasons[s]))
    #create graphs for spring presence/absence, spring first
    title <- paste0(season.titles[s], " Rufous Hummingbird Habitat Utilization - \n", var.titles[i])
    outfile <- paste0(fig.dir, season.titles[s], "_ruhu_habitat_utilization_", vars[i], ".pdf")
    spr.col <- c("gray30","cyan2")
    fall.col <- c("gray30", "orange")
    
    if (seasons[s] == "spr") {
      col <- spr.col
    } else {
      col <- fall.col
    }#end color if
    
    pdf(outfile, width=9, height=8)
    p <- ggplot(get(seasons[s]), aes(x=get(vars[i]),fill=presence)) + geom_density(alpha=.5)
    p <- p + scale_fill_manual( values = col)
    p <- p + theme_classic() + theme(text=element_text(size=14))
    p <- p + ggtitle(title) + xlab(xlab.titles[i])
    print(p)
    
    dev.off()
   }#end seasons loop
} #end vars loop


# Now, look at presence only and the difference between seasons on the same plot
for (j in 1:length(vars)) {
  #create graphs of presence spring vs fall
  title <- paste0("Spring vs. Fall Rufous Hummingbird Habitat Utilization - \n", var.titles[j])
  outfile <- paste0(fig.dir, "ruhu_springVsFall_habitat_utilization_", vars[j], ".pdf")
  pdf(outfile, width=9, height=8)
 
  p <- ggplot(pres, aes(x=get(vars[j]), fill=season)) + geom_density(alpha=.5)   
  p <- p + scale_fill_manual( values = c("orange", "cyan2"))
  p <- p + theme_classic() + theme(text=element_text(size=14))
  p <- p + ggtitle(title) + xlab(xlab.titles[j])
  print(p)
  dev.off()
} #end vars loop


###########################################################################################
#pairwise correlation plot
#first, we only want presence/absence and the environmental variables
keepcols <- c("Temp_sfc", "Total_precipitation_sfc", "EVI", "swrf", "lwrf", "u10m", "uplift", "SRTM_elev" )
df.cor <- sf[,names(sf) %in% keepcols]
df.cor <- df.cor[complete.cases(df.cor),]

out <- paste0(fig.dir, spp, "_variable_correlations.png")
png(out, width=12, height=12, units="in", res=600)
pairs(df.cor, diag.panel=panel.hist, upper.panel=panel.cor)
dev.off()


###########################################################################################
#Random forest - for variable importance
#TODO: need to filter NDVI/EVI for quality - just looking quickly at this now.
#first remove columns with bad preds (for some reason, rugosity5, slopeu, and slovev produced all NAs)
preds <- ann[,-c(4,7,12:14)]
preds <- na.omit(preds)

preds <- preds[,-1]
response <- preds$presence

#keep just the environmental vars here
preds.rf <- preds[,c(2,3,7,8,9,10,11,12,13,15,16,18,21,23,24,35,36,37)]
rf <- randomForest(preds.rf, response, na.rm=T, importance=T)

# Use variable importance plot based on mean decrease in GINI. Calle & Urrea (2010) show its rankings are more stable.
#   Explanation here: http://www.stat.berkeley.edu/~breiman/RandomForests/cc%5Fhome.htm#varimp
varImp <- importance(rf,type = 2)

# Plot results
rf.graph <- paste0(fig.dir, spp, "/", spp, "_varImp_allvars.pdf")
pdf(file = rf.graph, width=8, height=8)
par(mar=c(8,6,6,4))
varImpPlot(rf, type=2,pch=16, col="blue",main=paste(spp, "\nRandom Forest Variable Importance", sep=' '))
dev.off()


###########################################################################################
# Compare environmental data using Kolmogorov-Smirnoff (KS) tests
# KS test has less power to detect a shift in the median but more power to detect changes in 
# the shape of the distributions. Null hypothesis is that both groups were sampled from populations 
# with identical distributions. It tests for any violation of that null hypothesis -- 
# different medians, different variances, or different distributions.

critical_D <- function(n1, n2){
  # http://www.soest.hawaii.edu/wessel/courses/gg313/Critical_KS.pdf
  # assuming alpha = 0.05
  # values > Da are significant?
  Da <- 1.36 * sqrt((n1+n2)/(n1*n2))
  return(Da)
}

seasons <- c("spring", "fall")
vars <- c("SRTM_elev", "EVI", "lwrf", "swrf", "t10m") #add other vars as necessary

# Compare presence vs absence points in seasons and years for each variable
ks_pa <- data.frame("year"=1, "season"=NA, "var"=NA, "Dstat"=1, "Pvalue"=1)
for (y in unique(years)){
  for (s in unique(seasons)){
    pres <- sf[sf$present == 1 & sf$year == y & sf$season == s,]
    abs <- sf[sf$present == 0 & sf$year == y & sf$season == s,]
    for (var in unique(vars)){
      compare <- ggplot(pres, aes(x=get(var))) + geom_density(alpha=0.3, fill="red") + 
        geom_density(data=abs, alpha=0.3, fill="blue") + theme_classic() + xlab(var) + 
        ggtitle(paste(y,s,var, sep=" - "))
      print(compare)
      ks <- ks.test(pres[,colnames(pres) %in% var], abs[,colnames(abs) %in% var])
      ks_pa = rbind(ks_pa, c(y,s,var,round(as.numeric(ks$statistic),4),round(ks$p.value,4)))
    }  
  }
}
ks_pa <- ks_pa[-1,] #delete first row of dummy data


# compare years. Are there any years that really differ?
ks_yrs <- data.frame("year1"=1, "year2"=1, "var"=NA, "Dstat"=1, "pvalue"=1, "signif"=NA)
for (var in unique(vars)){
  for (y in unique(years)){
    if(y==2013)
      next
    for(ynext in unique(years)){
      if (y > ynext)
        next
      if(y==ynext)
        next
    pres <- sf[sf$present == 1 & sf$year == y,]
    abs <- sf[sf$present == 0 & sf$year == y,]
    pres2 <- sf[sf$present == 1 & sf$year == ynext,]
    abs2 <- sf[sf$present == 0 & sf$year == ynext,]
    compare <- ggplot(pres, aes(x=get(var))) + geom_density(alpha=0.3, fill="red") + 
      geom_density(data=pres2, alpha=0.3, fill="blue") + 
      geom_density(data=abs, alpha=0.3, fill="grey60",col="red") + 
      geom_density(data=abs2, alpha=0.3, fill="grey60",col="blue") + 
      theme_classic() + xlab(var) + 
      ggtitle(paste(y,ynext,var, sep=" - "))
    print(compare)
    ks <- ks.test(pres[,colnames(pres) %in% var], pres2[,colnames(pres2) %in% var])
    Da <- critical_D(nrow(pres), nrow(pres2))
    if(Da > ks$statistic){ sig <- "Y" }
    else{ sig <- "N" }  
    ks_yrs <- rbind(ks_yrs, c(y, ynext, var, round(as.numeric(ks$statistic),4), round(ks$p.value,4), sig))
  }}
}
ks_yrs <- ks_yrs[-1,]


# compare seasons Are hb selecting diff in diff seasons? How influenced by background diffs is this?
ks_season <- data.frame("year"=1, "var"=NA, "Dstat"=1, "pvalue"=1, "signif"=NA)
for (var in unique(vars)){
  for (y in unique(years)){
      spr_pres <- sf[sf$present == 1 & sf$year == y & sf$season == "spring",]
      fal_pres <- sf[sf$present == 1 & sf$year == y & sf$season == "fall",]
      compare <- ggplot(spr_pres, aes(x=get(var))) + geom_density(alpha=0.3, fill="cadetblue") + 
        geom_density(data=fal_pres, alpha=0.3, fill="orange") + 
        theme_classic() + xlab(var) + ggtitle(paste(y, var, sep=" - "))
      print(compare)
      ks <- ks.test(spr_pres[,colnames(spr_pres) %in% var], fal_pres[,colnames(fal_pres) %in% var])
      Da <- critical_D(nrow(pres), nrow(pres2))
      if(Da > ks$statistic){ sig <- "N" }
      else{ sig <- "Y" }  
      ks_season <- rbind(ks_season, c(y, var, round(as.numeric(ks$statistic),4), round(ks$p.value,4), sig))
    }
}
ks_season <- ks_season[-1,]



#---------------------------------------------------------------------------------
#         GLM test comparing distribution of envr. data with presence
#---------------------------------------------------------------------------------
seasons <- c("spring", "fall")

for (y in unique(years)){
  for (s in unique(seasons)){
    data <- sf[sf$year == y & sf$season == s,]
    
    fit <- glm(presence ~ swrf + lwrf + SRTM_elev + EVI + Temp_sfc, Total_precipitation_sfc, family="binomial", data=data)
    
    #check the residual deviance and degrees of freedom in the model (should not be significant)
    check <- 1 - pchisq(fit$deviance, fit$df.residual)
    if (check > 0.05) {
      print(paste(y, s, "WAS SUCCESSFUL! RESULTS BELOW:", sep=" "))
      print (summary(fit))
    }
    else { print("WE NEED A BETTER MODEL") }
  }
}


#---------------------------------------------------------------------------------------------------
#         GLMER test comparing distribution of envr. data with presence with year as random effect
#---------------------------------------------------------------------------------------------------

# get rid of columns that aren't needed
sf2 <- sf[complete.cases(sf[,-4]),]

#standardize the data
zscore = apply(sf2[,names(sf2)%in% vars ], 2, function(x) {
  y = (x - mean(x))/sd(x)
  return(y)
})

#make into data frame and and ID and time vars back in
zscore <- as.data.frame(zscore)
zscore$presence <- sf2$presence
zscore$year <- as.factor(sf2$year)
zscore$season <- sf2$season


fit <- glmer(presence ~ Temp_sfc + t10m + swrf + lwrf + EVI + Total_precipitation_sfc + uplift + SRTM_elev + 
               (1|year) + (1|season), data = sf, family = "binomial", 
                control = glmerControl(optimizer = "bobyqa"))
summary(fit)



#-------------------------- An interesting (or confusing?) way to plot the data using violin plots
tmp <- melt(sf[, c("present", "EVI", "SRTM_elev", "t10m", "swrf")], id.vars="present")
ggplot(tmp, aes(x = as.factor(present), y = value)) +
  geom_jitter(alpha = .1) +
  geom_violin(alpha = .75) +
  facet_grid(variable ~ .) +
  scale_y_sqrt()


###########################################################################################
# Correlate graphs, omitting NA values
ann.sub <- ann[,c(6,8,9,11,15:ncol(ann))]
corr.names <- c("elev","NDVI","EVI","Rugosity","uplift","lwrf","swrf","temp10m","uwind10m","vwind10m")
yy <- cor(ann.sub, use = "pairwise.complete.obs")

pres.sub <- pres[,c(6,8,9,11,15:ncol(pres))]
pp <- cor(pres.sub, use = "pairwise.complete.obs")

abs.sub <- abs[,c(6,8,9,11,15:ncol(abs))]
aa <- cor(abs.sub, use = "pairwise.complete.obs")

rownames(yy) <- corr.names
colnames(yy) <- corr.names
rownames(pp) <- corr.names
colnames(pp) <- corr.names
rownames(aa) <- corr.names
colnames(aa) <- corr.names

# Matrix plot of variable correlations ordered by principle components
circle.corr( yy, order = TRUE, bg = "white",
             col = colorRampPalette(c("blue","white","red"))(100) )

circle.corr( pp, order = TRUE, bg = "white",
             col = colorRampPalette(c("blue","white","red"))(100) )

circle.corr( aa, order = TRUE, bg = "white",
             col = colorRampPalette(c("blue","white","red"))(100) )




###########################################################################################
#Now try variable importance testing with the party package based on this blog http://alandgraf.blogspot.com/2012/07/random-forest-variable-importance.html and the Strobl et al. paper)
#WOW, this ran for almost 24 hrs and crashed my machine. Need to rethink.
# data.controls <- cforest_unbiased(ntree=1000, mtry=round(sqrt(ncol(preds))))
# cf1 <- cforest(present~.,data=preds,control=data.controls)
# varimp(cf1)
# varimp(cf1,conditional=TRUE)


#Variables can be considered informative and important if their variable importance value is above the absolute value of the lowest negative-scoring variable. 
#“The rationale for this rule of thumb is that the importance of irrelevant variables varies randomly around zero” (Strobl et al. 2009b: 342). 
# barplot(sort(varImp), horiz=TRUE, xlab="Variable Importance \n(predictors to right of red dashed are significant)", ))
# abline(v=abs(min(varImp)), col='red',lty='longdash', lwd=2)
# 
# #Variables can be considered informative and important if their variable importance value is above the absolute value of the lowest negative-scoring variable. 
# #“The rationale for this rule of thumb is that the importance of irrelevant variables varies randomly around zero” (Strobl et al. 2009b: 342). 
# dotplot(sort(data.cforest.varimp), xlab=”Variable Importance
#         in DATA\n(predictors to right of dashed vertical line are
#                   significant)”, panel = function(x,y){
#                     panel.dotplot(x, y, col=’darkblue’, pch=16, cex=1.1)
#                     panel.abline(v=abs(min(data.cforest.varimp)), col=’red’,
#                                  lty=’longdash’, lwd=2)