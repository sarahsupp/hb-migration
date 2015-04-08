#Code for eBird migration project
# (c) 2013 -2014 Sarah Supp 

library(ggmap)
library(maptools)
library(fields)
library(sp)
library(raster)
library(maps)
library(mapdata)
library(rgdal)
library(raster)
library(gamm4)

#set working directory
main = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration"
figpath = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/Figures"
gitpath = "C:/Users/sarah/Documents/GitHub/hb-migration"
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data"
setwd(wd)


#---------------------------------------------------------------------------------------
#              predict migration paths, dates, speed, and error across years
#---------------------------------------------------------------------------------------

# read in summary of effort data (Number of eBird checklists submitted per day per year)
effort = read.table("FAL_hummingbird_data/checklist_12_2004-2013wh_grp.txt", header=TRUE, as.is=TRUE)

# read in country outline shp files
USAborder = readShapePoly("borders/USA_adm/USA_adm0.shp")
Mexborder = readShapePoly("borders/MEX_adm/MEX_adm0.shp")
Canborder = readShapePoly("borders/CAN_adm/CAN_adm0.shp")

# read in the altitude layers, make rasters
elev = raster("alt_5m_bil/alt.bil")

# plot elev + map for extent
myext <- c(-175, -50, 15, 75)
plot.new()
#plot(elev, ext = myext, xlab="Longitude", ylab = "Latitude", xlim = c(-175,-50), ylim = c(15,75), col=gray(0:256/256))

borders <- function(){
  plot(USAborder, ext=myext, border="black", add=TRUE)
  plot(Mexborder, ext=myext, border="black", add=TRUE)
  plot(Canborder, ext=myext, border="black", add=TRUE)
}

plot(elev, ext=myext, addfun=borders, ylab="Latitude", xlab="Longitude", xlim = c(-175,-50), ylim = c(15,75), col=gray(0:256/256)) 

# read in the north america equal area hex grid map (FAL) and format for use
# other options include a quad map (terr_4h6/nw_vector_grid.shp) or a hexmap with land only (terr_4h6/terr_4h6.shp", sep="")
hexgrid = readShapePoly(paste(main, "/data/icosahedron_land_and_sea/icosahedron.shp", sep="")) #hex with land and sea, cropped to North America

# plot the hexgrid on a map of north america
plot(NA, NA, xlim=c(-175, -50), ylim=c(15, 75), xlab = "Longitude", ylab = "Latitude")
map('worldHires', c("usa", "canada", "mexico"), add=TRUE, fill=T, col="lightblue")
plot(hexgrid, add=T)

# make a North America base map
noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")

# map total birder effort 2004:2013
eft = PlotChecklistMap(effort, hexgrid, wd)
rm(eft)

# read in eBird data
files = list.files(pattern = "*.txt")

#for each eBird file, print the number sightings per year and per month.
#plot the locations of sightings on a map, color coded by month
for (f in 1:length(files)){
  
  require(ggmap)
  require(ggplot2)
  require(plyr)
  require(reshape2)
  require(Rmisc)
  require(sp)
  require(raster)
  require(segmented)
  source(paste(gitpath, "/migration-fxns.r", sep=""))
  
  humdat = read.table(files[f], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  names(humdat) = c("SCI_NAME", "PRIMARY_COM_NAME","YEAR", "DAY", "TIME", "GROUP_ID", "PROTOCOL_ID",
                    "PROJ_ID", "DURATION_HRS", "EFFORT_DISTANCE_KM", "EFFORT_AREA_HA", "NUM_OBSERVERS",
                    "LATITUDE", "LONGITUDE", "SUB_ID", "POLYFID", "MONTH")
  
  humdat$MONTH = factor(humdat$MONTH, levels=c(1:12), ordered=TRUE)
  
  #grab species name for setting directory paths and naming figures
  species = humdat$SCI_NAME[1] 
  species = gsub(" ","", species, fixed=TRUE)
  species = gsub("\"", "", species, fixed=TRUE)
  
  # set years of data to use - data after 2007 is more reliable
  years = c(2004:2013)

  #start a new directory
  dirpath = paste(figpath, "/", species, sep="")
     #dir.create(dirpath, showWarnings = TRUE, recursive = FALSE) #only need if directory did not previously exist
  
  #show how many records there are for the species across the years, write to txt file
  #yeartable = PlotRecords(humdat$YEAR, species)
  #ggsave(file=paste(dirpath, "/", "years", species, ".pdf", sep=""))
  
  write.table(yeartable, file = paste(dirpath, "/", species,".txt",sep=""), row.names=FALSE)
  
  #save a figure of the geographic number of checklists for the species, over all the years
  #count = PlotChecklistMap(humdat, hexgrid, dirpath)
  
  for (y in 1:length(years)){
    yrdat = humdat[which(humdat$YEAR == years[y]),]
    yreffort = effort[which(effort$YEAR == years[y]),]
    
    #plot frequency of sightings per month
    #monthtable = PlotRecords(yrdat$month, species)
    
    #get daily weighted mean location 
    meanlocs = AlternateMeanLocs(yrdat, species, hexgrid, yreffort)
    
    #use GAM model to predict daily location along a smoothing line
    preds = EstimateDailyLocs(meanlocs)
    
    #use gam approach to estimate rough starting points for segmentation from the mean loc latitude data
    startpoints = round(Est3MigrationDates(meanlocs))
    migration = startpoints
    
    #save a plot of the species migration route mapped onto continent with real observations
    pdf(file = paste(dirpath, "/trimmed-route_", species, years[y], ".pdf", sep=""), width = 7, height = 4.5)
    BasePlotMigration(preds, yrdat, migration, elev, USAborder, Mexborder, Canborder, myext)
    dev.off() 
    
    #save a plot of the species migration mapped onto an elevation raster
    pdf(file = paste(dirpath, "/elev-route_", species, years[y], ".pdf", sep=""), width = 7, height = 4.5)
    ElevPlotMigration(preds, yrdat, migration, elev, USAborder, Mexborder, Canborder, myext)
    dev.off() 
    
    #get Great Circle distances traveled each day between predicted daily locations
    dist = DailyTravel(preds, 4, 5, species, years[y], migration)
    
    #estimate migration speed for spring and fall
    speed = MigrationSpeed(dist, migration)
    
    #plot smoothed migration trajectory for the species and year
    #mig_path = PlotMigrationPath(preds, noam, species, years[y])
    #ggsave(mig_path, file=paste(dirpath, "/", "migration", species, years[y], ".pdf", sep=""))
    
#     #plot occurrences with lines showing beginning and end of migration
     PlotOccurrences(meanlocs, species, years[y], migration)
     ggsave(file=paste(dirpath, "/", "occurrences", species, years[y], ".tiff", sep=""), dpi=600)

    #Subset western species by flyway data (check bias in SE US data points) Sensu La Sorte et al in press - 
        #"The role of atmospheric conditions in the seasonal dynamics of North American migration flyways" - JOurnal of Biogeography
    if (species %in% c("Archilochusalexandri", "Selasphorusplatycercus", "Selasphorusrufus", "Selasphoruscalliope")) {
          
      #only use data west of the 103rd meridian (western flyway)
      west_yrdat = yrdat[which(yrdat$LONGITUDE <= -103),]
      
      #get daily weighted mean location 
      west_meanlocs = AlternateMeanLocs(west_yrdat, species, hexgrid, yreffort)
      
      #use GAM model to predict daily location along a smoothing line
      west_preds = EstimateDailyLocs(west_meanlocs)
      
      #use gam approach to estimate rough starting points for segmentation from the mean loc latitude data
      west_migration = round(Est3MigrationDates(west_meanlocs))
      
      #get Great Circle distances traveled each day between predicted daily locations
      west_dist = DailyTravel(west_preds, 4, 5, species, years[y], west_migration)
      
      #estimate migration speed for spring and fall
      west_speed = MigrationSpeed(west_dist, west_migration)
      
      #save a plot of the species migration route mapped onto continent with real observations
      pdf(file = paste(dirpath, "/WEST_trimmed-route_", species, years[y], ".pdf", sep=""), width = 7, height = 4.5)
      BasePlotMigration(west_preds, west_yrdat, west_migration, elev, USAborder, Mexborder, Canborder, myext)
      dev.off() 
      
      #save a plot of the species migration mapped onto an elevation raster
      pdf(file = paste(dirpath, "/WEST_elev-route_", species, years[y], ".pdf", sep=""), width = 7, height = 4.5)
      ElevPlotMigration(west_preds, west_yrdat, west_migration, elev, USAborder, Mexborder, Canborder, myext)
      dev.off() 
    }

    else {
      #create dummy variables for the eastern species. In this case, west == all variables
      west_preds = preds
      west_migration = migration
      west_speed = speed
    }

    #add year to preds, so we can save it to compare across years
    preds$year = years[y]
    west_preds$year = years[y]
    
    if (y == 1){
      pred_data = preds
      migdates = data.frame("spr_begin" = migration[[1]], "spr_end" = migration[[2]], 
                            "fal_begin" = migration[[2]], "fal_end" = migration[[3]],
                            "species" = species, "year" = years[y])
      migspeed = data.frame("spr" = speed[1], "fal" = speed[2], "species" = species, "year" = years[y])
      
      west_pred_data = west_preds
      west_migdates = data.frame("spr_begin" = west_migration[1], "spr_end" = west_migration[2], 
                                 "fal_begin" = west_migration[2], "fal_end" = west_migration[3],
                                 "species" = species, "year" = years[y])
      west_migspeed = data.frame("spr" = west_speed[1], "fal" = west_speed[2], "species" = species, "year" = years[y])
    }
    else{
      pred_data = rbind(pred_data, preds)
      dates = c(migration[[1]], migration[[2]], migration[[2]], migration[[3]], "species" = species, year = years[y])
      speed = c(speed, species, years[y])
      migdates = rbind(migdates, dates)
      migspeed = rbind(migspeed, speed)
      
      west_pred_data = rbind(west_pred_data, west_preds)
      west_dates = c(west_migration[[1]], west_migration[[2]], west_migration[[2]], west_migration[[3]], "species" = species, year = years[y])
      west_speed = c(west_speed, species, years[y])
      west_migdates = rbind(west_migdates, west_dates)
      west_migspeed = rbind(west_migspeed, west_speed)
    }
    
    if (y == length(years)){
      #write migration timing and speed data to file
      write.table(migdates, file = paste(getwd(), "/output_data/migration", species, ".txt", sep=""), 
                  append=FALSE, row.names=FALSE)
      write.table(migspeed, file = paste(getwd(), "/output_data/speed", species, ".txt", sep=""), 
                  append=FALSE, row.names=FALSE)
      write.table(pred_data, file = paste(getwd(), "/output_data/centroids", species, ".txt", sep=""), 
                  append=FALSE,row.names=FALSE)
      
      #write western flyway migration timing and speed data to file
      write.table(west_migdates, file = paste(getwd(), "/output_data/west_migration", species, ".txt", sep=""), 
                  append=FALSE, row.names=FALSE)
      write.table(west_migspeed, file = paste(getwd(), "/output_data/west_speed", species, ".txt", sep=""), 
                  append=FALSE, row.names=FALSE)
      write.table(west_pred_data, file = paste(getwd(), "/output_data/west_centroids", species, ".txt", sep=""), 
                  append=FALSE,row.names=FALSE)
    }
    
    rm(list=ls()[ls() %in% c("sitemap", "meanmap", "yrdat", "meanlocs", "migration", "preds", "dist", "mig_path", "speed",
                             "west_migration", "west_preds", "west_dist", "west_speed")])   # clears the memory of the map and year-level data
  }
  rm(list=ls()[!ls() %in% c("f", "files", "noam", "hexgrid", "effort", "pred_data", "migdates", 
                            "wd", "main", "gitpath","figpath", "USAborder", "Mexborder", "Canborder", "elev", "myext")])   # clears the memory of everything except the file list, iterator, and base map
}


#---------------------------------------------------------------------------------------
#                   compare migration rates and dates across species
#---------------------------------------------------------------------------------------
require(ggmap)
require(ggplot2)
require(plyr)
require(reshape2)
require(Rmisc)
require(sp)
require(raster)
require(gamm4)
require(chron)
require(fields)
require(RColorBrewer)
source(paste(gitpath, "/migration-fxns.r", sep=""))

setwd(wd)

# read in eBird data for migration speed, dates, and predicted path (should be in same order for species)
rfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = "west_speed", full.names=TRUE)
mfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = c("west_migration.*txt"), full.names=TRUE)
cfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = "west_centroids.*.txt", full.names=TRUE)

# species ordered by data files, body sizes from Dunning 2008, migration distance from Nature Serve centroids
species = c("Black-chinned", "Ruby-throated", "Calliope", "Broad-tailed", "Rufous")
mass = c(3.4, 3.1, 2.65, 3.55, 3.5)
distance = c(1721.49, 2765.86, 3252.82, 1737.89, 4102.58)
spdata = data.frame(species, mass, distance, "lat_r2"= NA, "lon_r2" = NA, "spr_speed" = NA,
                    "fal_speed" = NA, "spr_date" = NA, "fal_date" = NA)


#------------------------------------------ ANALYZE THE DATA -------------------------------------

#-------------------- 
#       model to test variance in lat and lon across years, with year as a random effect
#-------------------- get gamm4 results using all years (2004-2013) and later years (2008-2013)
for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  preds_sub = preds[which(preds$year > 2007),]
  print(preds[1,1]) #species name
  years = c(2004:2013)
  
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    between = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday <= dates[y,4]),]
    spring = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday < dates[y,2]),]
    fall = preds[which(preds$year == years[y] & preds$jday > dates[y,2] & preds$jday <= dates[y,4]),]
    
    if (y == 1){
      migpreds = between
      pred_spr = spring
      pred_fal = fall
    }
    else{
      migpreds = rbind(migpreds, between)
      pred_spr = rbind(pred_spr, spring)
      pred_fal = rbind(pred_fal, fall)
    }
  }
  
  # subset the seasonal data by the recent years (2008-2013)
  pred_spr_sub = pred_spr[which(pred_spr$year > 2007),]
  pred_fal_sub = pred_fal[which(pred_fal$year > 2007),]
  
  #gam analysis to determine variance across years, using all of the daily centroids (Jday 1-365) in each year
  lon_gam = gamm4(lon ~ s(jday, k=10), random = ~(1|year), data = preds_sub)
  print (paste("R2 for Longitude is:", round(summary(lon_gam$gam)$r.sq,4)))
  lat_gam = gamm4(lat ~ s(jday, k=10), random = ~(1|year), data=preds_sub)
  print (paste("R2 for Latitude is:", round(summary(lat_gam$gam)$r.sq,4)))
  
  #gam analysis to determine variance across years, using only spring centroids (start of spring to peak latitude date)
  lon_gam_spr = gamm4(lon ~ s(jday, k=10), random = ~(1|year), data=pred_spr_sub)
  print (paste("R2 for spring Longitude is:", round(summary(lon_gam_spr$gam)$r.sq,4)))
  lat_gam_spr = gamm4(lat ~ s(jday, k=10), random = ~(1|year), data=pred_spr_sub)
  print (paste("R2 for spring Latitude is:", round(summary(lat_gam_spr$gam)$r.sq,4)))
  
  #gam analsysis to determine variance across years, using only fall centroids (peak latitude to end of fall date)
  lon_gam_fal = gamm4(lon ~ s(jday, k=10), random = ~(1|year), data=pred_fal_sub)
  print (paste("R2 for fall Longitude is:", round(summary(lon_gam_fal$gam)$r.sq,4)))
  lat_gam_fal = gamm4(lat ~ s(jday, k=10), random = ~(1|year), data=pred_fal_sub)
  print (paste("R2 for fall Latitude is:", round(summary(lat_gam_fal$gam)$r.sq,4)))
  
  spdata[f,4] = summary(lat_gam$gam)$r.sq
  spdata[f,5] = summary(lon_gam$gam)$r.sq
}


#--------------------
#       get linear model results using years (2008-2013) and plot as barplot
#-------------------- 
lm_mig = data.frame("species"="Archilochusalexandri", "season" = "spring", "year" = 1, "lat_slope" = 1, "lat_r2" = 1, "lon_slope" = 1, "lon_r2" = 1)

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  species = preds[1,1] #species name
  years = c(2008:2013)
  
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    between = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday <= dates[y,4]),]
    spring = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday < dates[y,2]),]
    fall = preds[which(preds$year == years[y] & preds$jday > dates[y,2] & preds$jday <= dates[y,4]),]
    
    #get slope & r2 of spring and fall latitudinal migration
    spr_lm = LinearMigration(spring, years[y])
    spr_lm = data.frame("species" = species, "season"= "spring", spr_lm)
    fal_lm = LinearMigration(fall, years[y])
    fal_lm = data.frame("species" = species, "season"= "fall", fal_lm)
    lm_mig = rbind(lm_mig, spr_lm)
    lm_mig = rbind(lm_mig, fal_lm)
  }
}
lm_mig = lm_mig[-1,] #delete first row of dummy data

#plot the variance in estimated migration begin and end for all and for recent years
pdf(file = paste(figpath, "/linearslope_all_species.pdf", sep=""), width = 6, height = 5)

bxp_rate = ggplot(lm_mig, aes(season, abs(lat_slope), fill=season)) + geom_boxplot() + theme_classic() + 
  scale_fill_manual(values=c("cadetblue", "orange"), guide = "none") + 
  ylab("linear slope of seasonal migration") + theme(text = element_text(size=12)) + 
  scale_y_continuous(breaks = seq(0, 0.40, by = 0.10), limits = c(0,0.40)) + theme(text = element_text(size=12)) +
  facet_wrap(~species)

multiplot(bxp_rate, cols = 1)
dev.off()


#--------------------------- 
#       generate figures and table data for migration speed
#--------------------------- boxplots of migration speed for spring vs. fall for each species
rate = data.frame("spr" =1, "fal" = 1, "species" = "none", "year" = 1)
for (f in 1:length(rfiles)){
  sp_rate = read.table(rfiles[f], header=TRUE, sep=" ", fill=TRUE, comment.char="")
#  sp_rate = subset(sp_rate, year>2007) #use only years 2008-2013
  rate = rbind(rate, sp_rate)
  
  #print the mean and standard deviation of speed for each species (note: this includes all years, unless uncomment above line)
  print(sp_rate$species[1])
  print(paste("spring sd:", sd(sp_rate$spr)))
  print(paste("spring mean:", mean(sp_rate$spr)))
  print(paste("fall sd:", sd(sp_rate$fal)))
  print(paste("fall mean:", mean(sp_rate$fal)))
  print("")
  
  spdata[f,6] = mean(sp_rate$spr)
  spdata[f,7] = mean(sp_rate$fal)
}

rate = subset(rate, year > 2007) #subset to better-sampled years
r = melt(rate[,c(1,2,3,4)], id.vars = c("year", "species"))
names(r) = c("year", "species", "season", "rate")

#plot the variance in estimated migration speed for all and for recent years
pdf(file = paste(figpath, "/speed_all_species.pdf", sep=""), width = 15, height = 3)
  
bxp_speed = ggplot(r, aes(season, rate, fill=season)) + geom_boxplot() + theme_bw() + 
    scale_fill_manual(values=c("cadetblue", "orange"), guide="none") + ylab("km/day") + 
    scale_y_continuous(breaks = seq(0, 80, by = 20), limits = c(0,80)) + 
    facet_wrap(~species, nrow=1)
  
multiplot(bxp_speed, cols = 5)
dev.off() 


#--------------------------
#         generate figures and table data for migration dates
# ------------------------ Boxplots of the number of days +/- mean migration date, by species
dates = data.frame("spr_begin" = 1, "mid" = 1, "fal_end" = 1, "species" = "none", "year" = 1)

for (f in 1:length(mfiles)){
  sp_dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
# sp_dates = subset(sp_dates, year>2007)
  
  #Print species table data (note: this includes all years, unless uncomment above line)
  print(sp_dates$species[1])
  print(paste("spring begin median:", median(sp_dates$spr_begin)))
  print(paste("spring begin sd:", sd(sp_dates$spr_begin)))
  print(paste("lat peak median:", median(sp_dates$spr_end)))
  print(paste("lat peak sd:", sd(sp_dates$spr_end)))
  print(paste("fall end median:", median(sp_dates$fal_end)))
  print(paste("fall end sd:", sd(sp_dates$fal_end)))
  print("")
  
  #standardize dates, to get number of days +/- mean
  sp_dates = subset(sp_dates, year > 2007)
  
  spdata[f,8] = round(sd(sp_dates$spr_begin),4)
  spdata[f,9] = round(sd(sp_dates$fal_end), 4)
  
  sp_dates$spr_begin = sp_dates$spr_begin - mean(sp_dates$spr_begin)
  sp_dates$mid = sp_dates$spr_end - mean(sp_dates$spr_end)
  sp_dates$fal_end = sp_dates$fal_end - mean(sp_dates$fal_end)
  dates = rbind(dates, sp_dates[,c(1,7,4,5,6)])

}

dates = dates[-1,]
d = melt(dates, id.vars = c("species", "year"))
names(d) = c("species", "year", "season", "date")
  
#plot the variance in estimated migration begin and end for all and for recent years
pdf(file = paste(figpath, "/migdates_all_species.pdf", sep=""), width = 15, height = 3)
  
bxp_date = ggplot(d, aes(season, date, fill=season)) + geom_boxplot() + theme_bw() + 
  scale_fill_manual(values=c("cadetblue", "olivedrab3", "orange"), guide = "none") + 
  ylab("number of days -/+ mean date") + theme(text = element_text(size=12)) + 
  scale_y_continuous(breaks = seq(-40, 40, by = 20), limits = c(-40,40)) +
  facet_wrap(~species, nrow=1)
  
  multiplot(bxp_date)
  dev.off()


#------------------------------------------ PLOT THE DATA -------------------------------------


#------------------------------------
#       plot the data for species comparisons
#------------------------------------

ggplot(spdata, aes(distance, lat_r2)) + geom_point(size = mass) + xlab("total migration distance") +
  ylab("Latitude R2 by date") + stat_smooth(method = "lm") + theme_classic()

ggplot(spdata, aes(distance, lon_r2)) + geom_point(size = mass) + xlab("total migration distance") +
  ylab("Longitude R2 by date") + stat_smooth(method = "lm") + theme_classic()

ggplot(spdata, aes(distance, spr_speed)) + geom_point(size = mass) + xlab("total migration distance") +
  ylab("Population spring migration speed (km/day)") + 
  stat_smooth(method = "lm", col = "cadetblue", fill = "cadetblue", alpha = 0.2) + 
  theme_classic() #+ geom_text(label=species)

ggplot(spdata, aes(distance, fal_speed)) + geom_point(size = mass) + xlab("total migration distance") +
  ylab("Population fall migration speed (km/day)") + stat_smooth(method = "lm", col = "orange", fill = "orange", alpha = 0.2) +
  theme_classic() #+ geom_text(label=species)

ggplot(spdata, aes(distance, spr_date)) + geom_point(size = mass) + xlab("total migration distance") +
  ylab("sd in spring onset") + stat_smooth(method = "lm", col = "cadetblue", fill = "cadetblue", alpha = 0.2) +
  theme_classic()

ggplot(spdata, aes(distance, fal_date)) + geom_point(size = mass) + xlab("total migration distance") +
  ylab("sd in fall arrival") + stat_smooth(method = "lm", col = "orange", fill = "orange", alpha = 0.2) +
  theme_classic()


#---------------------------------
#       plot standard deviation in the daily centroid estimates across years 2008-2013
#---------------------------------
for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  species = preds[1,1] #species name
  years = c(2008:2013)
  
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    between = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday <= dates[y,4]),]

    if (y == 1){
      migpreds = between
    }
    else{
      migpreds = rbind(migpreds, between)
    }
  }
  
  #compare location across the years using mean and sd for 2008:2013
  jdays = sort(unique(migpreds$jday))
  patherr = data.frame("jday"=1,"meanlat"=1, "sdlat"=1, "meanlon"=1, "sdlon"=1)
  outcount = 1
  for (j in min(jdays):max(jdays)){
    tmp = migpreds[which(migpreds$jday == j),]
    meanlat = mean(tmp$lat)
    sdlat = sd(tmp$lat)
    meanlon = mean(tmp$lon)
    sdlon = sd(tmp$lon)
    patherr[outcount,] = c(j, meanlat, sdlat, meanlon, sdlon)
    outcount = outcount + 1
  }
    
  # plot the standard deviation in daily lat and lon across the 10 years, and across 6 most recent years
  pdf(file = paste(dirpath, "/ErrorinDailyLocs", species, ".pdf", sep=""), width = 5, height = 8)
  
  ymax = max(c(patherr$sdlat, patherr$sdlon),na.rm=TRUE) 
  sdlocs = ggplot(patherr, aes(jday, sdlat)) + geom_point(size=1) + theme_classic() + 
    ggtitle(paste(species, "sd in daily locs", min(years), "-", max(years))) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) + 
    scale_x_continuous(breaks = seq(0, 366, by = 25), limits = c(0, 366)) + 
    geom_point(aes(jday, sdlon), col = "indianred", size=1) + ylab("stdev daily lat (black) and lon (red)")

  lonlatday = ggplot(patherr, aes(sdlon, sdlat, col=jday)) + geom_point(size=1) + theme_classic() +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax))  +
    scale_x_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) 
 
  multiplot(sdlocs, lonlatday, cols = 1)
  dev.off()
}
  
  
#---------------------------------
#       plot all migration routes for a species on a map with elevation raster
#---------------------------------
pdf(file = paste(figpath, "/elev-route_summary_all_species.pdf", sep=""), width = 7, height = 4.5)

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  preds = subset(preds, year > 2007)
  elev = raster("alt_5m_bil/alt.bil")  #elevation layers
  myext = c(-175, -50, 15, 75) #set extent to North America
  species = preds[1,1] #species name
  
  AllMigration(preds, elev, myext, species)
}
dev.off() 


#---------------------------------
#       plot standard error in predicted centroids across years
#---------------------------------
pdf(file=paste(figpath, "/se_corlon-lat_all_species.pdf", sep=""), width = 10, height = 7)

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  species = preds[1,1] #species name
  years = c(2004:2013)

  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    between = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday <= dates[y,4]),]
  
    if (y == 1){ migpreds = between }
    else{ migpreds = rbind(migpreds, between) }
  }
  
  if (species == "Selasphoruscalliope") {by = 1}
  else {by = 0.5}
  
  ymax = max(c(migpreds$lat_se, migpreds$lon_se))
  latlon = ggplot(migpreds, aes(lon_se, lat_se)) + ggtitle(species) +
    xlab("estimated longitude standard error") + ylab("estimated latitude standard error") +
  geom_point(alpha = 0.5) + theme_classic() + facet_wrap(~year) +
  theme(text = element_text(size=12)) +
  scale_y_continuous(breaks = seq(0, ymax, by = by), limits = c(0, ymax)) +
  scale_x_continuous(breaks = seq(0, ymax, by = by), limits = c(0, ymax))
  multiplot(latlon, cols = 1)
}
dev.off()

#latlon = ggplot(migpreds2, aes(lon_se, lat_se)) + xlab("estimated longitude standard error") + 
#ylab("estimated latitude standard error") + geom_point(alpha=0.5, aes(col=cutoff)) + 
#theme_bw() + facet_grid(year~spname) + theme(text=element_text(size=12)) + guides(colour=FALSE)

#---------------------------------
#       plot standard error in predicted centroids across years, by julian date
#---------------------------------
pdf(file = paste(figpath, "/error_se-latlon_all_species.pdf", sep=""), width = 10, height = 4)

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = subset(dates, year > 2007)
  species = preds[1,1] #species name
  years = c(2008:2013)
  
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    between = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday <= dates[y,4]),]
    
    if (y == 1){ migpreds = between }
    else{ migpreds = rbind(migpreds, between) }
  }

  ymax = max(c(migpreds$lat_se, migpreds$lon_se)) 
  lat = ggplot(migpreds, aes(jday, lat_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
  geom_vline(xintercept = c(dates$spr_begin), col = "cadetblue") +
  geom_vline(xintercept = c(dates$fal_end), col = "orange") + ggtitle(species) +
  xlab("Julian Day") + ylab("latitude centroid standard error") +
  scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax))+
  theme(text = element_text(size=12))

  lon = ggplot(migpreds, aes(jday, lon_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
  geom_vline(xintercept = c(dates$spr_begin), col = "cadetblue") +
  geom_vline(xintercept = c(dates$fal_end), col = "orange") + ggtitle(species) +
  xlab("Julian Day") + ylab("longitude centroid standard error") +
  scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax)) +
  theme(text = element_text(size=12))

  multiplot(lat, lon, cols = 2)
}
dev.off()


#---------------------------------
#       plot predicted latitude and longitude as distribution across the years for all species
#---------------------------------
# save plots comparing daily lat and long and migration date across the years
pdf(file = paste(figpath, "/AllYears_lon-lat_all_species.pdf", sep=""), width = 8, height = 5)

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  preds = subset(preds, year > 2007)
  dates = subset(dates, year > 2007)
  species = preds[1,1] #species name
  
  yrlylat = ggplot(preds, aes(jday, lat, col=year)) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr_begin), col = "cadetblue") +
    geom_vline(xintercept = c(dates$spr_end), col = "olivedrab3") +
    geom_vline(xintercept = c(dates$fal_end), col = "orange") +
    scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
    theme(text = element_text(size=20)) + ggtitle(species)
  
  yrlylon = ggplot(preds, aes(jday, lon, col=year)) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr_begin), col = "cadetblue") +
    geom_vline(xintercept = c(dates$spr_end), col = "olivedrab3") +
    geom_vline(xintercept = c(dates$fal_end), col = "orange") +
    scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
    theme(text = element_text(size=20))

  multiplot(yrlylat, cols = 1)
  multiplot(yrlylon, cols = 1)
}
dev.off()

#----------------------------------------------------------
#       Makes Figure B3. Observations across the years (2004-2013)
#       For Appendix B
#----------------------------------------------------------

#NOTE!!: For this to work, make sure files and mfiles is in the same species-order!!
for (f in 1:length(files)) {
  
  #read in raw data
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = subset(dates, year > 2007)
  
  humdat = read.table(files[f], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  names(humdat) = c("SCI_NAME", "PRIMARY_COM_NAME","YEAR", "DAY", "TIME", "GROUP_ID", "PROTOCOL_ID",
                    "PROJ_ID", "DURATION_HRS", "EFFORT_DISTANCE_KM", "EFFORT_AREA_HA", "NUM_OBSERVERS",
                    "LATITUDE", "LONGITUDE", "SUB_ID", "POLYFID", "MONTH")
  
  humdat$MONTH = factor(humdat$MONTH, levels=c(1:12), ordered=TRUE)
  
  #grab species name for later identification
  species = humdat$SCI_NAME[1]
  species = gsub("\"", "", species, fixed=TRUE)
  years = c(2004:2013)
  
  for (y in 1:length(years)){
    dat = subset(humdat, YEAR==years[y])
    yrobs = data.frame(table(dat$DAY))
    yrobs$year = years[y]
    
    if(y == 1) { obs = yrobs}
    else{ obs = rbind(obs, yrobs)}
  }
    names(obs) = c("jday", "freq")
    obs$jday = as.numeric(obs$jday)
    obs$freq = as.numeric(obs$freq)
  ymax = max(obs$freq) + 25
  
  ggplot(obs, aes(jday, freq)) + geom_point(alpha=0.4) + xlab("Julian day of Year") + ylab("Number of Checklists") + 
    theme_classic() + theme(text = element_text(size=12)) +     
    geom_smooth(se=T, method='gam', formula=y~s(x, k=40), gamma=1.5) +
    scale_x_continuous(breaks = seq(0, 365, by = 50), limits = c(0,365)) + 
    scale_y_continuous(breaks = seq(0, ymax, by = 50)) +
    geom_vline(xintercept = dates[,1], col = "cadetblue", size = 0.5) +
    geom_vline(xintercept = dates[,4], col = "orange", size = 0.5) +
    geom_vline(xintercept = dates[,2], col = "olivedrab3", size = 0.5) 

  ggsave(file=paste(figpath, "/B3_occurrences", species, years[y], ".tiff", sep=""), dpi=600, height=3, width=4)
}


#-----------------------------------------------------------
#       Make synthetic panel figure for manuscript  - FIGURE 1
#       Make parts separately, and bring back together 
#         in Illustrator or Inkscape
#-----------------------------------------------------------
#reorder files to match manuscript ordering
files = c(files[4],files[1],files[2],files[3],files[5])
cfiles = c(cfiles[2],cfiles[1],cfiles[4],cfiles[3],cfiles[5])
mfiles = c(mfiles[2],mfiles[1],mfiles[4],mfiles[3],mfiles[5])

for (f in 1:length(files)) {
  
  #read in raw data
  humdat = read.table(files[f], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  names(humdat) = c("SCI_NAME", "PRIMARY_COM_NAME","YEAR", "DAY", "TIME", "GROUP_ID", "PROTOCOL_ID",
                    "PROJ_ID", "DURATION_HRS", "EFFORT_DISTANCE_KM", "EFFORT_AREA_HA", "NUM_OBSERVERS",
                    "LATITUDE", "LONGITUDE", "SUB_ID", "POLYFID", "MONTH")
  
  humdat$MONTH = factor(humdat$MONTH, levels=c(1:12), ordered=TRUE)
  
  #grab species name for later identification
  species = humdat$SCI_NAME[1]
  species = gsub("\"", "", species, fixed=TRUE)
   
  #plot the geographic number of checklists for the species, over all the years  
  #find the number of obs in each cell
  t = table(as.factor(humdat$POLYFID))
  
  #Merge the count data for the hexes with all the hexes in the map
  df = data.frame(POLYFID = names(t), count=as.numeric(t))
  df2 = data.frame(POLYFID = unique(hexgrid$POLYFID))
  df3 = merge(df2, df, all.x=TRUE)
  df3$species = species
  
  #combine all the datasets so we can get the colors on the same scale
  if (f == 1) {
    df3.1 = df3
    splist = species
  }
  
  else{
    df3.1 = rbind(df3.1, df3)
    splist = c(splist, species)
  }
}

# Set colors and legend scale by all species (e.g. some will be red, some will be shades of yellow)
  #matches colors with the number of observations
  yltord = colorRampPalette(brewer.pal(9, "YlOrRd"))(length(unique(df3.1$count)))
  cols = data.frame(id=c(NA,sort(unique(df3.1$count))), cols=yltord, stringsAsFactors=FALSE)
  df4 = merge(df3.1, cols, by.x="count", by.y="id")
  #set scale for legend
  vls = sort(unique(round(cols$id/500)*500))
  vls[1] = 1
  cols2 = rev(heat.colors(length(vls)))

#separate back out by species for plotting
for (s in 1:length(splist)){
  df4.1 = subset(df4, species == splist[s])
  df5 = merge(hexgrid, df4.1, by.x="POLYFID", by.y="POLYFID", all.x=TRUE)
  df5$cols = ifelse(is.na(df5$count), "white", df5$cols)  #hexes with no counts are white
  df5 = df5[order(df5$POLYFID),]

  #make a map with hexes colored by the number of times the species was observed in a given hex
  setEPS()
  postscript(file = paste(figpath, "/fig1_col1_", splist[s], ".eps", sep=""), width = 3, height = 3)
  plot(NA, NA, xlim = c(-140,-60), ylim=c(15,55),xlab="", ylab="", axes=FALSE)
    plot(hexgrid, col=df5$cols, border = "white", lwd = 0.25, las=1, add=TRUE)
    mtext(side=2,line=2, splist[s])
    map("worldHires", c("usa", "canada", "mexico"), add=TRUE, cex = 0.25)
  dev.off()
  
  #plot the legend separately
  setEPS()
  postscript(file = paste(figpath, "/fig1_col1_legend.eps", sep=""), width = 3, height = 3)
  plot(NA, NA, xlim = c(0,5), ylim=c(0,5), axes=FALSE, xlab = "", ylab="")
  legend.col(col = cols[,2], lev = sort(unique(df4$count)))
  dev.off()
}

for (f in 1:length(files)) {
  #read in the centroid and migration data (by western flyway), use data for years 2008-2013
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  preds_sub = preds[which(preds$year > 2007),]
  dates_sub = dates[which(dates$year > 2007),]
  species = preds[1,1]
  years = c(2008:2013)
  
# plot the migration routes in the next figure
preds_sub$month = as.factor(preds_sub$month)
cols3 = data.frame(id=c(sort(unique(preds_sub$month))), cols=tim.colors(length(unique(preds_sub$month))), stringsAsFactors=FALSE)
preds_sub = merge(preds_sub, cols3, by.x="month", by.y="id")
#set color scale
vls = sort(unique(round(cols3$id)))
vls[1] = 1
cols4 = tim.colors(length(vls))

setEPS()
postscript(file = paste(figpath, "/fig1_col2_", species, ".eps", sep=""), width = 3, height = 3)
plot(preds_sub$lon, preds_sub$lat, col=preds_sub$cols, pch=19, cex=0.25, ylab="", xlab="", xlim = c(-140, -60), ylim = c(15, 55),
     cex.lab = 1, cex.axis=1, axes=FALSE)
  map("worldHires", c("usa", "canada", "mexico"), add=TRUE)
  points(preds_sub$lon, preds_sub$lat, col=preds_sub$cols, pch=19, cex=0.25)
dev.off()

#plot legend separately  
setEPS()
postscript(file = paste(figpath, "/fig1_col2_legend.eps", sep=""), width = 3, height = 3)
plot(NA, NA, xlim = c(0,5), ylim=c(0,5), axes=FALSE, xlab = "", ylab="")
  legend("bottomleft", legend=vls, pch=22, pt.bg=cols4, pt.cex=0.75, cex=0.75, bty="n",
       col="black", title="month", x.intersp=1, y.intersp=0.5)  
dev.off()

# plot the latitudinal patterns with estimated dates of migration
setEPS()
postscript(file = paste(figpath, "/fig1_col4_", species, ".eps", sep=""), width = 3, height = 3)
latmin = ((min(preds_sub$lat)%/%5+1)*5)-5
latmax = ((max(preds_sub$lat)%/%5+1)*5) + 5 

par(bty="l")
plot(preds_sub$jday, preds_sub$lat, pch = 19, xlab = "", ylab = "", 
     cex.axis = 1, col = "grey20", cex = 0.25)
  abline(v=dates_sub$spr_begin, col="cadetblue")
  abline(v=dates_sub$spr_end, col = "olivedrab3")
  abline(v=dates_sub$fal_end, col = "orange")
dev.off()

}




#---------------------------------
#       plot predicted latitude with 95% confidence intervals
#---------------------------------
# save plots comparing daily lat and long and migration date across the years

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  preds = subset(preds, year > 2007)
  dates = subset(dates, year > 2007)
  species = preds[1,1] #species name
  
  #calculate 95% CI
  preds$lat_ucl = preds$lat + 1.96 * preds$lat_se
  preds$lat_lcl = preds$lat - 1.96 * preds$lat_se
  preds$lon_ucl = preds$lon + 1.96 * preds$lon_se
  preds$lon_lcl = preds$lon - 1.96 * preds$lon_se
  
  years = c(2008:2013)
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    between = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday <= dates[y,4]),]
    spring = preds[which(preds$year == years[y] & preds$jday >= dates[y,1] & preds$jday < dates[y,2]),]
    fall = preds[which(preds$year == years[y] & preds$jday > dates[y,2] & preds$jday <= dates[y,4]),]
    
    if (y == 1){
      migpreds = between
      pred_spr = spring
      pred_fal = fall
    }
    else{
      migpreds = rbind(migpreds, between)
      pred_spr = rbind(pred_spr, spring)
      pred_fal = rbind(pred_fal, fall)
    }
  }
     
  # plot spring and fall separately, color-coded to compare overlap
  seasons = ggplot(pred_spr, aes(lon, lat), col = year, group=year, order.by=jday) + theme_classic() + 
    geom_rect(data=pred_spr, aes(xmin=lon_lcl,ymin=lat_lcl,xmax=lon_ucl,ymax=lat_ucl, group=year), fill="cadetblue", alpha=0.02) +
    geom_path(aes(group=year)) + 
    geom_rect(data=pred_fal, aes(xmin=lon_lcl,ymin=lat_lcl,xmax=lon_ucl,ymax=lat_ucl, group=year), fill="orange", alpha=0.02) +
    geom_path(data=pred_fal, aes(lon, lat, group=year)) + 
    theme(text = element_text(size=12)) + xlab("longitude") + ylab("latitude")
  
  pdf(file = paste(figpath, "/fig1_col3_", species, ".pdf", sep=""), width = 3, height = 3)
  multiplot(seasons, cols = 1)
  dev.off()
}



#---------------------------------
#       Figure B1 for appendix, total eBird effort 2004-2013
#---------------------------------

#find the number of obs in each cell
t = table(as.factor(effort$POLYFID))

#Merge the count data for the hexes with all the hexes in the map
df = data.frame(POLYFID = names(t), count=as.numeric(t))
df2 = data.frame(POLYFID = unique(hexgrid$POLYFID))
df3 = merge(df2, df, all.x=TRUE)

#matches colors with the number of observations
#   # Set colors and legend scale by all species (e.g. some will be red, some will be shades of yellow)
yltord = colorRampPalette(brewer.pal(9, "YlOrRd"))(length(unique(df3$count)))
cols = data.frame(id=c(NA,sort(unique(df3$count))), cols=yltord, stringsAsFactors=FALSE)
df4 = merge(df3, cols, by.x="count", by.y="id")
df5 = merge(hexgrid, df4, by.x="POLYFID", by.y="POLYFID", all.x=TRUE)

#   #set scale for legend
vls = sort(unique(round(cols$id/500)*500))
vls[1] = 1
cols2 = rev(heat.colors(length(vls)))

#hexes with no counts are white
df5$cols = ifelse(is.na(df5$count), "white", df5$cols)
df5 = df5[order(df5$POLYFID),]

#make a map with hexes colored by the number of times the species was observed in a given hex
setEPS()
postscript(file = paste(figpath, "/FigureB1.eps", sep=""), width = 5.5, height = 4)
plot(NA, NA, xlim = c(-170,-50), ylim=c(15,75),xlab="", ylab="", axes=FALSE)
plot(hexgrid, col=df5$cols, border = "white", lwd = 0.25, las=1, add=TRUE)
map("worldHires", c("usa", "canada", "mexico"), add=TRUE, cex = 0.25)
dev.off()

#plot the legend separately
setEPS()
postscript(file = paste(figpath, "/figB1_col1_legend.eps", sep=""), width = 3, height = 3)
plot(NA, NA, xlim = c(0,5), ylim=c(0,5), axes=FALSE, xlab = "", ylab="")
legend.col(col = cols[,2], lev = sort(unique(df4$count)))
dev.off()


#--------------------------------------------------
#   wrt question about number of obs in eastern US
#--------------------------------------------------
for (f in 1:length(files)){
  humdat = read.table(files[f], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")

  names(humdat) = c("SCI_NAME", "PRIMARY_COM_NAME","YEAR", "DAY", "TIME", "GROUP_ID", "PROTOCOL_ID",
                  "PROJ_ID", "DURATION_HRS", "EFFORT_DISTANCE_KM", "EFFORT_AREA_HA", "NUM_OBSERVERS",
                  "LATITUDE", "LONGITUDE", "SUB_ID", "POLYFID", "MONTH")

  humdat$MONTH = factor(humdat$MONTH, levels=c(1:12), ordered=TRUE)

  #grab species name for later identification
  species = humdat$SCI_NAME[1]

  #subset eastern observations
  east = subset(humdat, LONGITUDE > -103)
  eastwinter = subset(east, MONTH %in% c(11, 12, 1, 2, 3))
  plot(table(east$YEAR), main=species)
  plot(table(eastwinter$YEAR), main=species)
}
  

#----------------------------------------------
#       Reviewer comment: detecting trends
#----------------------------------------------

#aggregate all the dates files
for (f in 1:length(mfiles)){
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  if (f == 1) { d = dates }
  else { d = rbind(d, dates) }
}
d = d[,-2]
d = subset(d, year > 2007)
names(d) = c("spring_begin", "peak_latitude", "autumn_end", "species", "year")

#plot the data
spring = ggplot(d, aes(year, spring_begin, group=species)) + geom_line(aes(linetype=species)) + 
  theme_bw() + ylab("begin spring migration")
breed = ggplot(d, aes(year, peak_latitude, group=species)) + geom_line(aes(linetype=species)) + 
  theme_bw() + ylab("reach peak latitude")
autumn = ggplot(d, aes(year, autumn_end, group=species)) + geom_line(aes(linetype=species)) + 
   theme_bw() + ylab("end autumn migration")

pdf(file = paste(figpath, "/migration_trends.pdf", sep=""), width = 6, height = 7)
  multiplot(spring, breed, autumn, cols=1)
dev.off()


library(plyr)
# Break up d by species, then fit the specified model to each piece and return a list
# SPRING
models <- dlply(d, "species", function(df) 
  lm(spring_begin ~ year, data = df))
# Apply coef to each model and return a data frame
ldply(models, coef)
# Print the summary of each model
l_ply(models, anova, .print = TRUE)
l_ply(models, summary, .print = TRUE)

# PEAK LAT
models <- dlply(d, "species", function(df) 
  lm(peak_latitude ~ year, data = df))
# Apply coef to each model and return a data frame
ldply(models, coef)
# Print the summary of each model
l_ply(models, anova, .print = TRUE)

# AUTUMN
models <- dlply(d, "species", function(df) 
  lm(autumn_end ~ year, data = df))
# Apply coef to each model and return a data frame
ldply(models, coef)
# Print the summary of each model
l_ply(models, anova, .print = TRUE)
