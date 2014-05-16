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
    
      #FIXME: What to do with RT - use data EAST of 103 meridian? Write second loop.
      
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
                                  #FIXME: Also analyze data using flyways
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

setwd(wd)

# read in eBird data for migration speed
rfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = "west_speed", full.names=TRUE)

for (f in 1:length(rfiles)){
  rate = read.table(rfiles[f], header=TRUE, sep=" ", fill=TRUE, comment.char="")
  species = rate$species[1]
  print(species)
  
  dirpath = paste(figpath, "/", species, sep="")
  
  rate_sub = rate[which(rate$year > 2007),]
  
  #plot the variance in estimated migration speed for all and for recent years
  pdf(file = paste(figpath, "/speed_", species, ".pdf", sep=""), width = 10, height = 4)
  
  r = melt(rate[,c(1,2,4)], id.vars = "year")
  names(r) = c("year", "season", "rate")
  bxp_speed = ggplot(r, aes(season, rate, fill=season)) + geom_boxplot() + theme_classic() + 
    scale_fill_manual(values=c("cadetblue", "orange")) + ylab("km/day")
  
  r_sub = r[which(r$year>2007),]
  bxp_speed_sub = ggplot(r_sub, aes(season, rate, fill=season)) + geom_boxplot() + theme_classic() + 
    scale_fill_manual(values=c("cadetblue", "orange")) + ylab("km/day")
  
  multiplot(bxp_speed, bxp_speed_sub, cols = 2)
  dev.off()
  
  #print the mean and standard deviation of speed
  print(paste("spring sd:", sd(rate_sub$spr)))
  print(paste("spring mean:", mean(rate_sub$spr)))
  print(paste("fall sd:", sd(rate_sub$fal)))
  print(paste("fall mean:", mean(rate_sub$fal)))
}


mfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = c("west_migration.*txt"), full.names=TRUE)

for (f in 1:length(mfiles)){
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates_sub = dates[which(dates$year > 2007),]
  species = dates$species[1]
  print(species)
  
  #plot the variance in estimated migration begin and end for all and for recent years
  pdf(file = paste(figpath, "/migdates_", species, ".pdf", sep=""), width = 10, height = 4)
  
  d = melt(dates[,c(1,2,4,6)], id.vars = "year")
  names(d) = c("year", "season", "date")
  d_sub = d[which(d$year>2007),]  
  
  bxp_date = ggplot(d, aes(season, date, fill=season)) + geom_boxplot() + theme_classic() + 
    scale_fill_manual(values=c("cadetblue", "olivedrab3", "orange"))
  
  bxp_date_sub = ggplot(d_sub, aes(season, date, fill=season)) + geom_boxplot() + theme_classic() + 
    scale_fill_manual(values=c("cadetblue", "olivedrab3", "orange"))
  
  multiplot(bxp_date, bxp_date_sub, cols = 2)
  dev.off()
  
  print(paste("spring begin median:", median(dates_sub$spr_begin)))
  print(paste("spring begin sd:", sd(dates_sub$spr_begin)))
  print(paste("lat peak median:", median(dates_sub$spr_end)))
  print(paste("lat peak sd:", sd(dates_sub$spr_end)))
  print(paste("fall end median:", median(dates_sub$fal_end)))
  print(paste("fall end sd:", sd(dates_sub$fal_end)))
  print("")
}


#read in all pred data, then re-analyze based on se results. 
#compare 2008-2013, test for impact of 2004-2007 years on overall distribution
#linear model on spring vs. fall in each year

# read in predicted data
cfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = "west_centroids.*.txt", full.names=TRUE)

for (f in 1:length(cfiles)){
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  species = preds[1,1]
  years = c(2004:2013)
  preds_sub = preds[which(preds$year > 2007),]
  dates_sub = dates[which(dates$year > 2007),]
  
  #set species-specific directory path for figures
  dirpath = paste(figpath, "/", species, sep="")
  
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    spr_begin = dates[y,1]
    spr_end = dates[y,2]
    fal_begin = dates[y,2]
    fal_end = dates[y,4]
    
    between = preds[which(preds$year == years[y] & preds$jday >= spr_begin & preds$jday <= fal_end),]
    #split spring and fall on estimated breeding bounds
    spring = preds[which(preds$year == years[y] & preds$jday >= spr_begin & preds$jday < spr_end),]
    fall = preds[which(preds$year == years[y] & preds$jday > fal_begin & preds$jday <= fal_end),]
    
    #get slope & r2 of spring and fall latitudinal migration
    #spr_lm = LinearMigration(spring, years[y])
    #fal_lm = LinearMigration(fall, years[y])
    
    if (y == 1){
      migpreds = between
      pred_spr = spring
      pred_fal = fall
      #lm_spr = spr_lm
      #lm_fal = fal_lm
    }
    else{
      migpreds = rbind(migpreds, between)
      pred_spr = rbind(pred_spr, spring)
      pred_fal = rbind(pred_fal, fall)
      #lm_spr = rbind(lm_spr, spr_lm)
      #lm_fal = rbind(lm_fal, fal_lm)
    }
  }
  
  #compare location across the years using mean and sd for ALL years
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
  
  #compare location across the years using mean and sd for 2008:2013
  migpreds_sub = migpreds[which(migpreds$year > 2007),]
  jdays = sort(unique(migpreds_sub$jday))
  patherr_sub = data.frame("jday"=1,"meanlat"=1, "sdlat"=1, "meanlon"=1, "sdlon"=1)
  outcount = 1
  for (j in min(jdays):max(jdays)){
    tmp = migpreds_sub[which(migpreds_sub$jday == j),]
    meanlat = mean(tmp$lat)
    sdlat = sd(tmp$lat)
    meanlon = mean(tmp$lon)
    sdlon = sd(tmp$lon)
    patherr_sub[outcount,] = c(j, meanlat, sdlat, meanlon, sdlon)
    outcount = outcount + 1
  }
  
  #----------- model to test variance in lat and lon across years, with year as a random effect

  lon.gam <- gamm4(lon ~ s(jday, k=10), random = ~(1|year), data = preds_sub, gamma = 1.5)
  print (paste("R2 for spring Longitude is:", round(summary(lon.gam$gam)$r.sq,4)))
  lat.gam <- gamm4(lat ~ s(jday, k=10), random = ~(1|year), data=preds_sub, gamma = 1.5)
  print (paste("R2 for spring Latitude is:", round(summary(lat.gam$gam)$r.sq,4)))
  
  #what is the variation in daily location across years using 2008: 2013 data?
  pred_spr_sub = pred_spr[which(pred_spr$year > 2007),]
  pred_fal_sub = pred_fal[which(pred_fal$year > 2007),]
  
  lon.gam.spr <- gamm4(lon ~ s(jday, k=10), random = ~(1|year), data=pred_spr_sub, gamma = 1.5)
    print (paste("R2 for spring Longitude is:", round(summary(lon.gam.spr$gam)$r.sq,4)))
  lat.gam.spr <- gamm4(lat ~ s(jday, k=10), random = ~(1|year), data=pred_spr_sub, gamma = 1.5)
   print (paste("R2 for spring Latitude is:", round(summary(lat.gam.spr$gam)$r.sq,4)))
  
  lon.gam.fal <- gamm4(lon ~ s(jday, k=10), random = ~(1|year), data=pred_fal_sub, gamma = 1.5)
  print (paste("R2 for fall Longitude is:", round(summary(lon.gam.fal$gam)$r.sq,4)))
  lat.gam.fal <- gamm4(lat ~ s(jday, k=10), random = ~(1|year), data=pred_fal_sub, gamma = 1.5)
  print (paste("R2 for fall Latitude is:", round(summary(lat.gam.fal$gam)$r.sq,4)))
}

  #----------- plot the data ------------
  
  # plot all the routes onto a single map
  #save a plot of the species migration mapped onto an elevation raster
  pdf(file = paste(dirpath, "/elev-route_summary_", species, ".pdf", sep=""), width = 7, height = 4.5)
  AllMigration(preds_sub, elev, myext, species)
  dev.off() 

  # compare the slope for lat and lon change in spring vs. fall
  pdf(file = paste(dirpath, "/slope_lon-lat", species, ".pdf", sep=""), width = 10, height = 4)
  
  lat_compare = ggplot(lm_spr, aes(year, abs(lat_slope))) +  geom_line() +
    geom_point(col="cadetblue", aes(size=lat_r2)) + geom_line(data=lm_fal,aes(year, abs(lat_slope))) + 
    geom_point(data=lm_fal, aes(year, abs(lat_slope), size=lat_r2), col="orange") + 
    theme_classic() + ylab("slope") + ggtitle("Latitude spring vs. fall")
  
  lon_compare = ggplot(lm_spr, aes(year, abs(lon_slope))) +  geom_line() +
    geom_point(col="cadetblue", aes(size=lon_r2)) + geom_line(data=lm_fal,aes(year, abs(lon_slope))) + 
    geom_point(data=lm_fal, aes(year, abs(lon_slope), size=lon_r2), col="orange") + 
    theme_classic() + ylab("slope") + ggtitle("Longitude spring vs. fall")
  
  multiplot(lat_compare, lon_compare, cols = 2)
  dev.off()
  
  # save plots comparing daily lat and long and migration date across the years
  pdf(file = paste(dirpath, "/AllYears_lon-lat", species, ".pdf", sep=""), width = 8, height = 10)
  
  yrlylon = ggplot(preds_sub, aes(jday, lon, col=year)) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates_sub$spr_begin), col = "cadetblue") +
    geom_vline(xintercept = c(dates_sub$spr_end), col = "olivedrab3") +
    geom_vline(xintercept = c(dates_sub$fal_end), col = "orange") +
    scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
    theme(text = element_text(size=20)) + ggtitle(species)
  
  yrlylat = ggplot(preds_sub, aes(jday, lat, col=year)) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates_sub$spr_begin), col = "cadetblue") +
    geom_vline(xintercept = c(dates_sub$spr_end), col = "olivedrab3") +
    geom_vline(xintercept = c(dates_sub$fal_end), col = "orange") +
    scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
    theme(text = element_text(size=20)) + ggtitle(species)
  
  multiplot(yrlylon, yrlylat, cols = 1)
  dev.off()
  
# save plots comparing daily LAT ONLY and migration date across the years
pdf(file = paste(dirpath, "/West_Years_lat", species, ".pdf", sep=""), width = 8, height = 5)

latmin = ((min(preds_sub$lat)%/%5+1)*5)-5
latmax = ((max(preds_sub$lat)%/%5+1)*5) + 5

yrlylat = ggplot(preds_sub, aes(jday, lat, col=year)) + geom_point(size=1) + theme_classic() +
  geom_vline(xintercept = c(dates_sub$spr_begin), col = "cadetblue") +
  geom_vline(xintercept = c(dates_sub$spr_end), col = "olivedrab3") +
  geom_vline(xintercept = c(dates_sub$fal_end), col = "orange") +
  scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
  scale_y_continuous(breaks = seq(latmin, latmax, by = 5)) +
  theme(text = element_text(size=12)) #+ ggtitle(species)

yrlylat

dev.off()

  #compare standard error in predicted centroids across years
  pdf(file = paste(dirpath, "/Error_selon-lat", species, ".pdf", sep=""), width = 10, height = 4)
  
  ymax = max(c(migpreds$lat_se, migpreds$lon_se)) 
  lat = ggplot(migpreds, aes(jday, lat_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr_begin), col = "cadetblue") +
    geom_vline(xintercept = c(dates$fal_end), col = "orange") + ggtitle(paste(species, "Latitude")) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax))
  
  lon = ggplot(migpreds, aes(jday, lon_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr_begin), col = "cadetblue") +
    geom_vline(xintercept = c(dates$fal_end), col = "orange") + ggtitle(paste(species, "Longitude")) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax))
  
  multiplot(lat, lon, cols = 2)
  dev.off()
  
  # plot the relationship between lon se and lat se for each year
  pdf(file = paste(dirpath, "/Error_corlon-lat", species, ".pdf", sep=""), width = 10, height = 7)
  
  ymax = max(c(migpreds$lat_se, migpreds$lon_se))
  latlon = ggplot(migpreds, aes(lon_se, lat_se)) + ggtitle(species) +
    geom_point(alpha = 0.5) + theme_classic() + facet_wrap(~year) +
    theme(text = element_text(size=12)) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) +
    scale_x_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax))
  multiplot(latlon, cols = 1)
  dev.off()
  
  # plot the standard deviation in daily lat and lon across the 10 years, and across 6 most recent years
  pdf(file = paste(dirpath, "/ErrorinDailyLocs", species, ".pdf", sep=""), width = 10, height = 8)
  
  ymax = max(c(patherr$sdlat, patherr$sdlon),na.rm=TRUE) 
  sdlocs = ggplot(patherr, aes(jday, sdlat)) + geom_point(size=1) + theme_classic() + ggtitle(paste(species, "sd in daily locs, 04-13")) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) + 
    scale_x_continuous(breaks = seq(0, 366, by = 25), limits = c(0, 366)) + 
    geom_point(aes(jday, sdlon), col = "indianred", size=1) + ylab("stdev daily lat (black) and lon (red)")
  sdlocs_sub = ggplot(patherr_sub, aes(jday, sdlat)) + geom_point(size=1) + theme_classic() + ggtitle(paste(species, "sd in daily locs, 08-13")) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) + 
    scale_x_continuous(breaks = seq(0, 366, by = 25), limits = c(0, 366)) + 
    geom_point(aes(jday, sdlon), col = "indianred", size=1) + ylab("stdev daily lat (black) and lon (red)")
  
  lonlatday = ggplot(patherr, aes(sdlon, sdlat, col=jday)) + geom_point(size=1) + theme_classic() +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax))  +
    scale_x_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) 
  lonlatday_sub = ggplot(patherr_sub, aes(sdlon, sdlat, col=jday)) + geom_point(size=1) + theme_classic() +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax))  +
    scale_x_continuous(breaks = seq(0, ymax, by = 0.5), limits = c(0, ymax)) 
  
  multiplot(sdlocs, lonlatday, sdlocs_sub, lonlatday_sub, cols = 2)
  dev.off()
  
  # save plots comparing spring vs fall migration routes across the years
  pdf(file = paste(dirpath, "/AllYears_sprVSfal", species, ".pdf", sep=""), width = 10, height = 4)
  
  sprplot = ggplot(pred_spr, aes(lon, lat, col=year, group=year)) + geom_point(size=1) + theme_classic() +
    theme(text = element_text(size=20)) + ggtitle(paste("spring-", species))
  
  falplot = ggplot(pred_fal, aes(lon, lat, col=year, group=year)) + geom_point(size=1) + theme_classic() +
    theme(text = element_text(size=20)) + ggtitle(paste("fall-", species))  
  
  multiplot(sprplot, falplot, cols = 2)
  dev.off()
}




#--------------------------------- Make synthetic panel figure for paper draft --------------------------
#get all the pathnames for the 5 species data, ordered by species so you can use the same iterator
cfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = "west_centroids.*.txt", full.names=TRUE)
mfiles = list.files(path = paste(getwd(), "/output_data/", sep=""), pattern = c("west_migration.*txt"), full.names=TRUE)
rtfiles = list.files(path = paste(getwd(), "/output_data/", sep = ""), 
                     pattern = "Archilochuscolubris*.txt",full.names=TRUE)
files = list.files(pattern = "*.txt")
files = files[c(1,3,2,5,4)]
cfiles = c(cfiles, rtfiles[1]) #centroid data
mfiles = c(mfiles, rtfiles[2]) #migration data


#Open pdf plotting window
setEPS()
postscript(file = paste(wd, "/Panel_figure.eps", sep=""), width = 7.5, height = 10)
par(mfrow=c(5,3), mai=c(0.2,0.5,0.2,0.5), oma = c(0, 2, 2, 2))

for (f in 1:length(files)) {
  
  #read in raw data
  humdat = read.table(files[f], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  names(humdat) = c("SCI_NAME", "PRIMARY_COM_NAME","YEAR", "DAY", "TIME", "GROUP_ID", "PROTOCOL_ID",
                    "PROJ_ID", "DURATION_HRS", "EFFORT_DISTANCE_KM", "EFFORT_AREA_HA", "NUM_OBSERVERS",
                    "LATITUDE", "LONGITUDE", "SUB_ID", "POLYFID", "MONTH")
  
  humdat$MONTH = factor(humdat$MONTH, levels=c(1:12), ordered=TRUE)
  
  #read in the centroid and migration data (by western flyway), use data for years 2008-2013
  preds = read.table(cfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  preds_sub = preds[which(preds$year > 2007),]
  dates_sub = dates[which(dates$year > 2007),]
  species = preds[1,1]
  years = c(2008:2013)
   
  #plot the geographic number of checklists for the species, over all the years  
  #find the number of obs in each cell
  t = table(as.factor(humdat$POLYFID))
  
  #Merge the count data for the hexes with all the hexes in the map
  df = data.frame(POLYFID = names(t), count=as.numeric(t))
  df2 = data.frame(POLYFID = unique(hexgrid$POLYFID))
  df3 = merge(df2, df, all.x=TRUE)
  
  #matches colors with the number of observations
  cols = data.frame(id=c(NA,sort(unique(df3$count))), cols=tim.colors(length(unique(df3$count))), stringsAsFactors=FALSE)
  df4 = merge(df3, cols, by.x="count", by.y="id")
  df5 = merge(hexgrid, df4, by.x="POLYFID", by.y="POLYFID", all.x=TRUE)
  df5$cols = ifelse(is.na(df5$count), "white", df5$cols)  #hexes with no counts are white
  df5 = df5[order(df5$POLYFID),]
  
  #set scale for legend
  vls = sort(unique(round(cols$id/500)*500))
  vls[1] = 1
  cols2 = tim.colors(length(vls))
  
  #make a map with hexes colored by the number of times the species was observed in a given hex
  plot(hexgrid, col=df5$cols, border="gray50", lwd=0.25, xlim=c(-170,-50), ylim=c(15,75), las=1)
  axis(side=1)
  axis(side=2, las=1)
  mtext(side=2,line=2, species)
  box()
  ## legend
  #legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1, cex=0.75, bty="",
  #       col="black", title="Number of checklists", x.intersp=1, y.intersp=0.5, bg="white")
  map("worldHires", c("usa", "canada", "mexico"), add=TRUE, cex = 0.5)

  
# plot the migration routes in the next figure
myext <- c(-175, -50, 15, 75) #extent for maps
preds_sub$month = as.factor(preds_sub$month)
cols3 = data.frame(id=c(sort(unique(preds_sub$month))), cols=tim.colors(length(unique(preds_sub$month))), stringsAsFactors=FALSE)
preds_sub = merge(preds_sub, cols3, by.x="month", by.y="id")
#set color scale
vls = sort(unique(round(cols3$id)))
vls[1] = 1
cols4 = tim.colors(length(vls))

plot(elev, ext=myext, ylab="", xlab="", xlim = c(-175, -50), ylim = c(15, 75),
     cex.lab = 1, cex.axis=1, col=gray(200:0/256))
points(preds_sub$lon, preds_sub$lat, col=preds_sub$cols, pch=19, cex=0.25)
#map("worldHires", c("usa", "canada", "mexico"), add=TRUE)
#legend("bottomleft", legend=vls, pch=22, pt.bg=cols4, pt.cex=0.75, cex=0.75, bty="n",
#       col="black", title="", x.intersp=0.5, y.intersp=0.75)  
  
# plot the latitudinal patterns with estimated dates of migration
latmin = ((min(preds_sub$lat)%/%5+1)*5)-5
latmax = ((max(preds_sub$lat)%/%5+1)*5) + 5 

plot(preds_sub$jday, preds_sub$lat, pch = 19, xlab = "", ylab = "", 
     cex.axis = 1, col = "grey20", cex = 0.5)
  abline(v=dates_sub$spr_begin, col="cadetblue")
  abline(v=dates_sub$spr_end, col = "olivedrab3")
  abline(v=dates_sub$fal_end, col = "orange")

}

dev.off()
