#Code for eBird migration project
# (c) 2013 -2014 Sarah Supp 

library(ggmap)
library(maptools)
library(sp)
library(raster)

#set working directory
main = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration"
figpath = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/Figures/"
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data"
setwd(wd)


#---------------------------------------------------------------------------------------
#              predict migration paths, dates, speed, and error across years
#---------------------------------------------------------------------------------------

# read in summary of effort data (Number of eBird checklists submitted per day per year)
effort = read.table("cell_effort_new.txt", header=TRUE, as.is=TRUE)

# read in the north america equal area hex grid map (F.A.L.)
#hexgrid = readShapePoly("/Volumes/Elements/eBird/terr_4h6/nw_vector_grid.shp") #quad map
hexgrid = readShapePoly("C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/terr_4h6/terr_4h6.shp") #hex map
  # crop to just North America, where the migratory species occur
  hexgrid = hexgrid[which(hexgrid$LATITUDE > 10 & hexgrid$LATITUDE <80 & 
                            hexgrid$LONGITUDE > -178 & hexgrid$LONGITUDE < -50),]

# make a North America base map
noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")

# read in eBird data
files = list.files(pattern = "*2013.txt")
files = files[c(3,5,6,8,9)]

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
  source("C:/Users/sarah/Documents/GitHub/hb-migration/migration-fxns.r")

  humdat = read.table(files[f], header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="") #quote="/"'"

  #keep only the columns that we need
  keepcols = c("COMMON.NAME", "SCIENTIFIC.NAME", "OBSERVATION.COUNT", "AGE.SEX", "COUNTRY",
               "COUNTRY_CODE", "STATE_PROVINCE", "COUNTY", "LATITUDE", "LONGITUDE",
               "OBSERVATION.DATE", "TIME.OBSERVATIONS.STARTED", "PROTOCOL.TYPE", "PROJECT.CODE",
               "DURATION.MINUTES", "EFFORT.DISTANCE.KM", "EFFORT.AREA.HA", "NUMBER.OBSERVERS")
  humdat = humdat[,which(names(humdat) %in% keepcols)]

  species = humdat$SCIENTIFIC.NAME[1]
  years = c(2004:2013)
  
  # separate year, month, day, and calculate julian date for later analysis and plotting
  year <-sapply(humdat$OBSERVATION.DATE,function(x){
    as.numeric(substring(x, 1, 4))
  })
  
  month <- sapply(humdat$OBSERVATION.DATE,function(x){
    as.numeric(substring(x, 6, 7))
  }) 
  
  day <- sapply(humdat$OBSERVATION.DATE,function(x){
    as.numeric(substring(x, 9, 10))
  })
  
  julian = sapply(humdat$OBSERVATION.DATE, function(x){
    julian(as.numeric(substring(x, 6, 7)), as.numeric(substring(x, 9, 10)), as.numeric(substring(x, 1, 4)), 
           origin. = c(1, 1, as.numeric(substring(x, 1, 4)))) + 1
  })
  
  humdat = cbind(humdat, year, month, day, julian)
    
  #start a new directory
  dirpath = paste(figpath, species, sep="")
 #   dir.create(dirpath, showWarnings = TRUE, recursive = FALSE)
  
  #show how many records there are for the species across the years, write to txt file
  yeartable = PlotRecords(humdat$year, species)
  write.table(yeartable, file = paste(dirpath, "/", species,".txt",sep=""), row.names=FALSE)
    
  #save a figure of the geographic number of checklists for the species, over all the years
  count = PlotChecklistMap(humdat, hexgrid, dirpath)
  
  for (y in 1:length(years)){
    yrdat = humdat[which(humdat$year == years[y]),]
    yreffort = effort[which(effort$YEAR == years[y]),]
    
    #plot frequency of sightings per month
    monthtable = PlotRecords(yrdat$month, species)
    
    #get daily weighted mean location
    altmeandat = AlternateMeanLocs(yrdat,species,hexgrid,yreffort)
    
    #grab dates for migration
    migration = GetMigrationDates(altmeandat)
    
    #use GAM model to predict daily location along a smoothing line
    preds = EstimateDailyLocs(altmeandat)
    
    #get Great Circle distances traveled each day between predicted daily locations
    dist = DailyTravel(preds, 4, 5, species, years[y], migration)
    ggsave(filename = paste(dirpath, "/", "distance", years[y], species,".jpeg",sep=""))
    
    #estimate migration speed for spring and fall
    speed = MigrationSpeed(dist, migration)
    
    #plot smoothed migration trajectory for the species and year
    mig_path = PlotMigrationPath(preds, noam, species, years[y])
    ggsave(mig_path, file=paste(dirpath, "/", "migration", species, years[y], ".pdf", sep=""))
    
    #plot occurrences with lines showing beginning and end of migration
    PlotOccurrences(altmeandat, species, migration[[1]], migration[[2]])
    ggsave(file=paste(dirpath, "/", "occurrences", species, years[y], ".pdf", sep=""))
    
    #add year to preds, so we can save it to compare across years
    preds$year = years[y]
    
    if (y == 1){
      pred_data = preds
      migdates = data.frame("spr" = migration[1], "fal" = migration[2])
      migspeed = data.frame("spr" = speed[1], "fal" = speed[2])
    }
    else{
      pred_data = rbind(pred_data, preds)
      migdates = rbind(migdates, migration)
      migspeed = rbind(migspeed, speed)
    }
    
    if (y == length(years)){
      #write migration timing and speed data to file
      write.table(migdates, file = paste(getwd(), "/output_data/migration", species, ".txt", sep=""), 
                  append=FALSE, row.names=FALSE)
      write.table(migspeed, file = paste(getwd(), "/output_data/speed", species, ".txt", sep=""), 
                  append=FALSE, row.names=FALSE)
      write.table(pred_data, file = paste(getwd(), "/output_data/centroids", species, ".txt", sep=""), 
                  append=FALSE,row.names=FALSE)
    }
    
    rm(list=ls()[ls() %in% c("sitemap", "meanmap", "yrdat", "altmeandat", "migration", "preds", "dist", "mig_path")])   # clears the memory of the map and year-level data
  }
  rm(list=ls()[!ls() %in% c("f", "files", "noam", "hexgrid", "effort", "pred_data", "migdates", "wd", "main")])   # clears the memory of everything except the file list, iterator, and base map
}


#---------------------------------------------------------------------------------------
#                   compare migration rates and dates across species
#---------------------------------------------------------------------------------------

# read in eBird data
files = list.files(path = main, pattern = "speed", recursive=TRUE, full.names=TRUE)
#grab just the strongly migratory species
files = files[c(1,2,7,8)]

for (f in 1:length(files)){
  rate = read.table(files[f], header=FALSE, sep=" ", quote="", fill=TRUE, comment.char="") #quote="/"'"
  names(rate) = c("spring", "fall")
  year = c(2004:2013)
  rate = cbind(rate, year)
  
  avg_speed = ggplot(rate, aes(year, spring)) + geom_point(col = "cadetblue") + 
    geom_point(aes(year, fall), col = "orange") + 
    theme_classic() + ylab("migration speed (km/day)") + theme(text = element_text(size=20)) + 
    geom_hline(yintercept = mean(rate$spring), col = "cadetblue", position="identity") +
    geom_hline(yintercept = mean(rate$fall), col = "orange", position="identity")
    
  print(avg_speed)
  
  print(paste("spring sd:", sd(rate$spring)))
  print(paste("spring mean:", mean(rate$spring)))
  print(paste("fall sd:", sd(rate$fall)))
  print(paste("fall mean:", mean(rate$fall)))
}

mfiles = list.files(path = main, pattern = c("migration.*txt"), recursive = TRUE, full.names=TRUE)
#grab just the strongly migratory species
mfiles = mfiles[c(1,2,7,8)]

for (f in 1:length(mfiles)){
  dates = read.table(mfiles[f], header=FALSE, sep=" ", as.is=TRUE, quote="", fill=TRUE, comment.char="") #quote="/"'"
  dates=dates[-1,]
  names(dates) = c("spring", "fall")
  dates$spring = as.numeric(dates$spring)
  dates$fall = as.numeric(dates$fall)
  year = c(2004:2013)
  dates = cbind(dates, year)
  dates=dates[which(dates$year >2007),]
  
  avg_date = ggplot(dates, aes(year, spring)) + geom_point(col = "cadetblue") + 
    geom_point(aes(year, fall), col = "orange") + 
    theme_classic() + ylab("migration date (start/end)") + theme(text = element_text(size=20)) + 
    geom_hline(yintercept = mean(dates$spring), col = "cadetblue", position="identity") +
    geom_hline(yintercept = mean(dates$fall), col = "orange", position="identity")
  print(avg_date)
  
  print(paste("spring:", sd(dates$spring)))
  print(paste("fall:", sd(dates$fall)))
}


#TODO: read in all pred data, then re-analyze based on se results. 
#compare 2008-2013, test for impact of 2004-2007 years on overall distribution
#linear model on spring vs. fall in each year
#focus on 5 migratory species

# read in predicted data
files = list.files(path = main, pattern = "centroids.*.txt", recursive=TRUE, full.names=TRUE)
mfiles = list.files(path = getwd(), pattern = c("migration.*txt"), recursive = TRUE, full.names=TRUE)

for (f in 1:length(files)){
  preds = read.table(files[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  dates = read.table(mfiles[f], header=TRUE, sep=" ", as.is=TRUE, fill=TRUE, comment.char="")
  years = c(2004:2013)
  species = preds[1,1]
  
  #set species-specific directory path for figures
  dirpath = paste(figpath, species, sep="")
  
  #grab only the predicted daily centroids from between the migration dates
  for (y in 1:length(years)){
    begin = dates[y,1]
    end = dates[y,2]
    med = median(c(dates[y,1], dates[y,2]))
    
    between = preds[which(preds$year == years[y] & preds$jday >= begin & preds$jday <= end),]
    spring = preds[which(preds$year == years[y] & preds$jday >= begin & preds$jday < med),]
    fall = preds[which(preds$year == years[y] & preds$jday <= end & preds$jday > med),]
    
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
  
  #----------- plot the data ------------
  
  # save plots comparing daily lat and long and migration date across the years
  pdf(file = paste(dirpath, "/AllYears_lon-lat", species, ".pdf", sep=""), width = 8, height = 10)
  
  yrlylon = ggplot(preds, aes(jday, lon, col=year)) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr), col = "cadetblue") +
    geom_vline(xintercept = c(dates$fal), col = "orange") +
    scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
    theme(text = element_text(size=20)) + ggtitle(species)
  
  yrlylat = ggplot(preds, aes(jday, lat, col=year)) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr), col = "cadetblue") +
    geom_vline(xintercept = c(dates$fal), col = "orange") +
    scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
    theme(text = element_text(size=20)) + ggtitle(species)
  
  multiplot(yrlylon, yrlylat, cols = 1)
  dev.off()
  
  #compare standard error in predicted centroids across years
  pdf(file = paste(dirpath, "/Error_selon-lat", species, ".pdf", sep=""), width = 10, height = 4)
  
  ymax = max(c(migpreds$lat_se, migpreds$lon_se)) 
  lat = ggplot(migpreds, aes(jday, lat_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr), col = "cadetblue") +
    geom_vline(xintercept = c(dates$fal), col = "orange") + ggtitle(paste(species, "Latitude")) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax))
  
  lon = ggplot(migpreds, aes(jday, lon_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
    geom_vline(xintercept = c(dates$spr), col = "cadetblue") +
    geom_vline(xintercept = c(dates$fal), col = "orange") + ggtitle(paste(species, "Longitude")) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax))
  
  multiplot(lat, lon, cols = 2)
  dev.off()
  
  # plot the relationship between lon se and lat se for each year
  pdf(file = paste(dirpath, "/Error_corlon-lat", species, ".pdf", sep=""), width = 10, height = 7)
  
  ymax = max(c(migpreds$lat_se, migpreds$lon_se))
  latlon = ggplot(migpreds, aes(lon_se, lat_se)) + ggtitle(species) +
    geom_point(aes(col = as.factor(month)), alpha = 0.5) + theme_classic() + facet_wrap(~year) +
    theme(text = element_text(size=20)) +
    scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax)) +
    scale_x_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax))
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
  
  sprplot = ggplot(pred_spr, aes(lon, lat, col=year)) + geom_point(size=1) + theme_classic() +
    theme(text = element_text(size=20)) + ggtitle(paste("spring-", species))
  
  falplot = ggplot(pred_fal, aes(lon, lat, col=year)) + geom_point(size=1) + theme_classic() +
    theme(text = element_text(size=20)) + ggtitle(paste("fall-", species))
  
  multiplot(sprplot, falplot, cols = 2)
  dev.off()
  
}
  
  