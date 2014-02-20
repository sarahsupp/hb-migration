#Code for eBird migration project
# (c) 2013 -2014 Sarah Supp 

library(ggmap)
library(maptools)
library(sp)
library(raster)

#set working directory
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data"
setwd(wd)

# read in summary of effort data (Number of eBird checklists submitted per day per year)
effort = read.table("cell_effort.txt", header=TRUE)

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
  dirpath = paste("C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/Figures/", species, sep="")
 #   dir.create(dirpath, showWarnings = TRUE, recursive = FALSE)
  
  #show how many records there are for the species across the years, write to txt file
  yeartable = PlotRecords(humdat$year, species)
  write.table(yeartable, file = paste(dirpath, "/", species,".txt",sep=""), row.names=FALSE)
  
  write.table(NULL, file = paste(dirpath, "/", "migration", species,".txt",sep=""), row.names=FALSE)
  
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
      write.table(migdates, file = paste(dirpath, "/", "migration", species, ".txt",sep=""), 
                  append=TRUE, col.names=FALSE, row.names=FALSE)
      write.table(migspeed, file = paste(dirpath, "/", "speed", species, ".txt",sep=""), 
                  append=TRUE, col.names=FALSE, row.names=FALSE)
      
      ggplot(pred_data, aes(jday, lon, col=year)) + geom_point() + theme_classic() +
        geom_vline(xintercept = c(migdates$spr), col = "cadetblue") +
        geom_vline(xintercept = c(migdates$fal), col = "orange") +
        scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
        theme(text = element_text(size=20)) + ggtitle(species)
      ggsave(filename = paste(dirpath, "/", "lon_allyears", species,".jpeg",sep=""))
      
      ggplot(pred_data, aes(jday, lat, col=year)) + geom_point() + theme_classic() +
        geom_vline(xintercept = c(migdates$spr), col = "cadetblue") +
        geom_vline(xintercept = c(migdates$fal), col = "orange") +
        scale_x_continuous(breaks = seq(0, 365, by = 30)) + 
        theme(text = element_text(size=20)) + ggtitle(species)
      ggsave(filename = paste(dirpath, "/", "lat_allyears", species,".jpeg",sep=""))
      
      pred_spr = pred_data[which(pred_data$jday > mean(migdates$spr) & pred_data$jday < median(c(migdates$spr, migdates$fal))),]
      pred_fal = pred_data[which(pred_data$jday < mean(migdates$fal) & pred_data$jday > median(c(migdates$spr, migdates$fal))),]
      
      ggplot(pred_spr, aes(lon, lat, col=year)) + geom_point() + theme_classic() +
        theme(text = element_text(size=20)) + ggtitle(paste("spring -", species))
      ggsave(filename = paste(dirpath, "/spr_allyears", species,".jpeg",sep=""))
     
      ggplot(pred_fal, aes(lon, lat, col=year)) + geom_point() + theme_classic() +
        theme(text = element_text(size=20)) + ggtitle(paste("fall -", species))
      ggsave(filename = paste(dirpath, "/fal_allyears", species,".jpeg",sep=""))
      
      #compare sd across years
      pdf(file = paste(dirpath, "/se_latANDlon", species, ".pdf", sep=""), width = 10, height = 4)
      
      ymax = ceiling(max(c(pred_data$lat_se, pred_data$lon_se)))  
       lat = ggplot(pred_data, aes(jday, lat_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
        geom_vline(xintercept = c(migdates$spr), col = "cadetblue") +
        geom_vline(xintercept = c(migdates$fal), col = "orange") + ggtitle(paste(species, "Latitude")) +
        scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax)) 
      
        lon = ggplot(pred_data, aes(jday, lon_se, col=as.factor(year))) + geom_point(size=1) + theme_classic() +
        geom_vline(xintercept = c(migdates$spr), col = "cadetblue") +
        geom_vline(xintercept = c(migdates$fal), col = "orange") + ggtitle(paste(species, "Longitude")) +
        scale_y_continuous(breaks = seq(0, ymax, by = 0.25), limits = c(0, ymax)) 
      
       multiplot(lat, lon, cols = 2)
      dev.off()
      
      #compare location across the years using mean and sd
      jdays = sort(unique(pred_data$jday))
      patherr = data.frame("jday"=1,"meanlat"=1, "sdlat"=1, "meanlon"=1, "sdlon"=1)
      outcount = 1
      for (j in 1:length(jdays)){
        tmp = pred_data[which(pred_data$jday == j),]
        meanlat = mean(tmp$lat)
        sdlat = sd(tmp$lat)
        meanlon = mean(tmp$lon)
        sdlon = sd(tmp$lon)
        patherr[outcount,] = c(j, meanlat, sdlat, meanlon, sdlon)
        outcount = outcount + 1
      }
      ggplot(patherr, aes(jday, sdlat)) + geom_point() + theme_classic()
      ggsave(filename = paste(dirpath, "/sdlat", species,".jpeg",sep=""))
      ggplot(patherr, aes(jday, sdlon)) + geom_point() + theme_classic()
      ggsave(filename = paste(dirpath, "/sdlon", species,".jpeg",sep=""))
      ggplot(patherr, aes(sdlon, sdlat, col=jday)) + geom_point() + theme_classic()
      ggsave(filename = paste(dirpath, "/sdlon-lat", species,".jpeg",sep=""))
    }
    
    rm(list=ls()[ls() %in% c("sitemap", "meanmap", "yrdat", "altmeandat", "migration", "preds", "dist", "mig_path")])   # clears the memory of the map and year-level data
  }
  rm(list=ls()[!ls() %in% c("f", "files", "noam", "hexgrid", "effort", "pred_data", "migdates")])   # clears the memory of everything except the file list, iterator, and base map
}


