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

  #hex grid map
  plot(hexgrid, xlim=c(-170,-50), ylim=c(10,80), col="lightblue", lwd=0.25, border="gray10")
  axis(side=1)
  axis(side=2, las=1)
  box()
  mtext("Longitude", side=1, cex=1.4, line=2.5)
  mtext("Latitude", side=2, cex=1.4, line=2.5)

# make a North America base map
noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")

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
    dir.create(dirpath, showWarnings = TRUE, recursive = FALSE)
  
  #show how many records there are for the species across the years, write to txt file
  yeartable = PlotRecords(humdat$year, species)
  write.table(yeartable, file = paste(dirpath, "/", species,".txt",sep=""), row.names=FALSE)
  
  write.table(NULL, file = paste(dirpath, "/", "migration", species,".txt",sep=""), row.names=FALSE)
  
  #save a figure of the geographic number of checklists for the species, over all the years
  count = PlotChecklistMap(humdat, hexgrid, dirpath)
  
  for (y in 1:length(years)){
    yrdat = humdat[which(humdat$year == years[y]),]
    
    #plot frequency of sightings per month
    PlotRecords(yrdat$month, species)
    
    #get daily mean location and sd 
    #meandat = MeanDailyLoc(yrdat, species)     #----may not need this version
    #cntrdat = FindCentroids(meandat, 7, 5, hexgrid)      #----may not need this version
    altmeandat = AlternateMeanLocs(yrdat,species,hexgrid)
    
    #grab dates for migration
    migration = GetMigrationDates(altmeandat)

    m = c(migration[[1]], migration[[2]])
    write.table(m, file = paste(dirpath, "/", "migration", species, ".txt",sep=""), 
                append=TRUE, col.names=FALSE, row.names=FALSE)
    
    #use GAM model to predict daily location along a smoothing line
    preds = EstimateDailyLocs(altmeandat)
    
    #get Great Circle distances traveled each day between predicted daily locations
    dist = DailyTravel(preds, 4, 5, species, year)
    
    #plot smoothed migration trajectory for the species and year
    mig_path = PlotMigrationPath(preds, noam, species, years[y])
    ggsave(mig_path, file=paste(dirpath, "/", "migration", species, years[y], ".pdf", sep=""))
    
    #plot occurrences with lines showing beginning and end of migration
    PlotOccurrences(altmeandat, species, migration[[1]], migration[[2]])
    
    rm(list=ls()[ls() %in% c("sitemap", "meanmap", "yrdat", "altmeandat", "migration", "preds", "dist", "mig_path")])   # clears the memory of the map and year-level data
  }
  rm(list=ls()[!ls() %in% c("f", "files", "noam", "hexgrid", "effort")])   # clears the memory of everything except the file list, iterator, and base map
}


