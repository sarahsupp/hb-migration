#Code for eBird migration project
# (c) 2013 -2014 Sarah Supp 

library(ggmap)
library(maptools)
library(sp)
library(raster)

#set working directory
#wd = "C:/Users/sarah/Documents/eBird_data/"
wd = "/Volumes/Elements/eBird/ebd_data/"
setwd(wd)

# make a North America base map
noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")

# read in the north ameica equal area hex grid map (F.A.L.)
#hexgrid = readShapePoly("/Volumes/Elements/eBird/terr_4h6/nw_vector_grid.shp") #quad map
hexgrid = readShapePoly("/Volumes/Elements/eBird/terr_4h6/terr_4h6.shp") #hex map
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

# read in eBird data
files = list.files(pattern = "*.txt")
files = files[6] #c(1,3,4,5,7,9)

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
  source("/Users/sarah/Documents/GitHub/hb-migration/migration-fxns.r")
  
  humdat = read.table(files[f], header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="") #quote="/"'"

  #keep only the columns that we need
  keepcols = c("COMMON.NAME", "SCIENTIFIC.NAME", "OBSERVATION.COUNT", "AGE.SEX", "COUNTRY",
               "COUNTRY_CODE", "STATE_PROVINCE", "COUNTY", "LATITUDE", "LONGITUDE",
               "OBSERVATION.DATE", "TIME.OBSERVATIONS.STARTED", "PROTOCOL.TYPE", "PROJECT.CODE",
               "DURATION.MINUTES", "EFFORT.DISTANCE.KM", "EFFORT.AREA.HA", "NUMBER.OBSERVERS")
  humdat = humdat[,which(names(humdat) %in% keepcols)]

  species = humdat$SCIENTIFIC.NAME[1]
  years = c(2004:2013)
  
  date = DateConvert(humdat$OBSERVATION.DATE)
  humdat = cbind(humdat, date)
    
  #start a new directory
  dirpath = paste("/Volumes/Elements/eBird/ebd_counts/", species, sep="")
    dir.create(dirpath, showWarnings = TRUE, recursive = FALSE)
  
  #show how many records there are for the species across the years, write to txt file
  yeartable = PlotRecords(humdat$year, species)
  write.table(yeartable, file = paste(dirpath, "/", species,".txt",sep=""), row.names=FALSE)
  
  #save a figure of the geographic number of checklists for the species, over all the years
  count = PlotChecklistMap(humdat, hexgrid, dirpath)
  
  for (y in 1:length(years)){
    yrdat = humdat[which(humdat$year == years[y]),]
    
    #plot frequency of sightings per month
    PlotRecords(yrdat$month, species)
    
    #get daily mean location and sd 
    meandat = MeanDailyLoc(yrdat, species)
    cntrdat = FindCentroids(meandat, 7, 5, hexgrid)
    altmeandat = AlternateMeanLocs(yrdat,species,hexgrid)
    
    #get Great Circle distances traveled each day between centroids
    dist = DailyTravel(cntrdat)
    
    #plot where species was sighted within each year
    sitemap = ggmap(noam) + geom_point(aes(LONGITUDE, LATITUDE, col=as.factor(month)), 
                                         data=yrdat) + ggtitle(paste(species, years[y], sep = " "))
    ggsave(sitemap, file=paste(dirpath, "/", species, years[y], ".pdf", sep=""))
    
    #plot species mean location on each julian day with the year
    meanmap = ggmap(noam) + geom_point(aes(meanlon, meanlat, col=as.factor(month)),
                                       data=meandat) + ggtitle(paste(species, years[y], "daily mean loc", sep = " "))
    ggsave(meanmap, file=paste(dirpath, "/", species, years[y],"meanlocs.pdf", sep=""))
    
    #plot mean latitude for each julian day, point size represents number of checklists
    meanlat = ggplot(meandat, aes(jday, meanlat, col=as.factor(month))) + geom_point(aes(size=count)) + 
      ggtitle(paste(species, years[y], "latitudinal migration", sep = " ")) + xlab("Julian Day") + ylab("Mean Latitude") +
      theme_bw() +  scale_x_continuous(breaks = seq(0, 365, by = 25)) +
      scale_y_continuous(breaks = seq(10, 80, by = 5)) + #stat_smooth(method="loess",na.rm=TRUE) #+ 
    geom_errorbar(aes(ymax=maxlat, ymin=minlat)) #+ geom_line(aes(col=as.factor(month)))
    
    ggsave(meanlat, file=paste(dirpath, "/", species, years[y], "meanlat.pdf", sep=""))
    
    rm(list=ls()[ls() %in% c("sitemap", "meanmap", "yrdat")])   # clears the memory of the map and year-level data
  }
  rm(list=ls()[!ls() %in% c("f", "files", "noam", "hexgrid")])   # clears the memory of everything except the file list, iterator, and base map
}


