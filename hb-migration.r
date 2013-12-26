#Code for eBird migration project
# (c) 2013 Sarah Supp 

library(reshape2)
library(plyr)
library(ggplot2)
library(ggmap)

#set working directory
#wd = "C:/Users/sarah/Documents/eBird_data/"
wd = "/Volumes/Elements/eBird/ebd_data/"
setwd(wd)

# make a North America base map
noam = get_map(location = "North America", zoom=3, maptype = "terrain", color = "bw")

# read in eBird data
files = list.files(pattern = "*.txt")
files = files[3]

for (f in 1:length(files)){
  
  require(ggmap)
  require(ggplot2)
  require(plyr)
  require(reshape2)
  
  humdat = read.table(files[f], header=TRUE, sep="\t", fill=TRUE)
    humdat$LONGITUDE = as.numeric(as.character(humdat$LONGITUDE))
    humdat$LATITUDE = as.numeric(as.character(humdat$LATITUDE))

  species = humdat$COMMON.NAME[1]
  years = c(2004:2013)
  
  date = DateConvert(humdat$OBSERVATION.DATE)
  humdat = cbind(humdat, date)
  
  #delete nonsensical data - where something else was recorded in date column
  humdat = humdat[which(humdat$year %in% years),]
  
  #show how many records there are for the species across the years
  PlotRecords(humdat$year, species)
  
  for (y in 1:length(years)){
    yrdat = humdat[which(humdat$year == years[y]),c()]
    
    #plot frequency of sightings per month
    PlotRecords(yrdat$month, species)
    
    #plot where species was sighted within each year
    sitemap = ggmap(noam) + geom_point(aes(LONGITUDE, LATITUDE, col=as.factor(month)), 
                                         data=yrdat) + ggtitle(species)
  }
  rm(list=ls()[!ls() %in% c("f", "files", "noam")])   # clears the memory of everything except the file list and base map
}


