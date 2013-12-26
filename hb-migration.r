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

# read in eBird data
bchu = read.table('ebd_bkchum_200401_201312_relNov-2013.txt',
                  header=TRUE, sep="\t", fill=TRUE)
ruhu = read.table('ebd_rufhum_200401_201312_relNov-2013.txt', 
                  header=TRUE, sep="\t", fill=TRUE, nrows=10)

files = list.files(pattern = "*.txt")

for (f in 1:length(files)){
  humdat = read.table(files[f], header=TRUE, sep="\t", fill=TRUE)

  species = humdat$COMMON.NAME[1]
  
  date = DateConvert(humdat$OBSERVATION.DATE)
  humdat = cbind(humdat, date)
  
  #delete nonsensical data - where something else was recorded in date column
  humdat = humdat[which(humdat$year %in% 2004:2013),]
  
  #show how many records there are for the species across the years
  PlotRecords(humdat$year, species)

}


