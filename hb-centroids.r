# Code for tracking hummingbird migration
# Modified from code for La Sorte et al. 2013. (Population-level scaling of 
# avian migration speed with body size and migration distance for powered fliers, Ecology)
# Sarah R. Supp (c) 2014

library(spaa)
library(fields)
library(colorRamps)
library(maptools)
library(fossil)
library(PBSmodelling)
library(aspace)
library(pastecs)
library(raster)
library(mgcv)
library(adehabitat)
library(SDMTools)
library(xlsx)
library(gamm4)

#Here are the dates for the weeks - the temporal design for SRD predictions.
# -----------------------
n.intervals.per.year = 365
# ------------
begin.jdate = 1
end.jdate = 365
year.list = c(2009)

n.intervals = n.intervals.per.year  + 1
jdate.ttt.seq = seq(from = begin.jdate, to = end.jdate, length = (n.intervals))

if (n.intervals > 1) {
  jdate.ttt.seq = round((jdate.ttt.seq[2:n.intervals] +
                            jdate.ttt.seq[1:(n.intervals-1)]) / 2)
}


year.seq = NULL
jdate.seq = NULL
for (iii.year in year.list){
  year.seq = c(year.seq,rep(iii.year, length(jdate.ttt.seq)))
  jdate.seq = c(jdate.seq, jdate.ttt.seq)
}

#change julian date into "YYYY-MM-DD" and "Mon-dd" strings
posix.dates = strptime( x=paste(2009,"-",jdate.ttt.seq, sep=""),"%Y-%j")
date.string = paste(months(posix.dates, abbreviate=T), posix.dates$mday,sep="_")

## import global equal area hex grid
setwd("/Volumes/Elements/eBird/terr_4h6")
hex.grd = readShapePoly("nw_vector_grid.shp")
hex.lon.lat = data.frame(POLYFID = hex.grd$POLYFID, LON = hex.grd$LONGITUDE, LAT = hex.grd$LATITUDE)
hex.lon.lat2 = unique(hex.lon.lat)

## Nature Serve centroids - a datafile with centroids already calculated? Is it just mean Lat and Long?
setwd("pathname/analysis")
ns.cntr = read.table("NatureServe_centroids2.txt", sep="\t", header=TRUE, as.is=TRUE)

## migration distance
sp2 = sort(unique(ns.cntr$scientific2))
dst = NULL

#use geodist to calculate Great Circle distance between points
for(i in 1:length(sp2)){
  ns.cntr2 <- ns.cntr[ns.cntr$scientific2 %in% sp2[i],]
  dst = c(dst, geodist(ns.cntr2[1,"lon"], ns.cntr2[1,"lat"], ns.cntr2[2,"lon"], ns.cntr2[2,"lat"]))
}

#TODO: should we include all jdates, e.g., jdate with no obs are simply NA?
mig.dst.all = data.frame(scientific2 = sp2, distance = dst, stringsAsFactors = FALSE)
mig.dst.all = na.omit(mig.dst.all)
mig.dst.all = mig.dst.all[mig.dst.all$distance > 0, ]    ## migration distance >0
sp = sort(unique(mig.dst.all$scientific2))




