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


# TODO: Code for migration
#   1. place all daily sightings into hexes
#   2. count which hex has the most
#   3. record central lon-lat for that hex - or record a weighted mean based on the centroid centers
#   4. Calculate great-circle distance between centroids for daily distance
#   5. Calculate 95% confidence bands (longitude) for occurence centroids
#   6. Define beginning of spring migration and end of fall migration
#   7. Link occurrence centroids with environmental data (NDVI, temp, precip, daylength)

## import global equal area hex grid (shared by FAL)
setwd("/Volumes/Elements/eBird/")

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

#To find the POLYFID for the hexes for each observation 
coords = yrdat[,c(10,9)]

# Matches observations with the polygon hexes in the map
ID <- over(SpatialPoints(coords), hexgrid)

#paste original lat-long with the cell info after matching
coords <- cbind(coords, ID) 

#find the number of obs in each cell, doesn't fill with zeroes
table(as.factor(coords2$POLYFID))
#find a way to save or match this up to plot num obs per cell on the map as colors? 
#Use this info to calculate a weighted mean for daily location?
#POLYFID is the ID for each polygon, or hex cell.


#-------------------------- NOTES below, non-functional so far


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




