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
library(insol)


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
t= table(as.factor(coords$POLYFID))

#Merge the count data for the hexes with all the hexes in the map
df = data.frame(POLYFID = names(t), count=as.numeric(t))
df2 = data.frame(POLYFID = unique(hexgrid$POLYFID))
df3 = merge(df2, df, all.x=TRUE)

cols = data.frame(id=c(NA,sort(unique(df3$count))), cols=tim.colors(length(unique(df3$count))), stringsAsFactors=FALSE)
df4 = merge(df3, cols, by.x="count", by.y="id")
df5 = merge(hexgrid, df4, by.x="POLYFID", by.y="POLYFID", all.x=TRUE)

#hexes with no counts are white
df5$cols <- ifelse(is.na(df5$count), "white", df5$cols)
df5 <- df5[order(df5$POLYFID),]



#make a map with hexes colored by the number of times the species was observed in a given hex
jpeg('ChecklistMap.jpg')
plot(hexgrid, col=df5$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1)
axis(side=1)
axis(side=2, las=1)
box()
mtext("Longitude", side=1, cex=1.4, line=2.5)
mtext("Latitude", side=2, cex=1.4, line=2.5)

## legend
vls = sort(unique(round(cols$id/100)*100))
vls[1] = 1
cols2 = tim.colors(length(vls))
legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1, cex=0.75, bty="n",
       col="black", title="Number of checklists", x.intersp=1, y.intersp=0.5, hjust=-1)

dev.off()

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




