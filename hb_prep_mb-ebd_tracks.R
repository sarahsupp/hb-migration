library(raster)
library(maptools)
library(rgdal)

##########################
# User input variables
#shapefiles associated with North America grid
grid.poly <- "C:/Share/tcormier/hummingbirds/migration_study/gridded_references/na_p3deg_albers_32k_unique_poly.shp"
grid.points <- "C:/Share/tcormier/hummingbirds/migration_study/gridded_references/na_p3deg_albers_32k_unique_points.shp"

#ebird daily observation data
ebd.shp <- "C:/Share/tcormier/hummingbirds/migration_study/data/ebird/ebd_rfuhum_2004_201312_relNove-2013_albers.shp"

#Output movebank tracks file (directory):
trackdir <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/track_csvs/"


##########################
#read in shapefiles (try first with maptools)
poly <- readShapePoly(grid.poly)
pts <- readShapePoints(grid.points)
ebd <- readShapePoints(ebd.shp)

#Intersect ebd with poly, which will give me the IDs of the points I need to submit to movebank.
#This is for model building, as we are only looking at areas (grid cells) in which the hb
#was observed. Later, for temporal prediction, we'll have to use the whole envelope, within which
#hbs were observed and not observed.

#in fact, when I do this "for real" later (soon?), I will just request the whole area (need to figure out
#how to determine what the "whole area" is) for every day and extract from that the ebird observations.  That way,
#I get all of the data in one swoop and don't ask for the same thing multiple times from movebank. But before that 
#happens, I need to tile this probably.  Also, haven't settled on a resoltion..though this res (32km) is finer
#than the hex analysis. 250m might be an unweildy amount of data.

#But I digress...this gives me a list of grid IDs that intersect all ebird obs (all years)
polyebd <- over(ebd, poly)

#Now get unique values
polyebd.un <- unique(polyebd$GRIDCODE)

#Pull out just 2012 for now
ebd2012 <- ebd[ebd$OBSERVAT_1 > "2011-12-31" & ebd$OBSERVAT_1 < "2013-01-01",] 
polyebd2012 <- over(ebd2012, poly)


#format as text file to submit to movebank
ebd.names <- c("timestamp", "location-long", "location-lat", "height-above-ellipsoid")
ebd.2012.csv <- as.data.frame(matrix(data=NA, nrow=length(polyebd2012$GRIDCODE),ncol=length(ebd.names),))

#head(ebd.2012.csv)
timestamp <- paste(ebd2012$OBSERVAT_1, "12:00:00.000")
ebd.2012.csv[,1] <- timestamp
ebd.2012.csv[,2] <- ebd2012$LONGITUDE
ebd.2012.csv[,3] <- ebd2012$LATITUDE
ebd.2012.csv[,4] <- ""

#Assign names last bc R doesn't like dashes in column names, so after I'm done
#fiddling with the columns, assign the names.
names(ebd.2012.csv) <- ebd.names

#write out tracks csv as well as an ID csv that I can use to join back to the grid after annotation.
out.tracks <- paste0(trackdir, unlist(strsplit(basename(ebd.shp), "\\."))[1], "_tracks.csv")
out.ids <- paste0(trackdir, unlist(strsplit(basename(ebd.shp), "\\."))[1], "_ids.csv")
write.csv(ebd.2012.csv, file=out.tracks, quote=F, row.names=F)
write.csv(polyebd2012$GRIDCODE, file=out.ids, quote=F, row.names=F)
