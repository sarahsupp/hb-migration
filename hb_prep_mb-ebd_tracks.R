#Purpose: Formats eBird data "tracks" for submission to Movebank.
#This script can format both present and absent data based on the format we 
#are using as of 6/25/14.
#
#Author: Tina Cormier
#
#Date: 06/25/2014
#
#Status: Runs but needs to be cleaned up!
###################################################
library(maptools)
library(raster)

#ebird daily observation data
#ebd.presfiles <- "C:/Share/tcormier/hummingbirds/migration_study/data/ebird/present_points/"
ebd.presfiles <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/ebird/present_points/"

#ebd.absfiles <- "C:/Share/tcormier/hummingbirds/migration_study/data/ebird/absent_points/"
ebd.absfiles <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/ebird/absent_points/"

#study area boundary (should be the western flyway, except for Ruby):
sa.file <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/boundaries/western_flyway_dissolve.shp"
  
#Output movebank tracks file (directory):
#trackdir <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/track_csvs/"
trackdir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/track_csvs/"
#lag in days
lag=0

#spp.list <- c("bchu", "bthu","cahu","rthu","ruhu")
spp.list <- c("rthu")
###################################################
 

#function to prep present data
#lag is in days - 0 if no lag.
prepPres <- function(spp, spp_tbl, lag, outdir) {
  #format as text file to submit to movebank
  ebd_names <- c("timestamp", "location-long", "location-lat", "height-above-ellipsoid")
  
  #Break down tbl by yr for pres data (already done for abs data)
  for (yr in c(2008:2013)) {
    ebd_yr <- spp_tbl[spp_tbl$YEAR == yr,]
    ebd <- as.data.frame(matrix(data=NA, nrow=nrow(ebd_yr),ncol=length(ebd_names),))
    
    #format date
    d <- paste(ebd_yr$DAY, ebd_yr$YEAR, sep="_")
    d2 <- as.Date(d, "%j_%Y")
    
    #incorporate lag if requested
    d2 <- d2-lag
    ts <- paste0(ebd_yr$TIME, ".000")
    dt <- paste(d2, ts)
    
    #fill in ebd df
    ebd[,1] <- dt
    ebd[,2] <- ebd_yr$LONGITUDE
    ebd[,3] <- ebd_yr$LATITUDE
    ebd[,4] <- ""
    
    #Assign names last bc R doesn't like dashes in column names, so after I'm done
    #fiddling with the columns, assign the names.
    names(ebd) <- ebd_names
    
    #write table
    outpres <- paste0(outdir,"/", spp, "_pres_", yr, "_lag",lag, ".csv") 
    write.csv(ebd, file=outpres, quote=F, row.names=F)
  }#end yr loop
}#end prepPres


#function to prep Absent data
#lag is in days - 0 if no lag.
#could consolidate pres and abs function into one with some logic. 
#Just rushing right now!
prepAbs <- function(spp, abs_tbl, lag, outdir, year) {
  #format as text file to submit to movebank
  ebd_names <- c("timestamp", "location-long", "location-lat", "height-above-ellipsoid")
  ebd <- as.data.frame(matrix(data=NA, nrow=nrow(abs_tbl),ncol=length(ebd_names),))
  
  #format date - 
  d <- paste(abs_tbl$DAY, abs_tbl$YEAR, sep="_")
  d2 <- as.Date(d, "%j_%Y")
  
  #incorporate lag if requested
  d2 <- d2-lag
  #Frank did not provide time with the absent points. Prob not important in this analysis.
  ts <- "12:00:00.000"
  dt <- paste(d2, ts)
  
  #fill in ebd df
  ebd[,1] <- dt
  ebd[,2] <- abs_tbl$LONGITUDE
  ebd[,3] <- abs_tbl$LATITUDE
  ebd[,4] <- ""
  
  #Assign names last bc R doesn't like dashes in column names, so after I'm done
  #fiddling with the columns, assign the names.
  names(ebd) <- ebd_names
  
  #write table
  outAbs <- paste0(outdir,"/", spp, "_abs_", year,"_lag",lag, ".csv") 
  write.csv(ebd, file=outAbs, quote=F, row.names=F)
}#end Abs function

#This can be a loop over multiple spp in spp.list
#spp <- spp.list[1]

for (spp in spp.list) {
  #prep and write out pres files
  pfile <- paste0(ebd.presfiles, spp, ".txt")
  tbl <- read.table(pfile, header=TRUE, sep=",", quote='"', fill=TRUE, as.is=TRUE, comment.char="")
  
  #first clip pres points to study area boundary (western flyway, except for rthu)
  
  
  #write out pres track files
  trackdir2 <- paste0(trackdir, spp, "/")
  prepPres(spp, tbl, 0, trackdir2)
  
  #prep and write out Abs files
  pattern=paste0(spp, "_[0-9]{4}_sub\\.csv")
  afiles <- list.files(path=ebd.absfiles, pattern=pattern, full.names=T)
  
  for (abs_y in afiles) {
    abs_tbl <- read.table(abs_y, header=TRUE, sep=",", quote='"', fill=TRUE, as.is=TRUE, comment.char="")
    year <- regexpr(pattern="[0-9]{4}", text=abs_y)
    year <- regmatches(abs_y, year)
    
    #Write out Abs track files
    prepAbs(spp, abs_tbl, 0, trackdir2, year)
  }
  
  #write out pres track files
  #prepPres(spp,tbl,0,trackdir)
}#end spp loop  


