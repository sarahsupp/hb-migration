# Subset absent data based on number of present points - for data reduction!
# 1. Clip by flyway

# Libraries needed for this code
library(maptools)
library(raster)
library(rgdal)
library(alphahull)
library(lubridate)
#
source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/migration-fxns.r")
source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/hb_RS_functions.R")
#source("/mnt/d/temp/tcormier/820_hummingbirds/migration_study/hb-migration/migration-fxns.r")
#source("/mnt/d/temp/tcormier/820_hummingbirds/migration_study/hb-migration/hb_RS_functions.R")
# for reproducibility
set.seed(55)

############### Inputs #####################
#wd = "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/ebird/"
#wd = "/mnt/d/temp/tcormier/820_hummingbirds/migration_study/data/ebird/"
setwd(wd)

#Directory containing eBird observation files
data.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/ebird/"

# Directory containing eastern and western flyway study area boundaries (all spp western flyway, except for Ruby, and Rufous is both): 
# SCRIPT LOGIC FOR WHICH FLYWAY TO USE DEPENDING ON SPP.
sa.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/boundaries/"
#sa.file <- "/mnt/d/temp/tcormier/820_hummingbirds/migration_study/boundaries/western_flyway_dissolve_envelope.shp"

#spp (this is how we look for input and name output files)
#spp.list <- c("bchu", "bthu","cahu","rthu","ruhu_east", "ruhu_west")
#spp.list <- c("ruhu")

#directory containing files with spring and fall migration dates (from Ecosphere paper)
mig.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/supp_migration/" 

############################################
# Sarah already grouped the present points by Group ID. Need to do the same for absent points
# sent by FAL (ONE-TIME) 

# List present and absent files
abs_files = list.files(path = sprintf("%sabsent_points/migrants/original/",data.dir), pattern = "*\\.txt", recursive=FALSE, full.names=TRUE)
pres_files = list.files(path = sprintf("%spresent_points/",data.dir), pattern = "*\\.txt", recursive=FALSE, full.names=TRUE)

spp <- "ruhu_west"
# For each species, lets bring in the eBird data, and clip it to the migration season. For this task, migration will be defined 
# within each year, using the dates from the analysis in Supp et al. 2015 (Ecosphere). We will begin one week prior to the onset 
# of spring migration date. Counting forward, we will end one week after the week that contains the end of fall migration. The data 
# will be returned to a new folder called "ebird_weeks" in the data directory, where it can be used later.
# loop over spp, putting ruhu in twice so we can split by eastern flyway and western. 

for (spp in spp.list) {
  print(paste0("working on ", spp))
  # Open present and absent data files for spp and make them spatial (important later for spatial subsetting).
  pfile <- grep(pattern=spp,x=pres_files,value=T)
  afile <- grep(pattern=spp,x=abs_files,value=T)
  pdata = read.table(pfile, header=TRUE, sep=",", quote='"', fill=TRUE, as.is=TRUE, comment.char="")
  adata = read.table(afile, header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  # Some structuring of the tables
  pdata$MONTH = factor(pdata$MONTH, levels=c(1:12), ordered=TRUE)
  
  # format date in absence data 
  d <- paste(adata$DAY, adata$YEAR, sep="_")
  adata$DATE <- as.Date(d, "%j_%Y")
  adata$MONTH <- month(adata$DATE)
  
  # Now subset by dates - starting with only years >2007
  pdata.subyr <- pdata[pdata$YEAR > 2007,]
  adata.subyr <- adata[adata$YEAR > 2007,]
  
  #Read in migration dates
  mig.files <- list.files(mig.dir, pattern="*.txt", full.names = T)
  mig.spp <- grep(pattern=spp,x=mig.files,value=T)
  mig.dates <- read.table(mig.spp, header=T, sep=" ")
  names(mig.dates) <- c("spr_begin", "peak_lat", "fal_begin", "fal_end", "species", "year")
  
  years = unique(pdata.subyr$YEAR)
  
  # trim data from the beginning of the year up to 1 week before the onset of spring migration and at least 1 
  # full week past the end of autumn migration. Start is defined as one week before the onset of spring migration;
  # End is defined as one week past the week that contains the end of autumn migration.
  # New columns will be added for increment (number of days since "start") and week (weeks since "start")
  sf.pres <- as.data.frame(matrix(data=NA, nrow=0, ncol=ncol(pdata.subyr)+2))
  names(sf.pres) <- c(names(pdata.subyr), "increment", "week")
  
  # This is from Sarah's code. *Does it make sense to do the same for the absence data
  # at this point as well? Think about it!
  for (yr in years){
    dat.yr <- pdata.subyr[pdata.subyr$YEAR == yr,]
    start <-  mig.dates$spr_begin[mig.dates$year == yr] - 7
    n.weeks <- floor(((mig.dates$fal_end[mig.dates$year == yr] - start)/7) + 2)
    end <- start + n.weeks * 7
    dat.sf <- dat.yr[dat.yr$DAY >= start & dat.yr$DAY <= end,]
    dat.sf$increment <- dat.sf$DAY - start + 1
    dat.sf$week <- ceiling(dat.sf$increment/7)
    
    sf.pres=rbind(sf.pres, dat.sf)
    
  }
  # Reset the rownames to be 1: length of sf.pres - otherwise causes indexing issues when doing spatial overlay.
  row.names(sf.pres) <- NULL
  #sf2 = sf[,c(1,2,13,14,3,4,18,5,17,19)] #take out all the extra columns that we don't really need here
  #write to file
  sf.pres.file <- paste(data.dir, "present_points//ebird_weeks/", spp, "_by_week.txt", sep="")
  write.table(sf.pres, file =sf.pres.file,  sep=",", row.names=FALSE, append=FALSE)
  #rm(c(dat, dat.yr, dat.sf, sf2))   
  
  # Now for spatial subsetting by flyway. Make eBird data spatial. 
  coordinates(sf.pres) <- ~LONGITUDE+LATITUDE
  #coordinates(adata) <- ~LONGITUDE+LATITUDE
  
  # Assign projection
  proj4string(sf.pres) <- CRS("+proj=longlat +datum=WGS84")
#   proj4string(adata) <- CRS("+proj=longlat +datum=WGS84")
  
  # flyway boundaries for the spp - will open both of them, since they are small. Otherwise, might consider adding 
  # logic for which one to open based on spp.
  sa.files <- c(paste0(sa.dir, "western_flyway_dissolve_envelope.shp"),paste0(sa.dir, "eastern_flyway_dissolve_envelope.shp"))
  sa.west <- readShapePoly (sa.files[1])
  proj4string(sa.west) <- CRS("+proj=longlat +datum=WGS84")
  sa.east <- readShapePoly (sa.files[2])
  proj4string(sa.east) <- CRS("+proj=longlat +datum=WGS84")
  
  #Can delete after you're sure you don't need this logic!
#   if (spp != "ruhu" & spp != "rthu") {
#     sa.file <- paste0(sa.dir, "western_flyway_dissolve_envelope.shp")
#     # study area boundary
#     sa <- readShapePoly(sa.file)
#     # Assign projection
#     proj4string(sa) <- CRS("+proj=longlat +datum=WGS84")
#   } else if (spp == "rthu") {
#       sa.file <- paste0(sa.dir, "eastern_flyway_dissolve_envelope.shp")
#       # study area boundary
#       sa <- readShapePoly(sa.file)
#   } else if (spp == "ruhu") {
#       # Need to clip by both flyways
#       sa.files <- c(paste0(sa.dir, "western_flyway_dissolve_envelope.shp"),sa.file <- paste0(sa.dir, "eastern_flyway_dissolve_envelope.shp"))
#       
#   }# end flyway if
  


  # Clip observations to flyway polygons - **ADD IN LOGIC FOR RUHU (NEED TO DO ONCE BY
  # WESTERN FLYWAY, THEN AGAIN BY EASTERN). sffw = spring/fall/flyway
  out_sffw.pres <- paste0(dirname(sf.pres.file), "/",unlist(strsplit(basename(sf.pres.file), "\\."))[1], "_clipflyway.csv")

  if (spp != "ruhu_east" & spp != "rthu") {
    ptsSA <- as.data.frame(PointsPolyIntersect(sf.pres, sa.west))
    #write out clipped present points as shapefiles
    write.csv(out_sffw.pres, file=out_pdata, quote=F, row.names=F)
    #clean up
    rm(pdata)
    # formatting the adata date in various ways (EDIT - DON'T DO ANYTHING WITH ABSENT DATA JUST YET - COMMENT OUT FOR NOW UNTIL
    # I DECIDE IF IT MAKES SENSE TO DO THIS NOW OR WAIT UNTIL THE ALPHA HULLS PART)
    #adata <- formatAbsent(adata)
    # Try clipping to SA
    #adata.clip <- PointsPolyIntersect(adata, sa.west)
    #rm(adata)
  } else if (spp == "ruhu_east" | spp == "rthu") {
      ptsSA <- as.data.frame(PointsPolyIntersect(sf.pres, sa.east)
      #write out clipped present points as shapefiles
      write.csv(pdata.clip, file=out_pdata, quote=F, row.names=F)
      #clean up
      rm(pdata)
      # formatting the adata date in various ways (EDIT - DON'T DO ANYTHING WITH ABSENT DATA JUST YET - COMMENT OUT FOR NOW UNTIL
      # I DECIDE IF IT MAKES SENSE TO DO THIS NOW OR WAIT UNTIL THE ALPHA HULLS PART)
      #adata <- formatAbsent(adata)
      # Try clipping to SA
      #adata.clip <- PointsPolyIntersect(adata, sa.east)
      #rm(adata)
  } else {
    print(paste("accepted spp include bchu,", "bthu,","cahu,","rthu,","ruhu_east,", "or ruhu_west. You entered:", spp))
  } #end clipping to flyways
  
  #Calculate alpha hulls for present points for each spp - adapted from Frank LaSorte's code.
  # read in clipped week data  #FIXME: reads in with extra quotes and X. in colname
  #clip <- read.table(paste0(data.dir, "ebird_weeks/", spp, "_by_week.txt"), header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  clip <- ptsSA
  # loop through years
  # loop through weeks
  # save alpha hulls, linked to each week and year to query the complete checklist.
  #    new dataframe should have columns for year, week, startofweekday, alphahull polygon coords.
  
  ## https://stat.ethz.ch/pipermail/r-sig-geo/2012-March/014409.html
  source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/ahull-to-polygon.r")
  #Convert sf.pres back to dataframe
  #sf.pres.df <- as.data.frame(sf.pres)
  #for the purposes of making the alphahull, need to remove duplicate locations.
  locs <- as.data.frame(cbind(sf.pres@coords[,1], y=sf.pres@coords[,2]))
  names(locs) <- c("LONGITUDE", "LATITUDE")
  uniq.locs <- subset(locs, !duplicated(locs))
  hull <- ahull(uniq.locs, alpha = 0.1)
  
  #plot for grins ***STOPPED HERE 12/19 - still need to write in loops - just testing here to get alpha hulls working!***
  plot(hull)

  phull <- ah2sp(hull)
  phull <- as(phull, "SpatialPolygons")
  proj4string(phull) <- CRS("+proj=longlat +datum=WGS84")
  writeOGR(phull, dsn = “.”, layer = “ahull”, driver = “ESRI Shapefile”)
  writeSpatialShape(phull, “ahull”)
  
  
  
  







  ## OLD SAMPLING ##
  #subset by date, then take random sample, as there are too many records to work with! 
  #Need to figure out appropriate number of samples (accounting for years, spatial stuff, etc.), 
  #but for now...
  #figure out how many records there are in pdata for each year between 2008 - 2013 +30% so we can remove those that are
  #outside of NA and remove duplicate group #s. Chose those yrs based on the first migration paper.

#   sampsize <- c(2008:2013)
#   names(sampsize) <- c(2008:2013)
#   
#   for (yr in (c(2008:2013))) {
#     print(yr)    
#     #first subset present points by year and clip present points by flyway (if file is too
#     #big to just clip the whole thing at once).
#     #if didn't clip first:    
#     pdata.sub <- pdata.clip[pdata.clip$YEAR == yr,] 
#     #pdata.clip <- clipPoints(pdata.sub, sa, "LONGITUDE", "LATITUDE")
#     #figure out sample size for this year
#     #sampsize[names(sampsize)==yr] <- round(nrow(pdata.yr[pdata.yr$YEAR==yr,])*.3 + nrow(pdata.yr[pdata.yr$YEAR==yr,]))
#     sampsize[names(sampsize)==yr] <- nrow(pdata.sub)
#     
# #     #write out clipped present points
# #     out_pdata <- paste0(dirname(pfile), "/",unlist(strsplit(basename(pfile), "\\."))[1], "_", yr, "_clipflyway.csv")
# #     write.csv(pdata.clip, file=out_pdata, quote=F, row.names=F)
#     
#     #some cleanup
#     rm(pdata.sub)
#     #rm(pdata.clip)
#     
#     #same with adata
#     adata.sub <- adata[adata$YEAR == yr,]
#     adata.clip <- clipPoints(adata.sub, sa, "LONGITUDE", "LATITUDE")
#     rm(adata.sub)
#     
#     #gets rid of duplicate records that are part of the same group. 
#     #Takes only the first record to move to analysis, keeps all single records that are not part of a group.
#     adata.dup <- GroupDuplicates(adata.clip)
#     rm(adata.clip)
#     
#     #write out csv of clipped absent points without duplicate groups so I never have to do this again!
#     adataclip.out <- paste0(dirname(afile), "/", unlist(strsplit(basename(afile), "\\."))[1], "_", yr, "_clipflyway.csv")
#     write.csv(adata.dup, file=adataclip.out, quote=F, row.names=F)
#     
#     #take random sample from adata.dup based on number of records in pdata for the given year
#     adata.samp <- adata.dup[sample(nrow(adata.dup), sampsize[names(sampsize)==yr]),]
#     rm(adata.dup)
#     
#     #write out clipped, sampled points
#     out_adata <- paste0(dirname(afile), "/",unlist(strsplit(basename(afile), "\\."))[1], "_", yr, "_sample.csv")
#     write.csv(adata.samp, file=out_adata, quote=F, row.names=F)
#     rm(adata.samp)   
#     }#end yr loop
#   
# #   #put years back together into one "absent" sampled, clipped file
# #   abs_samples <- list.files(path = paste0(dirname(afile), "/"), pattern = "*_sample.csv", full.names=T)
# #   abs.df <- data.frame()
# #   
# #   for (abs in abs_samples) {
# #     a <- read.csv(abs)
# #     abs.df <- rbind(abs.df, a)
# #   }#end putting years back together
# #   outAbs_samples <- paste0(dirname(afile), "/",unlist(strsplit(basename(afile), "\\."))[1], "_sample_allYrs.csv")
# #   write.csv(abs.df, file=outAbs_samples, row.names=F, quote=F)
# #}#end spp loop
#   
#   
