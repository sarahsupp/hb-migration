# Subset absent data based on number of present points - for data reduction!
# 1. Clip by flyway

library(sp)
library(maptools)
library(rgeos)
library(raster)
library(rgdal)
#
source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/migration-fxns.r")
source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/hb_RS_functions.R")
#source("/mnt/d/temp/tcormier/820_hummingbirds/migration_study/hb-migration/migration-fxns.r")
#source("/mnt/d/temp/tcormier/820_hummingbirds/migration_study/hb-migration/hb_RS_functions.R")
# for reproducibility
set.seed(55)

wd = "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/ebird/"
#wd = "/mnt/d/temp/tcormier/820_hummingbirds/migration_study/data/ebird/"
setwd(wd)

# study area boundary (should be the western flyway, except for Ruby): **THIS SHOULD BE JUST THE DIRECTORY!
# SCRIPT LOGIC FOR WHICH FLYWAY TO USE DEPENDING ON SPP.
sa.file <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/boundaries/western_flyway_dissolve_envelope.shp"
#sa.file <- "/mnt/d/temp/tcormier/820_hummingbirds/migration_study/boundaries/western_flyway_dissolve_envelope.shp"

#spp (this is how we look for input and name output files)
#spp.list <- c("bchu", "bthu","cahu","rthu","ruhu")
#spp.list <- c("ruhu")
spp <- "bchu"
# List present and absent files
abs_files = list.files(path = sprintf("%sabsent_points/migrants/original/",wd), pattern = "*\\.txt", recursive=FALSE, full.names=TRUE)
pres_files = list.files(path = sprintf("%spresent_points/",wd), pattern = "*\\.txt", recursive=FALSE, full.names=TRUE)

#directory containing files with spring and fall migration dates (from Ecosphere paper)
mig.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/supp_migration/" 

############# One-TIME operation on absent points - then can comment out ############
#Sarah already grouped the present points by Group ID. Need to do the same for absent points
#sent by FAL (ONE-TIME) 

# study area boundary
sa <- readShapePoly(sa.file)

# Assign projection
proj4string(sa) <- CRS("+proj=longlat +datum=WGS84")

spp <- "ruhu"

# loop over spp (note that for ruhu, need some logic to split eBird obs by western flyway,
# then by eastern flyway).
#for (spp in spp.list) {
  print(paste0("working on ", spp))
  pfile <- grep(pattern=spp,x=pres_files,value=T)
  afile <- grep(pattern=spp,x=abs_files,value=T)
  pdata = read.table(pfile, header=TRUE, sep=",", quote='"', fill=TRUE, as.is=TRUE, comment.char="")
  adata = read.table(afile, header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  # Make eBird data spatial.
  coordinates(pdata) <- ~LONGITUDE+LATITUDE
  coordinates(adata) <- ~LONGITUDE+LATITUDE

  # Assign projection
  proj4string(pdata) <- CRS("+proj=longlat +datum=WGS84")
  proj4string(adata) <- CRS("+proj=longlat +datum=WGS84")

  # Clip observations to flyway polygons - **ADD IN LOGIC FOR RUHU (NEED TO DO ONCE BY
  # WESTERN FLYWAY, THEN AGAIN BY EASTERN). 
  pdata.clip <- over(pdata, sa)
  ptid <- na.omit(pdata.clip) 
  pt.poly <- pdata[as.numeric(as.character(row.names(ptid))),]  

  #write out clipped present points as shapefiles
  out_pdata <- paste0(dirname(pfile), "/",unlist(strsplit(basename(pfile), "\\."))[1], "_clipflyway.shp")
  #write.csv(pdata.clip, file=out_pdata, quote=F, row.names=F)
  shapefile(pt.poly, filename=out_pdata)
  #rm(pdata)
  
  #formatting the adata date in various ways
  adata$DAY <- format(adata$DAY, format="%j")
  adata$YEAR <- format(adata$YEAR, format="%Y")
  dt <- paste(adata$DAY, adata$YEAR, sep="_")
  d2 <- as.Date(dt, "%j_%Y")
  adata$DATE <- d2
  adata$MONTH <- format(adata$DATE, "%m")
  adata$DOM <- format(adata$DATE, "%d")
  
  #subset by date, then take random sample, as there are too many records to work with! 
  #Need to figure out appropriate number of samples (accounting for years, spatial stuff, etc.), 
  #but for now...
  #figure out how many records there are in pdata for each year between 2008 - 2013 +30% so we can remove those that are
  #outside of NA and remove duplicate group #s. Chose those yrs based on the first migration paper.

  sampsize <- c(2008:2013)
  names(sampsize) <- c(2008:2013)
  
  for (yr in (c(2008:2013))) {
    print(yr)    
    #first subset present points by year and clip present points by flyway (if file is too
    #big to just clip the whole thing at once).
    #if didn't clip first:    
    pdata.sub <- pdata.clip[pdata.clip$YEAR == yr,] 
    #pdata.clip <- clipPoints(pdata.sub, sa, "LONGITUDE", "LATITUDE")
    #figure out sample size for this year
    #sampsize[names(sampsize)==yr] <- round(nrow(pdata.yr[pdata.yr$YEAR==yr,])*.3 + nrow(pdata.yr[pdata.yr$YEAR==yr,]))
    sampsize[names(sampsize)==yr] <- nrow(pdata.sub)
    
#     #write out clipped present points
#     out_pdata <- paste0(dirname(pfile), "/",unlist(strsplit(basename(pfile), "\\."))[1], "_", yr, "_clipflyway.csv")
#     write.csv(pdata.clip, file=out_pdata, quote=F, row.names=F)
    
    #some cleanup
    rm(pdata.sub)
    #rm(pdata.clip)
    
    #same with adata
    adata.sub <- adata[adata$YEAR == yr,]
    adata.clip <- clipPoints(adata.sub, sa, "LONGITUDE", "LATITUDE")
    rm(adata.sub)
    
    #gets rid of duplicate records that are part of the same group. 
    #Takes only the first record to move to analysis, keeps all single records that are not part of a group.
    adata.dup <- GroupDuplicates(adata.clip)
    rm(adata.clip)
    
    #write out csv of clipped absent points without duplicate groups so I never have to do this again!
    adataclip.out <- paste0(dirname(afile), "/", unlist(strsplit(basename(afile), "\\."))[1], "_", yr, "_clipflyway.csv")
    write.csv(adata.dup, file=adataclip.out, quote=F, row.names=F)
    
    #take random sample from adata.dup based on number of records in pdata for the given year
    adata.samp <- adata.dup[sample(nrow(adata.dup), sampsize[names(sampsize)==yr]),]
    rm(adata.dup)
    
    #write out clipped, sampled points
    out_adata <- paste0(dirname(afile), "/",unlist(strsplit(basename(afile), "\\."))[1], "_", yr, "_sample.csv")
    write.csv(adata.samp, file=out_adata, quote=F, row.names=F)
    rm(adata.samp)   
    }#end yr loop
  
#   #put years back together into one "absent" sampled, clipped file
#   abs_samples <- list.files(path = paste0(dirname(afile), "/"), pattern = "*_sample.csv", full.names=T)
#   abs.df <- data.frame()
#   
#   for (abs in abs_samples) {
#     a <- read.csv(abs)
#     abs.df <- rbind(abs.df, a)
#   }#end putting years back together
#   outAbs_samples <- paste0(dirname(afile), "/",unlist(strsplit(basename(afile), "\\."))[1], "_sample_allYrs.csv")
#   write.csv(abs.df, file=outAbs_samples, row.names=F, quote=F)
#}#end spp loop
  
  
