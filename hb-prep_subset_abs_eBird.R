# Cleans up the data needed for the hummingbird migration project. 
# combines presence and absence points into one file.

library(sp)
library(maptools)
library(rgeos)
#
source("C:/Share/tcormier/scripts/git_repos/hb_migration/migration-fxns.r")

set.seed(55)

wd = "C:/Share/tcormier/hummingbirds/migration_study/data/ebird/"
setwd(wd)

#study area boundary:
sa.file <- "C:/Share/tcormier/hummingbirds/boundaries/migration_boundaries/US_CAN_MEX_dissolve_singleFeature.shp"

#spp (this is how we look for input and name output files)
#spp.list <- c("bchu", "bthu","cahu","rthu","ruhu")
spp.list <- c("ruhu")
# List present and absent files
abs_files = list.files(path = sprintf("%sabsent_points/",wd), pattern = "*\\.txt", recursive=TRUE, full.names=TRUE)
pres_files = list.files(path = sprintf("%spresent_points/",wd), pattern = "*\\.txt", recursive=TRUE, full.names=TRUE)

############# One-TIME operation on absent points - then can comment out ############
#Sarah already grouped the present points by Group ID. Need to do the same for absent points
#sent by FAL (ONE-TIME) 

#study area boundary
sa <- readShapePoly(sa.file)

#spp <- "ruhu"
for (spp in spp.list) {
  print(paste0("working on ", spp))
  pfile <- grep(pattern=spp,x=pres_files,value=T)
  afile <- grep(pattern=spp,x=abs_files,value=T)
  pdata = read.table(pfile, header=TRUE, sep=",", quote='"', fill=TRUE, as.is=TRUE, comment.char="")
  adata = read.table(afile, header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  
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

  yrs <- c(2008:2013)
  names(yrs) <- c(2008:2013)
  for (yr in (c(2008:2013))) {
    print(yr)
    yrs[names(yrs)==yr] <- round(nrow(pdata[pdata$YEAR==yr,])*.3 + nrow(pdata[pdata$YEAR==yr,]))
      
    #take random sample from adata based on number of records in pdata, then clip to NA (study area poly)
    adata.sub <- adata[adata$YEAR == yr,]
    adata.samp <- adata.sub[sample(nrow(adata.sub), yrs[names(yrs)==yr]),]
    
    #Convert to spatial points df
    coordinates(adata.samp) <- ~LONGITUDE+LATITUDE
    
    #clip adata.sp to sa boundaries (north america) - reduce the number of points in analysis
    adata.na <- over(adata.samp, sa)
    adata.na2 <- cbind(adata.samp, adata.na)
    adata.na2 <- adata.na2[!is.na(adata.na2$ras),]
    
    #gets rid of duplicate records that are part of the same group. 
    #Takes only the first record to move to analysis, keeps all single records that are not part of a group.
    #Do this last bc I want to run this step on a smaller set of data bc it takes a long time.
    adata.na2 <- GroupDuplicates(adata.na2)
    
    out_adata <- paste0(dirname(afile), "/",unlist(strsplit(basename(afile), "\\."))[1], "_", yr, "_sub.csv")
    write.csv(adata.na2, file=out_adata, quote=F, row.names=F)
   
    
    }#end yr loop
}#end spp loop
  
  
