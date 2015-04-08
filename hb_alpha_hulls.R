# Compute weekly alpha hulls from ebird data
# The main code will be within a function s.t. it can be adapted for each species

# Kevin Guay (kguay@whrc.org)
# on 25 Feb 2015
# mo 06 Apr 2015

require(plyr)
require(rworldmap)

# set base path depending on OS
if(Sys.info()['sysname'] == 'Windows'){
  path.prefix <- 'C:/'
}else{
  path.prefix <- '/mnt/arctic/c/'
}

source(paste(path.prefix,'Share/tcormier/hummingbirds/scripts/R/functions.R', sep=''))

doy <- function(x,year){
  # convert DOY to Date
  as.Date(x-1,origin=as.Date(paste0(year, "-01-01")))
}

largest.poly <- function(x){
  # taken from: http://gis.stackexchange.com/questions/119624/extract-areas-of-multi-part-polygons-spatialpolygonsdataframe-r
  
  # total area
  sapply(slot(x, "polygons"), slot, "area")
  
  # get list of individual polys
  p <- lapply(x@polygons , slot , "Polygons")
  
  # areas of individual polygons
  size <- unlist(lapply(p[[1]], function(x) slot(x, "area")))
  # get the index of thepolygon with the largest area
  index.of.largest <- match(size[order(-size)][c(1)], size)
  
  SpatialPolygons(list(Polygons(list(Polygon(SpatialPoints(p[[1]][[index.of.largest]]))),1)))
}

get.3.week <- function(x){
  df <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
  
  # loop through 10 days before and 10 days after current date.
  for(i in -10:10){
    df <- rbind(df, data.frame(timestamp=paste(as.Date(x$DAY + i - 1, origin = paste(x$YEAR, "-01-01", sep='')), ' ', x$TIME, '.000', sep=''), 
                               location.long=x$LONGITUDE, 
                               location.lat=x$LATITUDE,
                               height.above.ellipsoid=''
    ))
  }
  # return dataframe
  df
}

humhul <- function(x, plot=F){
  require(alphahull)
  
  species.code <- x
  
  eb.prs.path <- paste(path.prefix, 'Share/tcormier/hummingbirds/migration_study/data/ebird/present_points/updated_20150307/eBird_checklists_2008-2014/aggregate_by_species/', sep='')
  eb.abs.path <- paste(path.prefix, 'Share/tcormier/hummingbirds/migration_study/data/ebird/absent_points/migrants/updated_20150307/', sep='')
  
  # 1. read in ebird prs data (from text or CSV file)
  eb.prs <- as.data.frame(sapply(   ldply(paste(eb.prs.path,list.files(eb.prs.path, pattern=paste('.*',species.code,'_humdat_west.txt$',sep='')), sep=''), read.table, header=T, stringsAsFactors=F)    , function(x) gsub("\"", "", x)))
  
  # read in ebird abs data
  eb.abs <- as.data.frame(sapply(   ldply(paste(eb.abs.path,list.files(eb.abs.path, pattern=paste('.*',species.code,'_humdat_abs.txt',sep='')), sep=''), read.table, header=T, stringsAsFactors=F)    , function(x) gsub("\"", "", x)))
  
  eb.abs <- ldply(paste(eb.abs.path,list.files(eb.abs.path, pattern=paste('.*',species.code,'_humdat_abs.txt$',sep='')), sep=''), read.table, header=T, stringsAsFactors=F)
  
  # add week column
  eb.prs <- cbind(eb.prs, data.frame(WEEK=as.numeric( format(doy(eb.prs$DAY,eb.prs$YEAR), "%U")) ))
  eb.abs <- cbind(eb.abs, data.frame(WEEK=as.numeric( format(doy(eb.abs$DAY,eb.abs$YEAR), "%U")) ))
  
  # look up start/ end date for spring/ fall migration
  # path: /mnt/arctic/c/Share/tcormier/hummingbirds/migration_study/data/supp_migration
  eb.dates <- read.delim(paste(path.prefix, 'Share/tcormier/hummingbirds/migration_study/data/ebird/present_points/updated_20150307/eBird_checklists_2008-2014/aggregate_by_species/', species.code, '_migration_west.txt', sep=''))
  
  if(plot){
    newmap <- getMap(resolution = "low")
    plot(newmap, xlim = c(-110, -103), ylim = c(20, 60), asp = 1)
  }
  
  for(ab in c('eb.prs.spr','eb.abs.spr')){
    assign()
  
  for(year in sort(unique(eb.pa.spr$YEAR))){
    # subset the dataframe based on start/ end dates for migration season
    eb.prs.spr <- subset(eb.prs, (eb.prs$DAY >= subset(eb.dates, year==y)$spr_begin & eb.prs$DAY < subset(eb.dates, year==year)$spr_end))
    eb.prs.fal <- subset(eb.prs, (eb.prs$DAY >= subset(eb.dates, year==y)$spr_end & eb.prs$DAY < subset(eb.dates, year==year)$fal_end))
    
    eb.abs.spr <- subset(eb.abs, (eb.abs$DAY >= subset(eb.dates, year==y)$spr_begin & eb.abs$DAY < subset(eb.dates, year==year)$spr_end))
    eb.abs.fal <- subset(eb.abs, (eb.abs$DAY >= subset(eb.dates, year==y)$spr_end & eb.abs$DAY < subset(eb.dates, year==year)$fal_end))
    
    # for each week:
    for(week in sort(unique(eb.pa.spr$WEEK))){
      print(paste(species, week))
      
      # 2.a. subset ebird observations by week.
      eb.pa.week <- subset(eb.pa.spr, (eb.pa.spr$SCI_NAME == species & eb.pa.spr$WEEK==week & eb.pa.spr$YEAR==year & eb.pa.spr$LONGITUDE <= -103))
      
      if(nrow(eb.prs.week) > 10){
        # make matrix with lat lon
        eb.pa.week.coords <- matrix(c(eb.pa.week$LONGITUDE, eb.pa.week$LATITUDE), nc=2)
        eb.pa.week.time <- matrix(c(eb.pa.week$YEAR, eb.pa.week$DAY, eb.pa.week$TIME), nc=3)
        
        eb.pa.week.coords <- as.data.frame(cbind(LONGITUDE=eb.pa.week$LONGITUDE, LATITUDE=eb.pa.week$LATITUDE))
        eb.pa.week.info <- as.data.frame(cbind(YEAR=eb.pa.week$YEAR, DAY=eb.pa.week$DAY, TIME=eb.pa.week$TIME))
        spdf <- SpatialPointsDataFrame(eb.prs.week.coords, eb.prs.week)
        
        # calculate alpha hull for those observations
        eb.pa.week.hull <- ahull(unique(eb.pa.week.coords), alpha=0.8)
        
        # convert alpha hull to shapefile (using function in functions.R)
        eb.pa.week.hull.shp <- ah2sp(eb.pa.week.hull)
        
        # calculate the largest alpha hull (optional)
        #eb.pa.week.hull.shp.largest <- largest.poly(eb.pa.week.hull.shp)
        
        # subset points that are in hull
        eb.pa.week.hull.pts <- SpatialPointsDataFrame(eb.pa.week.coords, eb.pa.week)
        eb.pa.week.hull.pts <- eb.pa.week.hull.pts[!is.na(over(eb.pa.week.hull.pts,as(eb.pa.week.hull.shp,"SpatialPolygons"))),]
        
        df <- rbind(df, get.3.week(eb.pa.week.hull.pts))
        
        if(plot){
          plot(eb.pa.week.hull.shp, fill=NULL, border=week, lwd=3, add=T)
          plot(eb.pa.week.hull.pts, col=week, pch=19, add=T)
        }
      } # end if
      
    } # end for week
  } # end for year
  
  # 3. WRITE MOVEBANK SUBMISSION FILE (CSV)
  # change column names of dataframe
  names(df) <- c("timestamp","location-long","location-lat","height-above-ellipsoid")
  # sort datafrmae by longitude
  df <- df[order(df[,2]),]
  species.name <- gsub(' ', '_', as.character(species))
  write.csv(df, paste(path.prefix, 'Share/tcormier/hummingbirds/movebank/submit/', output.folder, 'hb_mb_', year, '_', species.name, '.csv', sep=''), quote=F, row.names=F)
  
  
}


species <- c('bchu','bthu','cahu','rthu')

humhul('bchu')








