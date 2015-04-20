# Compute weekly alpha hulls from ebird data
# The main code will be within a function s.t. it can be adapted for each species

# Kevin Guay (kguay@whrc.org)
# on 13 Apr 2015
# mo 13 Apr 2015

# load required packages
require(plyr)
require(rworldmap)
require(alphahull)

# set base path depending on OS (Windows/ Unix)
if(Sys.info()['sysname'] == 'Windows'){
  path.prefix <- 'C:/'
}else{
  path.prefix <- '/mnt/arctic/c/'
}

# load functions from hb_alpha_function.R
source(paste(path.prefix, 'Share/kguay/hummingbird/scripts/hb_alpha_functions.R', sep=''))

# 0. SET UP ENVIRONMENT
# ---------------------

# set species code
species <- c('bchu','bthu','cahu','rthu')
species.code <- species[1]

# set path for presence data and load presence data into dataframe
prs.path <- paste(path.prefix, 'Share/tcormier/hummingbirds/migration_study/data/ebird/present_points/updated_20150307/eBird_checklists_2008-2014/aggregate_by_species/', sep='')
prs <- as.data.frame(sapply(ldply(paste(prs.path,list.files(prs.path, pattern=paste('t3.*',species.code,'_humdat_west.txt$',sep='')), sep=''), read.table, sep=',', quote = "\"", header=T) , function(x) gsub('\\\\', '', x)))
# convert column classes from factor to character or numeric
prs <- df.col.class(prs, c('c','c','n','n','c','c','c','c','c','c','c','c','n','n','c','c','n','n','n','n','n','n','c'))

# set path for absence data and load absence data into dataframe
abs.path <- paste(path.prefix, 'Share/tcormier/hummingbirds/migration_study/data/ebird/absent_points/migrants/updated_20150307/', sep='')
abs <- as.data.frame(ldply(paste(abs.path,list.files(abs.path, pattern=paste(species.code,'_humdat_abs.txt$',sep='')), sep=''), read.table, sep='\t', header=T))

# read file that determines start/end dates of migration for each species
eb.dates <- read.table(paste(path.prefix, 'Share/tcormier/hummingbirds/migration_study/data/ebird/present_points/updated_20150307/eBird_checklists_2008-2014/aggregate_by_species/', species.code, '_migration_west.txt', sep=''), header=T)

# init a map
newmap <- getMap(resolution = "low")

# 1. CALCULATE ALPHA HULLS
# ------------------------

# initialize dataframe for the output (move-bank csv)
prs.df <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
prs.df.m3 <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
prs.df.p3 <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
abs.df <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())

# itterate through years
for(year in unique(prs$YEAR)){
  
  # subset daraframes by year
  prs.year <- subset(prs, prs$YEAR==year)
  abs.year <- subset(abs, abs$YEAR==year)
  
  # loop through 1:N and calculate alpha hulls. Save to list to use in step 3.
  for(w in unique(prs.year$window)){
    
    # subset presence dataframe by window
    prs.window <- subset(prs.year, prs.year$window==w)
    # get the absense points that are from the same days as the presence window
    abs.window <- subset(abs.year, abs.year$DAY %in% unique(prs.window$DAY))
    
    # make sure that there are enough observations to make a hull (> 3)
    if(nrow(prs.window) > 5){
      
      # calculate alpha hull for points in prs.window
      
      prs.window.coords <- as.data.frame(cbind(LONGITUDE=prs.window$LONGITUDE, LATITUDE=prs.window$LATITUDE))
      prs.window.info <- as.data.frame(cbind(YEAR=prs.window$YEAR, DAY=prs.window$DAY, TIME=prs.window$TIME))
      
      abs.window.coords <- as.data.frame(cbind(LONGITUDE=abs.window$LONGITUDE, LATITUDE=abs.window$LATITUDE))
      abs.window.info   <- as.data.frame(cbind(YEAR=abs.window$YEAR, DAY=abs.window$DAY, TIME=as.character(abs.window$TIME)))
      # spdf <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
      
      # Itterate the value of alpha until it produces x% outliers
      # ---------------------------------------------------------
      # determine the ideal number of outliers
      outlier.goal <- nrow(prs.window.coords) %/% 5
      inside.goal <- nrow(prs.window.coords) - outlier.goal
      
      a <- 1.3
      go <- T # keep going until go == FALSE
      last <- 10
      last.last <- 10
      last.a <- 1
      last.last.a <- 1
      a.increment <- 0.1
      while(go){        
        # calculate alpha hull for those observations
        hull <- ahull(unique(prs.window.coords), alpha=a)
        
        # convert alpha hull to shapefile (using function in functions.R)
        hull.shp <- hb.ah2sp(hull)
        
        # sometimes alpha is too low to include any of the points. In this case, we need to increase a until a hull exists
        while(is.null(hull.shp)) {
          a <- a + 1
          hull <- ahull(unique(prs.window.coords), alpha=a)
          hull.shp <- hb.ah2sp(hull)
        }
        
        # extract points that fall within hull
        hull.pts <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
        hull.pts <- hull.pts[!is.na(over(hull.pts,as(hull.shp,"SpatialPolygons"))),]
        
        
        if(((nrow(hull.pts) == last.last) | (nrow(hull.pts) == last)) & (a.increment == 0.1)){
          a.increment <- 0.01
        } else if((nrow(hull.pts) == last.last) | (nrow(hull.pts) == last)){
          if(abs(last - inside.goal) < abs(last.last - inside.goal)){
            a <- last.a
          }
          else{
            a <- last.last.a
          }
          go <- F
        }
        
        # if the number of points in the hull does not equal the inside.goal, adjust alpha
        if((nrow(hull.pts) > (inside.goal+1)) & go){
          a <- a - a.increment
        }else if(nrow(hull.pts) < (inside.goal-1)){
          a <- a + a.increment
        }else{
          go <- F
        }
        
        last.last <- last
        last <- nrow(hull.pts)
        
        last.last.a <- last.a
        last.a <- a
        
        # print(paste('  ',inside.goal,last,a, a.increment, sep='  '))
      }
      
      # ABSENSE
      # -------
      # extract absense points from hull
      abs.hull <- SpatialPointsDataFrame(abs.window.coords, abs.window)
      abs.hull <- abs.hull[!is.na(over(abs.hull,as(hull.shp,"SpatialPolygons"))),]
      
      # select absense points that correspond with presence observations (day)
      for(day in unique(prs.window$DAY)){
        npoints <- length(subset(prs.window, prs.window$DAY == day))
        if(npoints < length(abs.hull)){
          abs.hull.smpl <- abs.hull[sample(1:length(abs.hull),npoints),]
        } else{
          abs.hull.smpl <- abs.hull
        }
        abs.df <- rbind(abs.df, data.frame(timestamp=paste(as.Date(abs.hull.smpl$DAY, origin = paste(year, "-01-01", sep='')), ' 12:00:00.000', sep=''), 
                                           location.long=as.num(abs.hull.smpl$LONGITUDE), location.lat=as.num(abs.hull.smpl$LATITUDE), height.above.ellipsoid='' ))
      }
      
      
      prs.df <- rbind(prs.df, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY), origin = paste(year, "-01-01", sep='')), ' 12:00:00.000', sep=''), 
                                         location.long=as.num(hull.pts$LONGITUDE), location.lat=as.num(hull.pts$LATITUDE), height.above.ellipsoid='' ))
      
      prs.df.m3 <- rbind(prs.df.m3, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY)-3, origin = paste(year, "-01-01", sep='')), ' 12:00:00.000', sep=''), 
                                               location.long=as.num(hull.pts$LONGITUDE), location.lat=as.num(hull.pts$LATITUDE), height.above.ellipsoid='' ))
      
      prs.df.p3 <- rbind(prs.df.p3, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY)+3, origin = paste(year, "-01-01", sep='')), ' 12:00:00.000', sep=''), 
                                               location.long=as.num(hull.pts$LONGITUDE), location.lat=as.num(hull.pts$LATITUDE), height.above.ellipsoid='' ))
      
      # plot(newmap, xlim = c(-110, -103), ylim = c(20, 60), asp = 1)
      # plot(hull.shp, add=T)
      # plot(hull.pts, add=T)
      # plot(abs.hull, add=T)
      
    }
  }
}

# 3. WRITE MOVEBANK SUBMISSION FILE (CSV)
# change column names of dataframe
names(df) <- c("timestamp","location-long","location-lat","height-above-ellipsoid")
# sort datafrmae by longitude
df <- df[order(df[,2]),]
write.csv(df, paste(path.prefix, 'Share/tcormier/hummingbirds/movebank/submit/test/hb_mb_', species.code, '.csv', sep=''), quote=T, row.names=F)


# EXTRAS
# calculate the number of polygons in a spatialPolygonDataFrame
# npoly <- length(hull.shp@polygons[[1]]@plotOrder)

