# Compute weekly alpha hulls from ebird data
# The main code will be within a function s.t. it can be adapted for each species

# Kevin Guay (kguay@whrc.org)
# on 13 Apr 2015
# mo 20 Apr 2015

# log
# version   m   d   y     notes
# -------   --  --  --    -----
# 1.0.0     4   20  15    working copy. submitted to repo

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

# function to write the dataframe to a file (sets column names)
write.df <- function(x, filename){
  # change column names of dataframe
  names(x) <- c("timestamp","location-long","location-lat","height-above-ellipsoid")
  write.csv(x, filename, quote=F, row.names=F)
}

# init a map
newmap <- getMap(resolution = "low")

# ------------------------------------------------------------------------------

# loop through species
for(species.code in c('bchu','bthu','cahu','ruhu')){ 
  
  # 0. SET UP ENVIRONMENT
  # ---------------------
  
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
  
  # 1. CALCULATE ALPHA HULLS
  # ------------------------
  
  # initialize dataframe for the output (move-bank csv)
  prs.df <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
  prs.df.m3 <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
  prs.df.p3 <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
  abs.df <- data.frame(timestamp=as.character(), location.long=as.numeric(), location.lat=as.numeric(), height.above.ellipsoid=as.character())
  
  # loop through years
  for(year in unique(prs$YEAR)){
    
    # subset daraframes by year
    prs.year <- subset(prs, prs$YEAR==year)
    abs.year <- subset(abs, abs$YEAR==year)
    
    # loop through 1:N and calculate alpha hulls. Save to list to use in step 3.
    for(w in sort.int(unique(prs.year$window))){
      
      # subset presence dataframe by window
      prs.window <- subset(prs.year, prs.year$window==w)
      # get the absense points that are from the same days as the presence window
      abs.window <- subset(abs.year, abs.year$DAY %in% unique(prs.window$DAY))
      
      # make sure that there are enough observations to make a hull (> 3)
      if(nrow(prs.window) > 10){
        
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
        } # end while go
        
        # ABSENSE
        # -------
        # extract absense points from hull
        abs.hull <- SpatialPointsDataFrame(abs.window.coords, abs.window)
        abs.hull <- abs.hull[!is.na(over(abs.hull,as(hull.shp,"SpatialPolygons"))),]
        
        # select absense points that correspond with presence observations (day)
        for(day in unique(prs.window$DAY)){
          if(length(abs.hull) != 0){
            npoints <- length(subset(prs.window, prs.window$DAY == day))
            if(npoints < length(abs.hull)){
              abs.hull.smpl <- abs.hull[sample(1:length(abs.hull),npoints),]
            } else{
              abs.hull.smpl <- abs.hull
            }
            abs.df <- rbind(abs.df, data.frame(timestamp=paste(as.Date(abs.hull.smpl$DAY, origin = paste(year, "-01-01", sep='')), ' ', as.character(abs.hull.smpl$TIME), '.000', sep=''), 
                                               location.long=as.num(abs.hull.smpl$LONGITUDE), location.lat=as.num(abs.hull.smpl$LATITUDE), height.above.ellipsoid='' ))
          }
        } # end for day
        
        # append each set of points to the corresponding data frame
        prs.df <- rbind(prs.df, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY), origin = paste(year, "-01-01", sep='')), ' ', as.character(hull.pts$TIME), '.000', sep=''), 
                                           location.long=as.num(hull.pts$LONGITUDE), location.lat=as.num(hull.pts$LATITUDE), height.above.ellipsoid='' ))
        
        prs.df.m3 <- rbind(prs.df.m3, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY)-3, origin = paste(year, "-01-01", sep='')), ' ', as.character(hull.pts$TIME), '.000', sep=''), 
                                                 location.long=as.num(hull.pts$LONGITUDE), location.lat=as.num(hull.pts$LATITUDE), height.above.ellipsoid='' ))
        
        prs.df.p3 <- rbind(prs.df.p3, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY)+3, origin = paste(year, "-01-01", sep='')), ' ', as.character(hull.pts$TIME), '.000', sep=''), 
                                                 location.long=as.num(hull.pts$LONGITUDE), location.lat=as.num(hull.pts$LATITUDE), height.above.ellipsoid='' ))
        
        # plot(newmap, xlim = c(-110, -103), ylim = c(20, 60), asp = 1)
        # plot(hull.shp, add=T)
        # plot(hull.pts, add=T)
        # plot(abs.hull, add=T)
        
      } # end if
    } # end window
  } # end year
  
  # 3. WRITE MOVEBANK SUBMISSION FILE (CSV)
  
  write.df(prs.df, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit/', species.code, '_prs_movebank_submit_fcl', '.csv', sep=''))
  write.df(prs.df.m3, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit/', species.code, '_prs_movebank_submit_min_3.csv', sep=''))
  write.df(prs.df.p3, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit/', species.code, '_prs_movebank_submit_pls_3.csv', sep=''))
  write.df(abs.df, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit/', species.code, '_abs_movebank_submit_fcl.csv', sep=''))
  
} # end species








