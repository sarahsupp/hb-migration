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

################################### I N P U T ###################################
# number of days in the time window
period <- 5
#################################################################################


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
for(species.code in c('bchu','bthu', 'cahu','ruhu')){ 
  
  meta.df <- data.frame(species=as.character(), year=as.numeric(), w=as.numeric(), alpha=as.numeric(), in.hull=as.numeric(), outlier=as.numeric(), total=as.numeric())
  
  # 0. SET UP ENVIRONMENT
  # ---------------------
  
  # set path for presence data and load presence data into dataframe
  prs.path <- paste(path.prefix, 'Share/kguay/hummingbird/data/ebird/presence/aggregate_by_species/t', period, '/', sep='')
  prs <- as.data.frame(sapply(ldply(paste(prs.path,list.files(prs.path, pattern=paste('.*',species.code,'_humdat_west.txt$',sep='')), sep=''), read.table, sep=',', quote = "\"", header=T) , function(x) gsub('\\\\', '', x)))
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
    print(year)
    # subset daraframes by year
    prs.year <- subset(prs, prs$YEAR==year)
    abs.year <- subset(abs, abs$YEAR==year)
    
    # loop through 1:N and calculate alpha hulls. Save to list to use in step 3.
    for(w in sort.int(unique(as.num(prs.year$window)))){
      print(w)
      # subset presence dataframe by window
      prs.window <- subset(prs.year, prs.year$window==w)
      
      print(paste(w,nrow(prs.window)))
      
      # get the absense points that are from the same days as the presence window
      abs.window <- subset(abs.year, abs.year$DAY %in% unique(prs.window$DAY))
      
      
      prs.window.coords <- as.data.frame(cbind(LONGITUDE=prs.window$LONGITUDE, LATITUDE=prs.window$LATITUDE))
      prs.window.info <- as.data.frame(cbind(YEAR=prs.window$YEAR, DAY=prs.window$DAY, TIME=prs.window$TIME))
      
      abs.window.coords <- as.data.frame(cbind(LONGITUDE=abs.window$LONGITUDE, LATITUDE=abs.window$LATITUDE))
      abs.window.info   <- as.data.frame(cbind(YEAR=abs.window$YEAR, DAY=abs.window$DAY, TIME=as.character(abs.window$TIME)))
      
      # make sure that there are enough observations to make a hull (> 3)
      if(nrow(prs.window) > 15){
        
        # spdf <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
        
        # Itterate the value of alpha until it produces x% outliers
        # ---------------------------------------------------------
        # determine the ideal number of outliers
        outlier.goal <- nrow(prs.window.coords) %/% 10
        inside.goal <- nrow(prs.window.coords) - outlier.goal
        
        a <- 1.3
        go <- T # keep going until go == FALSE
        a.increment <- 0.1
        rep.last <- 0
        history <- c(seq(0,1,1/30))
        history.a <- c(seq(0,1,1/30))
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
          # sometimes alpha is too low to include any of the points. In this case, we need to increase a until a hull exists
          while(hull.shp$Area < 0.01) {
            a <- a + 1
            hull <- ahull(unique(prs.window.coords), alpha=a)
            hull.shp <- hb.ah2sp(hull)
          }
          
          # extract points that fall within hull
          hull.pts <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
          
          hull.pts <- tryCatch({ hull.pts[!is.na(over(hull.pts,as(hull.shp,"SpatialPolygons"))),] },error=function(cond){return(NA)},warning=function(cond){return(NA)})
          
          if(is.na(hull.pts)){
            go <- F
            # make alpha equal to the alpha that produced the maximum in history
            a <- history.a[2]
            
            # calculate alpha hull for those observations
            hull <- ahull(unique(prs.window.coords), alpha=a)
            
            # convert alpha hull to shapefile (using function in functions.R)
            hull.shp <- hb.ah2sp(hull)
            
            # extract points that fall within hull
            hull.pts <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
            
            hull.pts <- hull.pts[!is.na(over(hull.pts,as(hull.shp,"SpatialPolygons"))),]
            
          }
          
          # in.hull <- apply(prs.window.coords, 1, iah, hull)
          # hull.pts.2 <- subset(hull.pts, in.hull)
          
          if( ( (history[1] == history[3]) & (history[1] != history[2]) & (history[2] != history[3]) ) | (length(unique(history[1:30]))==1) | (nrow(hull.pts) < history[1]) ){
            a.increment <- 0.01
            rep.last <- rep.last + 1
            if(rep.last > 5 | (nrow(hull.pts) < history[1])){
              if((abs(history[1] - inside.goal) < abs(history[2] - inside.goal))  ){
                a <- history.a[1]
              }
              else{
                a <- history.a[2]
              }
              go <- F
            }
          }
          
          # if the number of points in the hull does not equal the inside.goal, adjust alpha
          if((nrow(hull.pts) > (inside.goal+1)) & go){
            a <- a - a.increment
          }else if((nrow(hull.pts) < (inside.goal-1) )  & go){
            a <- a + a.increment
          }else{
            go <- F
          }
          
          history <- c(nrow(hull.pts), history)
          history.a <- c(a, history.a)
          
          #print(paste('  ',inside.goal,history[1],a, a.increment, rep.last, '\t', history[2], history[3], history[4], history[20], (length(unique(history[1:30]))==1 & !is.na(history[1])), nrow(hull.pts) < history[1], sep='  '))
        } # end while go
        
      } # end if
      else{
        pts.orig <- unique(prs.window.coords)
        pts <- pts.orig
        #x <- abs(pts[1,1] - pts[2,1])
        #y <- abs(pts[1,2] - pts[2,2])
        for(i in 1:20){
          pts <- rbind(pts, pts.orig + runif(nrow(pts.orig)*2,-1, 1))
        }
        
        a <- 5
        
        # calculate alpha hull for those observations
        hull <- ahull(unique(pts), alpha=a)
        
        # convert alpha hull to shapefile (using function in functions.R)
        hull.shp <- hb.ah2sp(hull)
        
        # extract points that fall within hull
        hull.pts <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
        
        hull.pts <- hull.pts[!is.na(over(hull.pts,as(hull.shp,"SpatialPolygons"))),]
        
      }
      # ABSENSE
      # -------
      # extract absense points from hull
      abs.hull <- SpatialPointsDataFrame(abs.window.coords, abs.window)
      abs.hull <- abs.hull[!is.na(over(abs.hull,as(hull.shp,"SpatialPolygons"))),]
      
      # select absense points that correspond with presence observations (day)
      for(day in unique(prs.window$DAY)){
        if(length(abs.hull) != 0){
          npoints <- nrow(subset(prs.window, prs.window$DAY == day))
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
                                         location.long=as.num(hull.pts@coords[,1]), location.lat=as.num(hull.pts@coords[,2]), height.above.ellipsoid='' ))
      
      prs.df.m3 <- rbind(prs.df.m3, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY)-3, origin = paste(year, "-01-01", sep='')), ' ', as.character(hull.pts$TIME), '.000', sep=''), 
                                               location.long=as.num(hull.pts@coords[,1]), location.lat=as.num(hull.pts@coords[,2]), height.above.ellipsoid='' ))
      
      prs.df.p3 <- rbind(prs.df.p3, data.frame(timestamp=paste(as.Date(as.num(hull.pts$DAY)+3, origin = paste(year, "-01-01", sep='')), ' ', as.character(hull.pts$TIME), '.000', sep=''), 
                                               location.long=as.num(hull.pts@coords[,1]), location.lat=as.num(hull.pts@coords[,2]), height.above.ellipsoid='' ))
      
      
      total <- SpatialPointsDataFrame(prs.window.coords, prs.window.info)
      outl <- total[is.na(over(total,as(hull.shp,"SpatialPolygons"))),]    
      inh <- total[!is.na(over(total,as(hull.shp,"SpatialPolygons"))),]        
      
      
      meta.df <- rbind(meta.df, data.frame(species=species.code, year=year, w=w, alpha=a, in.hull=length(hull.pts), outlier=length(outl), total=nrow(unique(prs.window.coords))))
      
#       png(paste('~/Desktop/hb_hull_plots/hull_plot_',w,'.png', sep=''))
#       
#       plot(newmap, xlim = c(-115, -108), ylim = c(30,55), asp = 1)
#       plot(hull.shp, col='lightyellow', lwd=1, add=T)
#       title(main=w)
#       plot(outl, pch=4, add=T)
#       plot(hull.pts, col='darkgreen', pch=18, add=T)
#       
#       text(-127,35,paste('window:',w), pos=4) 
#       text(-127,34,paste('alpha:',round(a, digits=3)), pos=4) 
#       text(-127,33,paste('in hull:',length(inh)), pos=4)
#       text(-127,32,paste('outlier:',length(outl)), pos=4)
#       text(-127,31,paste('total:',length(total)), pos=4)   
#       text(-127,30,paste('outliers:',round(length(outl)/(length(total))*100, digits=0), '%' ), pos=4)               
#       
#       dev.off()
    } # end window
  } # end year
  
  # 3. WRITE MOVEBANK SUBMISSION FILE (CSV)
  
  write.df(prs.df, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit_csv/', species.code, '_prs_movebank_submit_fcl', '.csv', sep=''))
  write.df(prs.df.m3, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit_csv/', species.code, '_prs_movebank_submit_min_3.csv', sep=''))
  write.df(prs.df.p3, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit_csv/', species.code, '_prs_movebank_submit_pls_3.csv', sep=''))
  write.df(abs.df, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit_csv/', species.code, '_abs_movebank_submit_fcl.csv', sep=''))
  write.csv(meta.df, paste(path.prefix, 'Share/kguay/hummingbird/movebank/submit_csv/meta/', species.code, '_meta.csv', sep=''), quote=F, row.names=F)

} # end species

# 
# require(geosphere)
# 
# data.df <- as.data.frame(hull.pts@coords)
# coords <- hull.pts@coords
# # Make a list, each element of the list has one SpatialPolygons polygon (one circle)
# # Here, coords is a data.frame with the lon and lat of each point.
# list.poly=list()
# for (i in 1:length(coords[,1])) {
#   myCircle = data.frame(matrix(NA,361,2))
#   for (angle in 0:360) { myCircle[(angle+1),] = destPoint(coords[i,],angle,100000) }
#   list.poly[i] = SpatialPolygons(list(Polygons(list(Polygon(SpatialPoints(myCircle))),1))) }
# 
# #----------------------------
# # Get the ID names of each polygon...intially they are all the same, but each polygon needs a unique ID:
# IDs <- sapply(list.poly, function(x) slot(slot(x, "polygons")[[1]], "ID"))
# # Give each polygon a unique ID  (1: number of list elements)
# Spol1 <- SpatialPolygons(lapply(1:length(list.poly), function(i) {
#   Pol <- slot(list.poly[[i]], "polygons")[[1]]
#   slot(Pol, "ID") <- as.character(i)
#   Pol } ) )
# # Get the IDs, and set the row.names of the polygons as these IDs
# IDs = sapply(slot(Spol1, "polygons"), function(x) slot(x, "ID"))
# row.names(Spol1) = IDs
# #As a check, this should be true now:
# length(unique(IDs)) == length(list.poly)
# 
# 
# #----------------------------
# # Build a data.frame to attach to the SpatialPolyons
# data.df = data.frame(IDs);  row.names(data.df)=IDs
# spol.df = SpatialPolygonsDataFrame(Spol1,data = data.df)
# projection(spol.df) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# #assign the final variable, "polys.wgs84.35m" as the spatial polygons dataframe:
# polys.wgs84.35m=spol.df
# 
# unionSpatialPolygons(spol.df, IDs)
# 
# 




