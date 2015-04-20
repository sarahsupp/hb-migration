# From: https://stat.ethz.ch/pipermail/r-sig-geo/2012-March/014409.html

# convert alpha hull (from {alphahull} package) to spatial polygon data frame
hb.ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
  require(alphahull)
  require(maptools)
  if (class(x) != "ahull"){
    stop("x needs to be an ahull class object")
  }
  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)
  # Remove all cases where the coordinates are all the same      
  xdf <- subset(xdf,xdf$r > 0)
  res <- NULL
  if (nrow(xdf) > 0){
    # Convert each arc to a line segment
    linesj <- list()
    prevx<-NULL
    prevy<-NULL
    j<-1
    for(i in 1:nrow(xdf)){
      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)
      # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2),0)
      # Calculate coordinates from arc() description for ipoints along the arc.
      angles <- anglesArc(v, theta)
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)
      # Check for line segments that should be joined up and combine their coordinates
      if (is.null(prevx)){
        prevx<-x
        prevy<-y
      } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
          if (i == nrow(xdf)){
          #We have got to the end of the dataset
            prevx<-append(prevx,x[2:ipoints])
            prevy<-append(prevy,y[2:ipoints])
            prevx[length(prevx)]<-prevx[1]
            prevy[length(prevy)]<-prevy[1]
            coordsj<-cbind(prevx,prevy)
            colnames(coordsj)<-NULL
            # Build as Line and then Lines class
            linej <- Line(coordsj)
            linesj[[j]] <- Lines(linej, ID = as.character(j))
          } else {
            prevx<-append(prevx,x[2:ipoints])
            prevy<-append(prevy,y[2:ipoints])
          }
        } else {
      # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
          prevx[length(prevx)]<-prevx[1]
          prevy[length(prevy)]<-prevy[1]
          coordsj<-cbind(prevx,prevy)
          colnames(coordsj)<-NULL
      # Build as Line and then Lines class
          linej <- Line(coordsj)
          linesj[[j]] <- Lines(linej, ID = as.character(j))
          j<-j+1
          prevx<-NULL
          prevy<-NULL
        }
    }
    # Promote to SpatialLines
    lspl <- SpatialLines(linesj)
    # Convert lines to polygons
    # Pull out Lines slot and check which lines have start and end points that are the same
    lns <- slot(lspl, "lines")
    polys <- sapply(lns, function(x) { 
      crds <- slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    }) 
    # Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
    # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }  
  return(res)
}

# convert date (month, day, year) to day of year (doy)
hb.doy <- function(x,year){
  # convert DOY to Date
  as.Date(x-1,origin=as.Date(paste0(year, "-01-01")))
}

# calculate the largest polygon in a spatial polygon data frame
hb.largest.poly <- function(x){
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

# get dates for a three week window around the current day (x)
hb.get.3.week <- function(x){
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

# convert factor to numeric
as.num <- function(x){
  as.numeric(as.character(x))
}

# convert columns in a dataframe to character or numeric types based on input string
df.col.class <- function(x,y){
  if(ncol(x) != length(y)){
    print(paste('length of arguements not equal: x=', ncol(x), ' y=', length(y)))
  }
  else{
    for(i in 1:ncol(x)){
      x[,i] <- as.character(x[,i])
      if(y[i] == 'n'){
        x[,i] <- as.numeric(x[,i])
      }
    }
  }
  x
}
