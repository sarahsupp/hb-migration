# Define funcitons to use in hb-migration code

DateConvert = function(date){
  #convert eBird column OBSERVATION.DATE into year, month and day columns
  
  year <-sapply(date,function(x){
    as.numeric(substring(x, 1, 4))
  })
  
  month <- sapply(date,function(x){
    as.numeric(substring(x, 6, 7))
  }) 
  
  day <- sapply(date,function(x){
    as.numeric(substring(x, 9, 10))
  })
  
  newdate = as.data.frame(cbind(year,month,day), row.names=FALSE)

  newdate = JulianDay(newdate)

  return (newdate)
}

JulianDay = function(daymoyr){
  ## Add a column for julian day
  require(chron)
  for (row in 1:nrow(daymoyr)){
    line = daymoyr[row,]
    #Julian date ranges 1:365
    daymoyr$julian[row] = julian(line$month, line$day, line$year, origin. = c(month=1, day=1, year=line$year)) + 1
  }
  return(daymoyr)
}


PlotRecords = function(yeardata, species) {
  #plots the number of records across the years
  
  t = as.data.frame(table(yeardata))
  names(t) = c('year', 'count')
 
  barplot = ggplot(data = t, aes(year, count)) + geom_bar(fill="cadetblue2") + theme_bw() + 
    ggtitle(species)
  
  print(barplot)
  
  return(t)
}


MeanDailyLoc = function(yeardata, species) {
  # input is yearly data frame for a species
  # makes a new dataframe with mean julian day lat and long (and sd)
  require(Rmisc)
  require(chron)
  
  #count the number of days in the year
  year = dat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  julian = seq(1:numdays)
  
  df = data.frame("species"=NA, "jday"=1, "month"=1, "count"=1, "meanlat"=1, "sdlat"=1, "meanlon"=1, 
                  "sdlon"=1, "ucilat"=1, "lcilat"=1, "minlat"=1, "maxlat"=1)
  outcount = 1
  
  for (d in 1:length(julian)) {
    daydat = yeardata[which(yeardata$julian == julian[d]),]
    count = nrow(daydat)
    month = as.numeric(months(julian[d]))
    lon = mean(daydat$LONGITUDE)
    lat = mean(daydat$LATITUDE)
    if(count > 0){
      sdlon = sd(daydat$LONGITUDE)
      sdlat = sd(daydat$LATITUDE)
      minlat = min(daydat$LATITUDE)
      maxlat = max(daydat$LATITUDE)
      if(count > 2){
        ucilat = as.numeric(CI(daydat$LATITUDE)[1], ci=0.95)
        lcilat = as.numeric(CI(daydat$LATITUDE)[2], ci=0.95)
        }
      else {
        ucilat = NA
        lcilat = NA
      }
    }
    else{
      sdlon = NA
      sdlat = NA
      minlat = NA
      maxlat = NA
      ucilat = NA
      lcilat = NA
    }
      df[outcount,1] = c(as.character(species))
    df[outcount,2:12] = c(julian[d], month, count, lat, sdlat, lon, sdlon, ucilat, lcilat, minlat, maxlat)
    outcount = outcount + 1
  }
  return(df)
}


DailyTravel = function(meanlocs){
  #use geodist to calculate Great Circle distance between daily location centers
  require(spaa)
  require(mgcv)
  dst=c(NA)
  
  for(i in 1:nrow(meanlocs)){
    if (i < nrow(meanlocs)){
      dist = geodist(meanlocs[i,7], meanlocs[i,5], meanlocs[i+1,7], meanlocs[i+1,5])
      if (is.nan(dist) == TRUE) {
        dist = 0  #appaarently geodist doesn't calculate zero change in location, records as NaN
      }
      dst = c(dst, dist)
    }
  }
  distdat = cbind(meanlocs,dst)
  
  print (ggplot(distdat, aes(jday, dst)) + geom_line(size=1) + theme_bw() + xlab("Julian Day") + 
    ylab("Distance Traveled (km)") + ggtitle(distdat[1,1]))
  
  #TODO: use La Sorte method to fit gamm and determine threshold for beginning 
  # of spring and end of fall migration for each species and year; Fig A4 in 2013 article
  print (ggplot(distdat, aes(jday, count)) + geom_point() + theme_bw() + xlab("Julian Day") + 
           ylab("Number of observances") + ggtitle(distdat[1,1]) + 
           geom_smooth(se=T, method='gam', formula=y~s(x), color='indianred'))
  
    return (dst)
}
  

DailyCentroid = function(daydat, hexgrid){
  #idenifies the centroid of location for each daily record for the species
  # place each lon-lat (row) into one of the hexes
  # replace the specific lat-long with the central coordinates for that hex
  # calculate a weighted mean location for the day
  # record: number of rows, number of hexes, weighted mean lon and lat
  numobs = nrow(daydat)
  
  meanlon = mean(daydat$LONGITUDE)
  meanlat = mean(daydat$LATITUDE)
  
  #To find the POLYFID for the hexes for each observation, pull out mean lon and lat
  coords = data.frame("mean_lon" = meanlon, "mean_lat" = meanlat)
  
  # Matches observations with the polygon hexes in the map
  ID <- over(SpatialPoints(coords), hexdat)
  coords <- cbind(coords, ID) 
  
  centroiddat = cbind(numobs, coords)
  return (centroiddat)
}


YearlyCentroid = function(dat, hexgrid) {
  #creates a new dataframe with centroid lon-lat for each julian day in a year
  year = dat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  jdate = seq(1:numdays)
  centroids = data.frame(julian=NA, meanlon=NA, meanlat=NA, JOIN_COUNT=NA, AREA=NA, PERIMETER=NA, BOB_=NA, BOB_ID=NA, ID=NA, POLYFID=NA, LONGITUDE=NA, LATITUDE=NA, numobs=NA)
  outcount = 1
  
  for (j in seq_along(jdate)){
    tmp = dat[which(dat$julian == j),]
    #if the daily data exists record weighted mean lon-lat
    if (nrow(tmp) > 0) {
      locs = DailyCentroid(tmp, hexgrid)
      head(locs)
      result = cbind(j, locs)
    }
    #if no observations for the date, record NA
    else {
     result = cbind(j, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    }
    centroids[j,] = result
    outcount = outcount + 1
  }
  return (centroids)
}

PlotChecklistMap = function(humdat, hexdat, dirpath){
  require(fields)
  #plots the hexmap with the number of checklist over the entire time period color coded. 
  #saves a jpg file in the species directory
  
  #To find the POLYFID for the hexes for each observation, pull out lon and lat
  coords = humdat[,c(10,9)]
  names(coords) = c("real_lon", "real_lat") #rename so we can tell real from hex coordinates later
  
  # Matches observations with the polygon hexes in the map
  ID <- over(SpatialPoints(coords), hexdat)
  coords <- cbind(coords, ID) 
  
  #find the number of obs in each cell
  t= table(as.factor(coords$POLYFID))
  
  #Merge the count data for the hexes with all the hexes in the map
  df = data.frame(POLYFID = names(t), count=as.numeric(t))
  df2 = data.frame(POLYFID = unique(hexdat$POLYFID))
  df3 = merge(df2, df, all.x=TRUE)
  
  #matches colors with the number of observations
  cols = data.frame(id=c(NA,sort(unique(df3$count))), cols=tim.colors(length(unique(df3$count))), stringsAsFactors=FALSE)
  df4 = merge(df3, cols, by.x="count", by.y="id")
  df5 = merge(hexdat, df4, by.x="POLYFID", by.y="POLYFID", all.x=TRUE)
  
  #hexes with no counts are white
  df5$cols <- ifelse(is.na(df5$count), "white", df5$cols)
  df5 <- df5[order(df5$POLYFID),]
  
  vls = sort(unique(round(cols$id/100)*100))
  vls[1] = 1
  cols2 = tim.colors(length(vls))
  
  #make a map with hexes colored by the number of times the species was observed in a given hex
  png(paste(dirpath,"/", "ChecklistMap.png", sep=""))
  plot(hexdat, col=df5$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1)
  axis(side=1)
  axis(side=2, las=1)
  box()
  mtext("Longitude", side=1, cex=1.4, line=2.5)
  mtext("Latitude", side=2, cex=1.4, line=2.5)
  ## legend
  legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1, cex=0.75, bty="n",
         col="black", title="Number of checklists", x.intersp=1, y.intersp=0.5)
  dev.off()  
  
  return(df5)
}
