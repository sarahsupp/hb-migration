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
  julian = sort(unique(yeardata$julian))
  
  df = data.frame("species"=NA, "jday"=1, "month"=1, "count"=1, "meanlat"=1, "sdlat"=1, "meanlon"=1, "sdlon"=1)
  outcount = 1
  
  for (d in 1:length(julian)) {
    daydat = yeardata[which(yeardata$julian == julian[d]),]
    count = nrow(daydat)
    month = daydat$month[1]
    lat = mean(daydat$LATITUDE)
    sdlat = sd(daydat$LATITUDE)
    lon = mean(daydat$LONGITUDE)
    sdlon = sd(daydat$LONGITUDE)
    df[outcount,1] = c(as.character(species))
    df[outcount,2:8] = c(julian[d], month, count, lat, sdlat, lon, sdlon)
    outcount = outcount + 1
  }
  return(df)
}

DailyTravel = function(meanlocs){
  #use geodist to calculate Great Circle distance between daily location centers
  require(spaa)
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
  
  #TODO: get a better fit line here
  #TODO: use La Sorte method to fit gamm and determine threshold for beginning 
  # of spring and end of fall migration for each species and year; Fig A4 in 2013 article
  print (ggplot(distdat, aes(jday, count)) + geom_point() + theme_bw() + xlab("Julian Day") + 
           ylab("Number of observances") + ggtitle(distdat[1,1]) + stat_smooth(col="indianred"))
  
    return (dst)
}
  

DailyCentroid = function(daydat, hexgrid){
  #idenifies the centroid of location for each daily record for the species
  # place each lon-lat (row) into one of the hexes
  # replace the specific lat-long with the central coordinates for that hex
  # calculate a weighted mean location for the day
  # record: number of rows, number of hexes, weighted mean lon and lat
  numobs = nrow(datday)
  numcells = length(unique(hexID))
  centroiddat = c(numobs, numcells, centr_lon, centr_lat)
  return (centroiddat)
}

YearlyCentroid = function(yrdat, hexgrid) {
  #creates a new dataframe with the daily centroid location
  
  year = yrdat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  jdate = seq(1:numdays)
  cntr_lon = NULL
  cntr_lat = NULL
  
  for (j in seq_along(jdate)){
    tmp = yrdat[which(yrdat$julian == j),]
    #if the daily data exists record weighted mean lon-lat
    if (nrow(tmp) > 0) {
      locs = DailyCentroid(tmp, hexgrid)
      cntr_lon = c(cntr_lon, )
      cntr_lat = c(cntr_lat, )
    }
    #if no observations for the date, record NA
    else {
      cntr_lon = c(cntr_lon, NA)
      cntr_lat = c(cntr_lat, NA)
    }
  }
  centroids = cbind(jdate, cntr_lon, cntr_lat)
  return (centroids)
}
