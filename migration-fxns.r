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
  
  newdate = as.data.frame(cbind(year,month,day))
  names(newdate) = c("year", "month", "day")

  newdate = JulianDay(newdate)

  return (newdate)
}

JulianDay = function(daymoyr){
  ## Add a column for julian day
  require(chron)
  for (row in 1:nrow(daymoyr)){
    line = daymoyr[row,]
    daymoyr$julian[row] = julian(line$month, line$day, line$year, origin. = c(month=1, day=1, year=line$year))
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
