# Define functions to use in hb-migration code
require(chron)
require(mgcv)
require(Rmisc)
require(spaa)
require(ggplot2)
require(SDMTools)


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
    
  julian = sapply(date, function(x){
    julian(as.numeric(substring(x, 6, 7)), as.numeric(substring(x, 9, 10)), as.numeric(substring(x, 1, 4)), 
    origin. = c(1, 1, as.numeric(substring(x, 1, 4)))) + 1
  })
  
  newdate = as.data.frame(cbind(year, month, day, julian), row.names=FALSE)

  return (newdate)
}


PlotRecords = function(yeardata, species) {
  #plots the number of records across the years
  
  t = as.data.frame(table(yeardata))
  names(t) = c('year', 'count')
 
  barplot = ggplot(data = t, aes(year, count)) + geom_bar(fill="cadetblue2", stat="identity") + theme_bw() + 
    ggtitle(species) + theme(text = element_text(size=20)) + xlab("time") + ylab("number checklists")
  
  print(barplot)
  
  return(t)
}


AlternateMeanLocs = function(dat, species, hexdat, yreffort) {
  #input is a yearly data frame for a species
  #output is a new dataframe with mean lat and long, calculated from hexagon centers, after placing daily observed
  #lat and long onto the isacohedron equal-area cells map, sensu La Sorte et al. 2013
  df = data.frame("jday"=1, "month"=1, "numobs"=1, "numcells"=1, "centerlon"=1, "centerlat"=1, "sdlon"=1, "sdlat"=1)
  outcount = 1
  
  #count the number of days in the year
  year = dat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  julian = seq(1:numdays)
  
  #To find the POLYFID for the hexes for each observation, identify observed lon and lat  
  # Matches observations with the polygon hexes in the map
  ID <- over(SpatialPoints(dat[,c(10,9)]), hexdat)  #TODO: Check with TAC that this is working correctly based on map projection
    names(ID) = c("JOIN_COUNT", "AREA", "PERIMETER", "BOB_", "BOB_ID", "ID", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE")
  coords <- cbind(dat, ID) 
  
  #aggregate daily info by mean centroid location
  dailyHexes = count(coords, vars=c("julian", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE"))
  
  #calculate weighted mean (weights based on number of obs per hex) lat and long
  for (j in 1:length(julian)){
    jdata = dailyHexes[which(dailyHexes$julian == j),]
    jeffort = yreffort[which(yreffort$DAY == j),]
    jdata = merge(jdata, jeffort, by.x = "POLYFID", by.y = "POLYFID")  #TODO - DOESN'T ALWAYS MATCH UP - WHY ARE SOME MISSING? (e.g. 2012 Ar. alex., day 100)
    if (nrow(jdata) > 0){                                               # POLYFID c(8378, 8880, 9386, 9051) is in jdata, but not in jeffort - why?
      numcells = nrow(jdata)
      numobs = sum(jdata$freq)
      mo = as.numeric(months(j))
      wtmean_lon = wt.mean(jdata$HEX_LONGITUDE, jdata$freq/jdata$COUNT)
      wtmean_lat = wt.mean(jdata$HEX_LATITUDE, jdata$freq/jdata$COUNT)
      lon_sd <- wt.sd(jdata$HEX_LONGITUDE, jdata$freq/jdata$COUNT)
      lat_sd <- wt.sd(jdata$HEX_LATITUDE, jdata$freq/jdata$COUNT)
      #wtmean_lon = weighted.mean(jdata$HEX_LONGITUDE, jdata$freq, na.rm=TRUE)
      #wtmean_lat = weighted.mean(jdata$HEX_LATITUDE, jdata$freq, na.rm=TRUE)
      df[outcount,] = c(j, mo, numobs, numcells, wtmean_lon, wtmean_lat,lon_sd, lat_sd)
      outcount = outcount + 1
    }
    else {
      mo = as.numeric(months(j))
      df[outcount,] = c(j, mo, 0, 0, NA, NA, NA, NA)
      outcount = outcount + 1
    }
  }
  return(df)
}


DailyTravel = function(meanlocs, loncol, latcol, species, year, migr_dates){
  #use geodist to calculate Great Circle distance between daily location centers
  require(spaa)
  require(mgcv)
  dst=c(NA)
  
  for(i in 1:nrow(meanlocs)){
    if (i < nrow(meanlocs)){
      dist = geodist(meanlocs[i,loncol], meanlocs[i,latcol], meanlocs[i+1,loncol], meanlocs[i+1,latcol])
      if (is.nan(dist) == TRUE) {
        dist = 0  #appaarently geodist doesn't calculate zero change in location, records as NaN
      }
      dst = c(dst, dist)
    }
  }
  distdat = cbind(meanlocs,dst)
  
  subdistdat = distdat[which(distdat$jday >= migr_dates[1] & distdat$jday <= migr_dates[2]),]

  print(ggplot(subdistdat, aes(jday, dst)) + geom_line(size=1, col = "#4daf4a") + theme_bw() + xlab("Julian Day") + 
          ylab("Distance Traveled (km)") + ggtitle(paste(species, year, sep = " ")) +
          theme(text = element_text(size=20)) +    
          #geom_vline(xintercept = c(migr_dates[1], migr_dates[2]), linetype = "dashed", size = 1) +
          geom_vline(xintercept = median(migr_dates), col = "indianred", linetype = "dashed"))

  return(distdat)
}

GetBreedingDates = function(dat, migration_dates){
  # uses a generalized additive model (GAM) on latitude to estimate the end of spring and 
  # begin of fall migration, which defines the time at the breeding grounds
  
  #get median migration date
  med = median(migration_dates)
  
  #GAM model on predicted latitude of centroids by julian date
  gam1 = gam(lat ~ s(jday, k = 40), data = dat, gamma = 1.5) 
  xpred = data.frame(jday = c(1:max(dat$jday)))
  dpred = predict(gam1, newdata=xpred, type="response", se.fit=TRUE)
  
  lat_threshold = min(dpred$se.fit[c(med-30:med+30)]*2.56 + dpred$fit[c(med-30:med+30)])
  spring_index = med-30:med
  fall_index = med:med+30
  spring_max = spring_index[which.max(dpred$fit[spring_index])]
  fall_max = fall_index[which.max(dpred$fit[fall_index])]
  
  #identify beginning of spring migration
  tst = 1000
  spring_index2 = spring_max
  while(tst > lat_threshold){
    tst = dpred$fit[spring_index2]
    if(spring_index2 == 1) break
    spring_index2 = spring_index2 - 1
  }
  spring = spring_index2 + 1
  
  #identify end of fall migration
  tst <- 1000
  fall_index2 = fall_max
  while(tst > lat_threshold){
    tst = dpred$fit[fall_index2]
    if(fall_index2==365) break
    fall_index2 <- fall_index2 + 1
  }
  fall <- fall_index2 - 1
  
  dates = c(spring, fall)
  return(dates)
}

GetMigrationDates = function(data) {
  # uses a generalized additive model (GAM) to define the start of spring and end of fall migration
  
  #GAM model on number of cells with observations by julian date
  gam1 = gam(numcells ~ s(jday, k = 40), data = data, gamma = 1.5) 
  xpred = data.frame(jday = c(1:max(data$jday)))
  dpred = predict(gam1, newdata=xpred, type="response", se.fit=TRUE)
  
  ## cutoff based on 2 SE for spring and fall combined, following La Sorte et al. 2013 methods
  # Spring migration should be between 11 Jan and 9 July
  # Fall migration should be between 8 August and 21 Dec
  spring_threshold = min(dpred$se.fit[c(1:120)]*2.56 + dpred$fit[c(1:120)])
  fall_threshold = min(dpred$se.fit[c(280:365)]*2.56 + dpred$fit[c(280:365)])
  spring_index = 11:190
  fall_index = 220:355
  spring_max = spring_index[which.max(dpred$fit[spring_index])]
  fall_max = fall_index[which.max(dpred$fit[fall_index])]
  
  #identify beginning of spring migration
  tst = 1000
  spring_index2 = spring_max
  while(tst > spring_threshold){
    tst = dpred$fit[spring_index2]
    if(spring_index2 == 1) break
    spring_index2 = spring_index2 - 1
  }
  spring = spring_index2 + 1
  
  #identify end of fall migration
  tst <- 1000
  fall_index2 = fall_max
  while(tst > fall_threshold){
    tst = dpred$fit[fall_index2]
    if(fall_index2==365) break
    fall_index2 <- fall_index2 + 1
  }
  fall <- fall_index2 - 1
  
  dates = c(spring, fall)
  return(dates)
}


PlotOccurrences = function(data, species, spring, fall) {
  #takes a dataframe with julian dates, and the number of cells, and returns a plot with lines fitted
  # for a gam and for the onset of spring and end of fall migration
  # sensu La Sorte and Fink code from 2013 paper
  
  med = median(c(spring,fall))
  occ = ggplot(data, aes(jday, numcells)) + geom_point() + xlab("julian day") + ylab("number of cells") + 
    geom_smooth(se=T, method='gam', formula=y~s(x, k=40), gamma=1.5) + 
    theme_bw() + ggtitle(species) + scale_x_continuous(breaks = seq(0, 365, by = 25)) +
    scale_y_continuous(breaks = seq(0, 80, by = 5)) + 
    geom_vline(xintercept = spring, col = "cadetblue", size = 1) +
    geom_vline(xintercept = fall, col = "orange", size = 1) +
    geom_vline(xintercept = med, col = "indianred", linetype = "dashed") +
    theme(text = element_text(size=20))
  
  print(occ)
  return(occ)
}


EstimateDailyLocs = function(dat) {
  #input a dataframe with the mean daily locations for the year, and uses a GAM smoothing function
  # to estimate daily occurrence lat and long separately. Binds the fitted lat and long values together,
  # with standard errors. Returns a dataframe.
      #check choice of k and gamma in the GAM function - same as FAL 2013
  #find the best fit line for the data, for longitude and latitude separately
  lon_gam <- gam(centerlon ~ s(jday, k=10), data = dat, gamma = 1.5)
  lat_gam <- gam(centerlat ~ s(jday, k=10), data = dat, gamma = 1.5)
  
  #predict values along the smoothing line
  xpred <- data.frame(jday=sort(unique(dat$jday)))
  lonpred <- predict(lon_gam, newdata = xpred, type="response", se.fit=T)
  latpred <- predict(lat_gam, newdata = xpred, type="response", se.fit=T)
  
  #bring the data back together
  preds =  data.frame(spname = species, jday = xpred$jday, month = dat$month, lon = lonpred$fit, 
                      lat = latpred$fit, lon_se = lonpred$se.fit, lat_se = latpred$se.fit)
  
  return(preds)
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
  df5$cols = ifelse(is.na(df5$count), "white", df5$cols)
  df5 = df5[order(df5$POLYFID),]
  
  vls = sort(unique(round(cols$id/500)*500))
  vls[1] = 1
  cols2 = tim.colors(length(vls))
  
  #make a map with hexes colored by the number of times the species was observed in a given hex
  pdf(paste(dirpath,"/", "ChecklistMap.pdf", sep=""))
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


PlotMigrationPath = function(dat, map, species, year) {
  # makes a nice figure showing the migration trajectory of the GAM predicted daily locations, with se for longitude
  dat$month = as.factor(dat$month)
  
  map = ggmap(map) + geom_line(data = dat, aes(lon, lat, col=month)) + 
    geom_errorbarh(data = dat, aes(xmin = lon - lon_se, xmax = lon + lon_se, col = month)) + 
    xlab("Longitude") + ylab("Latitude") + ggtitle(paste(species, year, sep = " "))
  
  print(map)
  return(map)
}


PlotMeanLatitude = function(dat, species, year){
  #plot mean latitude for each julian day, point size represents number of checklists, returns the map object
  meanlat = ggplot(dat, aes(jday, meanlat, col=as.factor(month))) + geom_point(aes(size=count)) + 
    ggtitle(paste(species, year, "latitudinal migration", sep = " ")) + xlab("Julian Day") + ylab("Mean Latitude") +
    theme_bw() +  scale_x_continuous(breaks = seq(0, 365, by = 25)) +
    scale_y_continuous(breaks = seq(10, 80, by = 5)) +
    geom_smooth(se=T, method='gam', formula=y~s(x), color='indianred')
  return(meanlat)
}


PlotAllPoints = function (dat, map, species, year){
  #plots all observed locations for hummingbird sightings
  dat$month = as.factor(dat$month)
  sitemap = ggmap(map) + geom_point(data=dat, aes(LONGITUDE, LATITUDE, col=month)) + 
            ggtitle(paste(species, year, sep = " "))
  print(sitemap)
  return(sitemap)
}

MigrationSpeed = function(dat, migration){
  # estimates daily migration speed for spring and fall, separately
  # takes the top 5 fastest migration speeds for each time period, and assigns the median as the migration speed
  # migration has two elements, beginning of spring migration and end of fall migration
  med = median(migration)
  springdat = dat[which(dat$jday >= migration[1] & dat$jday < med),]
  falldat = dat[which(dat$jday <= migration[2] & dat$jday > med),]
  
  springspeed = median(sort(springdat$dst, decreasing=TRUE)[1:5])
  fallspeed = median(sort(falldat$dst, decreasing=TRUE)[1:5])

  return(c(springspeed,fallspeed))
}


FindMismatch = function(dat, species, hexdat, yreffort, map) {
  #input is a yearly data frame for a species
  #output is a new dataframe with mean lat and long, calculated from hexagon centers, after placing daily observed
  #lat and long onto the isacohedron equal-area cells map, sensu La Sorte et al. 2013

  #count the number of days in the year
  year = dat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  julian = seq(1:numdays)
  
  #To find the POLYFID for the hexes for each observation, identify observed lon and lat  
  # Matches observations with the polygon hexes in the map
  ID <- over(SpatialPoints(dat[,c(10,9)]), hexdat)  #TODO: Check with TAC that this is working correctly based on map projection
  names(ID) = c("JOIN_COUNT", "AREA", "PERIMETER", "BOB_", "BOB_ID", "ID", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE")
  coords <- cbind(dat, ID) 
  
  #aggregate daily info by mean centroid location
  dailyHexes = count(coords, vars=c("julian", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE"))
  
  #calculate weighted mean (weights based on number of obs per hex) lat and long
  missingdat = data.frame("julian"=1, "POLYFID"=1, "H"=1, "numcells"=1, "centerlon"=1, "centerlat"=1, "sdlon"=1, "sdlat"=1)
  outcount = 1
  
  for (j in 1:length(julian)){
    jdata = dailyHexes[which(dailyHexes$julian == j),]
    jeffort = yreffort[which(yreffort$DAY == j),]
    hexes = unique(jeffort$POLYFID)
    misses = jdata[which(!jdata$POLYFID %in% hexes),]
    if(j == 1){
      missing = misses
    }
    else {
      missing = rbind(missing, misses)
    }
    }
  missing = missing[complete.cases(missing),]
  
  errmap = ggmap(map) + geom_point(data=missing, aes(HEX_LONGITUDE, HEX_LATITUDE)) + 
    ggtitle(paste(species, year, sep = " "))
  
  print(errmap)
}
