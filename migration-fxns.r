# Define functions to use in hb-migration code
require(chron)
require(mgcv)
require(Rmisc)
require(spaa)
require(ggplot2)


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
 
  barplot = ggplot(data = t, aes(year, count)) + geom_bar(fill="cadetblue2") + theme_bw() + 
    ggtitle(species) + theme(text = element_text(size=20)) + xlab("time") + ylab("number checklists")
  
  print(barplot)
  
  return(t)
}


AlternateMeanLocs = function(dat, species, hexdat) {
  #input is a yearly data frame for a species
  #output is a new dataframe with mean lat and long, calculated from hexagon centers, after placing daily observed
  #lat and long onto the isacohedron equal-area cells map, sensu La Sorte et al. 2013
  df = data.frame("jday"=1, "month"=1, "numobs"=1, "numcells"=1, "centerlon"=1, "centerlat"=1)
  outcount = 1
  
  #count the number of days in the year
  year = dat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  julian = seq(1:numdays)
  
  #To find the POLYFID for the hexes for each observation, identify observed lon and lat  
  # Matches observations with the polygon hexes in the map
  ID <- over(SpatialPoints(dat[,c(10,9)]), hexdat)
    names(ID) = c("JOIN_COUNT", "AREA", "PERIMETER", "BOB_", "BOB_ID", "ID", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE")
  coords <- cbind(dat, ID) 
  
  #aggregate daily info by mean centroid location
  dailyHexes = count(coords, vars=c("julian", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE"))
  
  #calculate weighted mean (weights based on number of obs per hex) lat and long
  for (j in 1:length(julian)){
    jdata = dailyHexes[which(dailyHexes$julian == j),]
    if (nrow(jdata) > 0){
      numcells = nrow(jdata)
      numobs = sum(jdata$freq)
      mo = as.numeric(months(j))
      wtmean_lon = weighted.mean(jdata$HEX_LONGITUDE, jdata$freq, na.rm=TRUE)
      wtmean_lat = weighted.mean(jdata$HEX_LATITUDE, jdata$freq, na.rm=TRUE)
      df[outcount,] = c(j, mo, numobs, numcells, wtmean_lon, wtmean_lat)
      outcount = outcount + 1
    }
    else {
      mo = as.numeric(months(j))
      df[outcount,] = c(j, mo, 0, 0, NA, NA)
      outcount = outcount + 1
    }
  }
  return(df)
}


MeanDailyLoc = function(dat, species) {
  # input is yearly data frame for a species
  # makes a new dataframe with mean julian day lat and long, sd, and 95% confidence intervals
  require(Rmisc)
  require(chron)
  
  #count the number of days in the year
  year = dat$year[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  julian = seq(1:numdays)
  
  df = data.frame("species"=NA, "jday"=1, "month"=1, "count"=1, "meanlat"=1, "sdlat"=1, "meanlon"=1, 
                  "sdlon"=1, "ucilon"=1, "lcilon"=1, "minlat"=1, "maxlat"=1)
  outcount = 1
  
  for (d in 1:length(julian)) {
    daydat = dat[which(dat$julian == julian[d]),]
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
        ucilon = as.numeric(CI(daydat$LONGITUDE)[1], ci=0.95)
        lcilon = as.numeric(CI(daydat$LONGITUDE)[2], ci=0.95)
        }
      else {
        ucilon = NA
        lcilon = NA
      }
    }
    else{
      sdlon = NA
      sdlat = NA
      minlat = NA
      maxlat = NA
      ucilon = NA
      lcilon = NA
    }
      df[outcount,1] = c(as.character(species))
    df[outcount,2:12] = c(julian[d], month, count, lat, sdlat, lon, sdlon, ucilon, lcilon, minlat, maxlat)
    outcount = outcount + 1
  }
  return(df)
}


DailyTravel = function(meanlocs, loncol, latcol, species, year){
  #use geodist to calculate Great Circle distance between daily location centers
  require(spaa)
  require(mgcv)   #TODO: GEt some different values when use spaa vs gmt package to calculate great circle distance. Track down issue
  #require(gmt)
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
  
  print (ggplot(distdat, aes(jday, dst)) + geom_line(size=1, col = "#4daf4a") + theme_bw() + xlab("Julian Day") + 
           ylab("Distance Traveled (km)") + ggtitle(paste(species, year, sep = " ")))
  
  return(distdat)
}


GetMigrationDates = function(data) {
  # uses a generalized additive model (GAM) to define the start of spring and end of fall migration
  
  #GAM model on number of cells with observations by julian date
  gam1 = gam(numcells ~ s(jday, k = 40), data = data, gamma = 1.5) 
  xpred = data.frame(jday = c(1:max(data$jday)))
  dpred = predict(gam1, newdata=xpred, type="response", se.fit=T)
  
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
  spring = spring_index2+1
  
  #identify end of fall migration
  tst <- 1000
  fall_index2 = fall_max
  while(tst > fall_threshold){
    tst = dpred$fit[fall_index2]
    if(fall_index2==365) break
    fall_index2 <- fall_index2 + 1
  }
  fall <- fall_index2-1
  
  dates = list(spring=spring, fall=fall)
  return(dates)
}


PlotOccurrences = function(data, species, spring, fall) {
  #takes a dataframe with julian dates, and the number of cells, and returns a plot with lines fitted
  # for a gam and for the onset of spring and end of fall migration
  # sensu La Sorte and Fink code from 2013 paper
  
  occ = ggplot(data, aes(jday, numcells)) + geom_point() + xlab("julian day") + ylab("number of cells") + 
    geom_smooth(se=T, method='gam', formula=y~s(x, k=40), gamma=1.5, col='#7b3294', fill='#af8dc3') + 
    theme_bw() + ggtitle(species) + scale_x_continuous(breaks = seq(0, 365, by = 25)) +
    scale_y_continuous(breaks = seq(0, 80, by = 5)) + 
    geom_vline(xintercept = c(spring, fall), col = "#008837", linetype = "dashed", size = 1) +
    theme(text = element_text(size=20))
  
  print(occ)
  return(occ)
}


EstimateDailyLocs = function(dat) {
  #input a dataframe with the mean daily locations for the year, and uses a GAM smoothing function
  # to estimate daily occurrence lat and long separately. Binds the fitted lat and long values together,
  # with standard errors. Returns a dataframe.
      #TODO: check choice of k and gamma in the GAM function
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
  

FindCentroids = function(dat, loncol, latcol, hexgrid){
  # matches centroids with the mean lat and lons
  
  #To find the POLYFID for the hexes for each observation, pull out mean lon and lat
  coords = dat[,c(loncol,latcol)]
    names(coords) = c("MeanLon", "MeanLat")
    coords = coords[complete.cases(coords),] #next step can't deal with NA values
  jday = dat[complete.cases(dat$meanlon),]$jday
  
  # Matches observations with the polygon hexes in the map
  ID <- over(SpatialPoints(coords), hexdat)
  coords <- cbind(jday, coords, ID) 
  
  df = merge(dat, coords, by = "jday", all.x = TRUE)
  df = df[,-which(names(df) %in% c("MeanLon", "MeanLat"))] #delete repetitive columns
  
  return(df)
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
