# Define functions to use in hb-migration code
require(chron)
require(mgcv)
require(Rmisc)
require(spaa)
require(ggplot2)
require(SDMTools)


PlotRecords = function(yeardata, species) {
  #plots the number of records across the years
  
  t = as.data.frame(table(yeardata))
  names(t) = c('year', 'count')
 
  barplot = ggplot(data = t, aes(year, count)) + geom_bar(fill="cadetblue2", stat="identity") + theme_bw() + 
    ggtitle(species) + theme(text = element_text(size=20)) + xlab("time") + ylab("number checklists") + 
    theme(axis.text.x=element_text(angle = -60, hjust = 0))
  
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
  year = dat$YEAR[1]
  numdays = as.numeric(as.POSIXlt(paste(year, "-12-31", sep = "")) - as.POSIXlt(paste(year, "-01-01", sep="")) + 1)
  julian = seq(1:numdays)

  #merge data on POLYFID so it includes the center of each hex cell (LONGITUDE.y, LATITUDE.y)
  dat_merged = merge(dat, hexdat, by.x = "POLYFID", by.y = "POLYFID")
  
  #aggregate daily info by mean centroid location
  dailyHexes = count(dat_merged, vars=c("DAY", "POLYFID", "LONGITUDE.y", "LATITUDE.y"))
  
  #calculate weighted mean (weights based on number of obs per hex and total effort per hex) LAT and LON
  for (j in 1:length(julian)){
    jdata = dailyHexes[which(dailyHexes$DAY == j),]
    jeffort = yreffort[which(yreffort$DAY == j),]
    jdata = merge(jdata, jeffort, by.x = c("POLYFID", "DAY"), by.y = c("POLYFID", "DAY"))
    if (nrow(jdata) > 0){                                              
      numcells = nrow(jdata)
      numobs = sum(jdata$freq)
      mo = as.numeric(months(j))
      wtmean_lon = wt.mean(jdata$LONGITUDE.y, jdata$freq/jdata$COUNT)
      wtmean_lat = wt.mean(jdata$LATITUDE.y, jdata$freq/jdata$COUNT)
      lon_sd <- wt.sd(jdata$LONGITUDE.y, jdata$freq/jdata$COUNT)
      lat_sd <- wt.sd(jdata$LATITUDE.y, jdata$freq/jdata$COUNT)
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
  
  subdistdat = distdat[which(distdat$jday >= migr_dates[[1]] & distdat$jday <= migr_dates[[3]]),]

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
  
  # cutoff based on 2 SE for the breeding period, defined as +/- 30 days from the median date,
  # based on start of spring and end of fall migration in GetMigrationDates function.
  lat_threshold = min(dpred$se.fit[c((med-30):(med+30))]*2.56 + dpred$fit[c((med-30):(med+30))])
  spring_index = (med-30):med
  fall_index = med:(med+30)
  spring_max = spring_index[which.max(dpred$fit[spring_index])]
  fall_max = fall_index[which.max(dpred$fit[fall_index])]
  
  #identify end of spring migration
  tst = 1000
  spring_index2 = spring_max
  while(tst > lat_threshold){
    tst = dpred$fit[spring_index2]
    if(spring_index2 == 1) break
    spring_index2 = spring_index2 - 1
  }
  spring = spring_index2 + 1
  
  #identify beginning of fall migration
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


PlotOccurrences = function(data, species, year, migration) {
  #takes a dataframe with julian dates, and the number of cells, and returns a plot with lines fitted
  # for a gam and for the onset of spring and end of fall migration
  # sensu La Sorte and Fink code from 2013 paper + new segmented approach
  
  occ = ggplot(data, aes(jday, numcells)) + geom_point() + xlab("julian day") + ylab("number of cells") + 
    geom_smooth(se=T, method='gam', formula=y~s(x, k=40), gamma=1.5) + 
    theme_bw() + ggtitle(paste(species, year, sep=" ")) + scale_x_continuous(breaks = seq(0, 365, by = 25)) +
    scale_y_continuous(breaks = seq(0, 80, by = 5)) + 
    geom_vline(xintercept = migration[1], col = "cadetblue", size = 1) +
    geom_vline(xintercept = migration[3], col = "orange", size = 1) +
    geom_vline(xintercept = migration[2], col = "olivedrab3", size = 1) +
    theme(text = element_text(size=12))+
    theme(axis.text.x=element_text(angle = -60, hjust = 0))
  
  print(occ)
  return(occ)
}


EstimateDailyLocs = function(dat) {
  #input a dataframe with the mean daily locations for the year, and uses a GAM smoothing function
  # to estimate daily occurrence lat and long separately. Binds the fitted lat and long values together,
  # with standard errors. Returns a dataframe.
      #check choice of k and gamma in the GAM function - same as FAL 2013
  #find the best fit line for the data, for longitude and latitude separately
  lon_gam = gam(centerlon ~ s(jday, k=20), data = dat, gamma = 1.5) # TODO: consider upping gamma? and changing basis? (adaptive spline or penalized spline?)
  lat_gam = gam(centerlat ~ s(jday, k=20), data = dat, gamma = 1.5)
  
  #predict values along the smoothing line
  xpred = data.frame(jday=sort(unique(dat$jday)))
  lonpred = predict(lon_gam, newdata = xpred, type="response", se.fit=T)
  latpred = predict(lat_gam, newdata = xpred, type="response", se.fit=T)
  
  #bring the data back together
  preds =  data.frame(spname = species, jday = xpred$jday, month = dat$month, lon = lonpred$fit, 
                      lat = latpred$fit, lon_se = lonpred$se.fit, lat_se = latpred$se.fit)
  
  return(preds)
}
  

PlotChecklistMap = function(humdat, hexdat, dirpath){
  require(fields)
  #plots the hexmap with the number of checklist over the entire time period color coded. 
  #saves a pdf file in the species directory

  #find the number of obs in each cell
  t = table(as.factor(humdat$POLYFID))
 
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
  
  #set scale for legend
  vls = sort(unique(round(cols$id/500)*500))
  vls[1] = 1
  cols2 = tim.colors(length(vls))
  
  #make a map with hexes colored by the number of times the species was observed in a given hex
  setEPS()
  postscript(paste(dirpath,"/", "ChecklistMap.eps", sep=""), width=10, height=7)
  plot(hexdat, col=df5$cols, border="gray50", lwd=0.25, xlim=c(-170,-50), ylim=c(15,75), las=1)
  axis(side=1)
  axis(side=2, las=1)
  box()
  mtext("Longitude", side=1, cex=1.4, line=2.5)
  mtext("Latitude", side=2, cex=1.4, line=2.5)
  ## legend
  legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1, cex=0.75, bty="",
         col="black", title="Number of checklists", x.intersp=1, y.intersp=0.5, bg="white")
  map("worldHires", c("usa", "canada", "mexico"), add=TRUE)
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
    scale_y_continuous(breaks = seq(15, 75, by = 5)) +
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
  # breeding has two elements, end of spring migration and beginning of fall migration from breeding grounds

  springdat = dat[which(dat$jday >= migration[[1]] & dat$jday <= migration[[2]]),]
  falldat = dat[which(dat$jday >= migration[[2]] & dat$jday <= migration[[3]]),]
  
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
  missing = missing[complete.cases(missing),] #TODO: Are the missing points because of some mismatch in what is being counted towards effort?
  nadat = missing[!complete.cases(missing),] #TODO: Are the NA points coming from the coast? TAC found some coastal mismatch - worth trying to put into nearest neighbor cell?
  
  errmap = ggmap(map) + geom_point(data=missing, aes(HEX_LONGITUDE, HEX_LATITUDE)) + 
    ggtitle(paste(species, year, sep = " "))
  
  print(errmap)
  return(list(missing, nadat))
}


LinearMigration = function(seasondat, year){
  # takes the data from a single season and runs a linear regression on the lat or lon
  # returns slope and r2 fit
  
  lm1 = lm(lat ~ jday, data = seasondat)
  lm2 = lm(lon ~ jday, data = seasondat)
  
  dat = data.frame("year" = year, "lat_slope" = lm1$coef[[2]], "lat_r2" = summary(lm1)$r.squared, 
                   "lon_slope" = lm2$coef[[2]], "lon_r2" = summary(lm2)$r.squared)
  return(dat)
}


GroupDuplicates = function(humdat) { 
  #gets rid of duplicate records that are part of the same group. 
  #Takes only the first record to move to analysis, keeps all single records that are not part of a group.
  gid = sort(unique(humdat$GROUP_ID))
  #gid = gid[which(gid!=NA)]
  
  #open start a new dataframe with same columns as the main dataframe
  keep = humdat[1,]
  out = 0
  
  for (g in 1:length(gid)){
    out = out + 1
    tmp = humdat[which(humdat$GROUP_ID == gid[g]),]
    #record the first line of the data (assume the group has the same information)
    if (nrow(tmp) == 1) { 
      keep[out,] = tmp
    }
    else{
      keep[out,] = tmp[1,]
    }
  }
  
  keepnongroup = humdat[which(is.na(humdat$GROUP_ID)),]
  keep = rbind(keep, keepnongroup)
  return(keep)
}


Est3MigrationDates = function(dat){
  #takes in predicted centroids for migration path, and estimates the beginning of spring migration,
  # the end of fall migration, and the date where the species reaches maximum latitude.
  
  #GAM model on predicted latitude of centroids by julian date
  gam1 = gam(centerlat ~ s(jday, k = 40), data = dat, gamma = 1.5) 
  xpred = data.frame(jday = c(1:max(dat$jday)))
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
  spring_begin = spring_index2 + 1
  
  #identify end of fall migration
  tst <- 1000
  fall_index2 = fall_max
  while(tst > fall_threshold){
    tst = dpred$fit[fall_index2]
    if(fall_index2==365) break
    fall_index2 <- fall_index2 + 1
  }
  fall_end <- fall_index2 - 1
  
  max_lat =   xpred$jday[which.max(dpred$fit[xpred$jday])]
  
  dates = c(spring_begin, max_lat, fall_end)
  return(dates)
}


Est4MigrationDates = function(dat){
  #takes in predicted centroids for migration path, and estimates the beginnig and end of spring migration,
  # and the beginning and end of fall migration
  
  #GAM model on predicted latitude of centroids by julian date
  gam1 = gam(centerlat ~ s(jday, k = 40), data = dat, gamma = 1.5) 
  xpred = data.frame(jday = c(1:max(dat$jday)))
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
  spring_begin = spring_index2 + 1
  
  #identify end of fall migration
  tst <- 1000
  fall_index2 = fall_max
  while(tst > fall_threshold){
    tst = dpred$fit[fall_index2]
    if(fall_index2==365) break
    fall_index2 <- fall_index2 + 1
  }
  fall_end <- fall_index2 - 1
  
  med = median(c(spring_begin, fall_end))
  # cutoff based on 2 SE for the breeding period, defined as +/- 30 days from the median date,
  # based on start of spring and end of fall migration in GetMigrationDates function
  lat_threshold = min(dpred$se.fit[c((med-30):(med+30))]*2.56 + dpred$fit[c((med-30):(med+30))])
  spring_index = (med-30):med
  fall_index = med:(med+30)
  spring_max = spring_index[which.max(dpred$fit[spring_index])]
  fall_max = fall_index[which.max(dpred$fit[fall_index])]
  
  #identify end of spring migration
  tst = 1000
  spring_index2 = spring_max
  while(tst > lat_threshold){
    tst = dpred$fit[spring_index2]
    if(spring_index2 == 1) break
    spring_index2 = spring_index2 - 1
  }
  spring_end = spring_index2 + 1
  
  #identify beginning of fall migration
  tst <- 1000
  fall_index2 = fall_max
  while(tst > lat_threshold){
    tst = dpred$fit[fall_index2]
    if(fall_index2==365) break
    fall_index2 <- fall_index2 + 1
  }
  fall_begin <- fall_index2 - 1
  
  dates = c(spring_begin, spring_end, fall_begin, fall_end)
  return(dates)
}


#Base plot migration route, trimmed by migration dates
BasePlotMigration = function(preds, yrdat, migration, elev, USAborder, Mexborder, Canborder, myext){
  
  predpts = preds[which(preds$jday > migration[1] & preds$jday < migration[3]),]
  
  predpts$month = as.factor(predpts$month)
  cols3 = data.frame(id=c(sort(unique(predpts$month))), cols=tim.colors(length(unique(predpts$month))), stringsAsFactors=FALSE)
  predpts = merge(predpts, cols3, by.x="month", by.y="id")
  #set color scale
  vls = sort(unique(round(cols3$id)))
  vls[1] = 1
  cols4 = tim.colors(length(vls))
  
  plot(predpts$lon, predpts$lat, xlim = c(-175, -50), ylim = c(15, 75), 
       xlab = "longitude", ylab = "latitude", main = "summarized migration route", cex.lab = 1.5, cex.axis = 1.5)
  plot(USAborder, ext=myext, border="gray30", add=TRUE)
  plot(Mexborder, ext=myext, border="gray30", add=TRUE)
  plot(Canborder, ext=myext, border="gray30", add=TRUE)
  points(yrdat$LONGITUDE, yrdat$LATITUDE, pch=16, col = "grey60", cex = 0.5)
  points(predpts$lon, predpts$lat, pch=19, col = predpts$cols, cex = 0.5)
  #legend("bottomleft", legend=vls, pch=22, pt.bg=cols4, pt.cex=1.25, cex=1.25, bty="n",
  #       col="black", title="", x.intersp=0.5, y.intersp=0.25)
}


ElevPlotMigration = function(preds, yrdat, migration, elev, USAborder, Mexborder, Canborder, myext) {
  
  predpts = preds
  
  predpts$month = as.factor(predpts$month)
  cols3 = data.frame(id=c(sort(unique(predpts$month))), cols=tim.colors(length(unique(predpts$month))), stringsAsFactors=FALSE)
  predpts = merge(predpts, cols3, by.x="month", by.y="id")
  #set color scale
  vls = sort(unique(round(cols3$id)))
  vls[1] = 1
  cols4 = tim.colors(length(vls))
  
  fun <- function(){
    plot(USAborder, ext=myext, border="gray30", add=TRUE)
    plot(Mexborder, ext=myext, border="gray30", add=TRUE)
    plot(Canborder, ext=myext, border="gray30", add=TRUE)
    #points(yrdat$LONGITUDE, yrdat$LATITUDE, pch=16, col = "black", cex = 0.25)
    points(predpts$lon, predpts$lat, col=predpts$cols, pch=19, cex=0.5)
  }
  
  plot.new()
  
  plot(elev, ext=myext, addfun=fun, ylab="Latitude", xlab="Longitude", main="summarized migration route", 
       cex.lab = 1.5, cex.axis=1.5, col=gray(10:256/256))
  #legend("bottomleft", legend=vls, pch=22, pt.bg=cols4, pt.cex=1.25, cex=1.25, bty="n",
   #      col="black", title="", x.intersp=0.5, y.intersp=0.25)
}


AllMigration = function(preds, elev, myext, species) {
  
  predpts = preds
  
  predpts$month = as.factor(predpts$month)
  cols3 = data.frame(id=c(sort(unique(predpts$month))), cols=tim.colors(length(unique(predpts$month))), stringsAsFactors=FALSE)
  predpts = merge(predpts, cols3, by.x="month", by.y="id")
  #set color scale
  vls = sort(unique(round(cols3$id)))
  vls[1] = 1
  cols4 = tim.colors(length(vls))

  plot(elev, ext=myext, ylab="Latitude", xlab="Longitude", main=species, 
       cex.lab = 1, cex.axis=1, col=gray(10:256/256))
  points(predpts$lon, predpts$lat, col=predpts$cols, pch=19, cex=0.25)
  legend("bottomleft", legend=vls, pch=22, pt.bg=cols4, pt.cex=0.75, cex=0.75, bty="n",
        col="black", title="", x.intersp=0.5, y.intersp=0.75)
}


# plot a color bar for legend
# http://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/
legend.col <- function(col, lev){
  opar <- par
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  xx <- rep(box.cx, each = 2)
  par(xpd = TRUE)
  for(i in 1:n){
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = 0.25)
  par <- opar
}

