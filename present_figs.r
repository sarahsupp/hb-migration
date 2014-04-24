# figures for AMOVEE talk 

#to describe methods for population-level migration tracking
j100 = subset(yrdat, julian == 100)

plot(j100$LONGITUDE, j100$LATITUDE, col = "red", pch = 16, add=TRUE, 
     xlim = c(-130, -65), ylim = c(25, 50), 
     xlab = "longitude", ylab = "latitude", main = "observations on day 100", cex.lab = 2, cex.axis = 2)
map("usa", add=T)
plot(hexgrid, add=T)

ID = over(SpatialPoints(j100[,c(10,9)]), hexgrid) 
names(ID) = c("JOIN_COUNT", "AREA", "PERIMETER", "BOB_", "BOB_ID", "ID", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE")
coords = cbind(j100, ID) 

t= table(as.factor(coords$POLYFID))

#Merge the count data for the hexes with all the hexes in the map
df = data.frame(POLYFID = names(t), count=as.numeric(t))
df2 = data.frame(POLYFID = unique(hexgrid$POLYFID))
df3 = merge(df2, df, all.x=TRUE)

#matches colors with the number of observations
cols = data.frame(id=c(NA,sort(unique(df3$count))), cols=tim.colors(length(unique(df3$count))), stringsAsFactors=FALSE)
df4 = merge(df3, cols, by.x="count", by.y="id")
df5 = merge(hexgrid, df4, by.x="POLYFID", by.y="POLYFID", all.x=TRUE)

#hexes with no counts are white
df5$cols = ifelse(is.na(df5$count), "white", df5$cols)
df5 = df5[order(df5$POLYFID),]
#set color scale
vls = sort(unique(round(cols$id/3)*3))
vls[1] = 1
cols2 = tim.colors(length(vls))

plot(j100$LONGITUDE, j100$LATITUDE, col = "white", pch = 16, add=TRUE, 
     xlim = c(-130, -65), ylim = c(25, 50), 
     xlab = "longitude", ylab = "latitude", main = "observations on day 100", cex.lab = 2, cex.axis = 2)
plot(hexgrid, col=df5$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1, add=T)
legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1.5, cex=1.5, bty="n",
       col="black", title="# obs", x.intersp=1, y.intersp=0.5)
map("usa", add=T)

#PLOT WHERE CENTROID WILL OCCUR

#aggregate daily info by mean centroid location
dailyHexes = count(coords, vars=c("julian", "POLYFID", "HEX_LONGITUDE", "HEX_LATITUDE"))

jdata = dailyHexes
jeffort = yreffort[which(yreffort$DAY == 100),]
jdata = merge(jdata, jeffort, by.x = "POLYFID", by.y = "POLYFID")                                           
  numcells = nrow(jdata)
  numobs = sum(jdata$freq)
  mo = as.numeric(months(100))
  day = as.numeric(days(100))
  wtmean_lon = wt.mean(jdata$HEX_LONGITUDE, jdata$freq/jdata$COUNT)
  wtmean_lat = wt.mean(jdata$HEX_LATITUDE, jdata$freq/jdata$COUNT)
  lon_sd = wt.sd(jdata$HEX_LONGITUDE, jdata$freq/jdata$COUNT)
  lat_sd = wt.sd(jdata$HEX_LATITUDE, jdata$freq/jdata$COUNT)
  meanloc = c(wtmean_lon, wtmean_lat, lon_sd, lat_sd)

df5$cols2 = df5$cols
df5$cols2[df5$cols2!="white"] <- "grey60"

plot(j100$LONGITUDE, j100$LATITUDE, col = "white", pch = 16, 
     xlim = c(-130, -65), ylim = c(25, 50), 
     xlab = "longitude", ylab = "latitude", main = "observations on day 100", cex.lab = 2, cex.axis = 2)
map("usa", add=T)
plot(hexgrid, col=df5$cols2, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1, add=T)
points(meanloc[1], meanloc[2], col = "red", pch=19, cex = 2)


#PLOT SMOOTHED LOCATION AND ERROR BARS
p100 = preds[preds$jday==100,]

plot(p100$lon, p100$lat, col = "white", pch = 16, 
     xlim = c(-130, -65), ylim = c(25, 50), 
     xlab = "longitude", ylab = "latitude", main = "observations on day 100", cex.lab = 2, cex.axis = 2)
map("usa", add=T)
points(meanloc[1], meanloc[2], col = "grey60", pch=19, cex = 2)
points(p100$lon, p100$lat, col = "red", pch=19, cex=2)

#plot entire location and error bars
preds$month = as.factor(preds$month)
cols3 = data.frame(id=c(sort(unique(preds$month))), cols=tim.colors(length(unique(preds$month))), stringsAsFactors=FALSE)
preds = merge(preds, cols3, by.x="month", by.y="id")
#set color scale
vls = sort(unique(round(cols3$id)))
vls[1] = 1
cols4 = tim.colors(length(vls))

plot(preds$lon, preds$lat, xlim = c(-130, -65), ylim = c(25, 50), 
     xlab = "longitude", ylab = "latitude", main = "summarized migration route", cex.lab = 2, cex.axis = 2)
map("usa", add=TRUE)
points(yrdat$LONGITUDE, yrdat$LATITUDE, pch=16, col = "grey60")
points(preds$lon, preds$lat, pch=19, col = preds$cols)
legend("bottomleft", legend=vls, pch=22, pt.bg=cols4, pt.cex=1.5, cex=1.5, bty="n",
       col="black", title="", x.intersp=0.5, y.intersp=0.25)

