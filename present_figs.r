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

vls = sort(unique(round(cols$id/5)*5))
vls[1] = 1
cols2 = tim.colors(length(vls))

plot(j100$LONGITUDE, j100$LATITUDE, col = "white", pch = 16, add=TRUE, 
     xlim = c(-130, -65), ylim = c(25, 50), 
     xlab = "longitude", ylab = "latitude", main = "observations on day 100", cex.lab = 2, cex.axis = 2)
plot(hexgrid, col=df5$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1, add=T)
legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1.5, cex=1.5, bty="n",
       col="black", title="# obs", x.intersp=1, y.intersp=0.5)
map("usa", add=T)

PLOT WHERE CENTROID WILL OCCUR, 

PLOT SMOOTHED LOCATION AND ERROR BARS

