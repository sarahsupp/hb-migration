# This code is to bind together the effort tables for the migration project and plot the results
# In this case, effort represents the total number of submitted checklists per hex cell
# The data can be aggregated by day, year, or hex cell. 
# Data was obtained from FAL
#(c) 2014 Sarah Supp


#set working directory
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data/effort_2004-2013"
setwd(wd)

# read in the north america equal area hex grid map (FAL)
hexgrid = readShapePoly("C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/terr_4h6/terr_4h6.shp") #hex map
# crop to just North America, where the migratory species occur
hexgrid = hexgrid[which(hexgrid$LATITUDE > 10 & hexgrid$LATITUDE <80 & 
                          hexgrid$LONGITUDE > -178 & hexgrid$LONGITUDE < -50),]

#hex grid map
pdf("HexgridMap.pdf")
plot(hexgrid, xlim=c(-170,-50), ylim=c(10,80), col="lightblue", lwd=0.25, border="gray10")
axis(side=1)
axis(side=2, las=1)
box()
mtext("Longitude", side=1, cex=1.4, line=2.5)
mtext("Latitude", side=2, cex=1.4, line=2.5)
dev.off()

# read in effort data, number of checklists submitted per hex cell
files = list.files(pattern = "*.txt")

#import and rbind all the files together

for (f in 1:length(files)){
  
  effort = read.table(files[f], sep="\t", header=T, as.is=TRUE)
  
  if (f == 1) {
    effort_all = effort
  }
  else{
   effort_all = rbind(effort_all, effort) 
  }
}

#aggregate by daily effort in each year
effort_day <- aggregate(effort_all$COUNT, list(year=effort_all$YEAR, day=effort_all$DAY), sum)
names(effort_day)[[3]] <- "effort"

#aggregate by total effort in each hex cell
effort_cell <- aggregate(effort_all$COUNT, list(year=effort_all$YEAR, POLYFID=effort_all$POLYFID), sum)
names(effort_cell)[[3]] <- "effort"

years = c(2004:2013)

for (y in 1:length(years)){
  eft = effort_cell[which(effort_cell$year == years[y]),c(2,3)]
  df = merge(hexgrid, eft, by.all.x=TRUE)
  #matches colors with the number of observations
  cols = data.frame(id=c(NA,sort(unique(df$effort))), cols=tim.colors(length(unique(df$effort))), 
                    stringsAsFactors=FALSE)
  df = merge(df,cols, by.x="effort", by.y="id")
  #hexes with no counts are white
  df$cols = ifelse(is.na(df$effort), "white", df$cols)
  df = df[order(df$POLYFID),]
  
  vls = sort(unique(round(cols$id/500)*500))
  vls[1] = 1
  cols2 = tim.colors(length(vls))
  
  #make a map with hexes colored by the number of times the species was observed in a given hex
  pdf(paste("TotalEffort", years[y], ".pdf", sep=""))
  plot(hexdat, col=df$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1)
  axis(side=1)
  axis(side=2, las=1)
  box()
  mtext("Longitude", side=1, cex=1.4, line=2.5)
  mtext("Latitude", side=2, cex=1.4, line=2.5)
  mtext(paste("Number of eBird checklists - ", years[y], sep=""), side = 3, cex = 1.4, line=1)
  ## legend
  legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1, cex=0.75, bty="n",
         col="black", title="Number of checklists", x.intersp=1, y.intersp=0.5)
  dev.off()
}

write.table(effort_all, file = "cell_effort.txt", row.names=FALSE)
