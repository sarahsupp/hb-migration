# This code is for the eBird observer effort tables for the migration project and to plot the results
# In this case, effort represents the total number of submitted checklists per hex cell
# The data can be aggregated by day, year, or hex cell. 
# Data was obtained from FAL, 2014/02/21
# NOTE: This is only complete through 30 Nov 2013
#(c) 2014 Sarah Supp


#set working directory
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data/"
setwd(wd)

#read in checklist effort data from FAL
eft = read.table("icosohedron_ebd_effort_2004-2013.txt", sep=" ", header=T, as.is=TRUE)

#convert to year and julian date & bind into new dataframe
YEAR = sapply(eft$OBSERVATION.DATE,function(x){
  as.numeric(substring(x, 1, 4))
})

DAY = sapply(eft$OBSERVATION.DATE, function(x){
  julian(as.numeric(substring(x, 6, 7)), as.numeric(substring(x, 9, 10)), as.numeric(substring(x, 1, 4)), 
         origin. = c(1, 1, as.numeric(substring(x, 1, 4)))) + 1
})

eft_new = data.frame(YEAR, DAY, "POLYFID" = eft$POLYFID, "COUNT" = eft$COUNT)


#aggregate by daily effort in each year
effort_day = aggregate(eft_new$COUNT, list(year=eft_new$YEAR, day=eft_new$DAY), sum)
names(effort_day)[3] = "effort"

#aggregate by total effort in each hex cell in each year
effort_cell = aggregate(effort_all$COUNT, list(year=effort_all$YEAR, POLYFID=effort_all$POLYFID), sum)
names(effort_cell)[3] = "effort"

effort_cellTot = aggregate(effort_all$COUNT, list(POLYFID=effort_all$POLYFID), sum)
names(effort_cellTot)[2] = "effort"

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
  pdf(paste("TotalEffort_new", years[y], ".pdf", sep=""))
  plot(hexgrid, col=df$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1)
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


#------------ Plot the total effort for 2004-2013
eft = effort_cellTot
df = merge(hexgrid, eft, by.all.x=TRUE)
#matches colors with the number of observations
cols = data.frame(id=c(NA,sort(unique(df$effort))), cols=tim.colors(length(unique(df$effort))), 
                  stringsAsFactors=FALSE)
df = merge(df,cols, by.x="effort", by.y="id")
#hexes with no counts are white
df$cols = ifelse(is.na(df$effort), "white", df$cols)
df = df[order(df$POLYFID),]

vls = sort(unique(round(cols$id/2000)*2000))
vls[1] = 1
cols2 = tim.colors(length(vls))

#make a map with hexes colored by the number of times the species was observed in a given hex
pdf("TotalEffort_new_2004-2013.pdf")
plot(hexgrid, col=df$cols, border="gray10", lwd=0.25, xlim=c(-170,-50), ylim=c(10,80), las=1)
axis(side=1)
axis(side=2, las=1)
box()
mtext("Longitude", side=1, cex=1.4, line=2.5)
mtext("Latitude", side=2, cex=1.4, line=2.5)
mtext("Number of eBird checklists 2004-2013", side = 3, cex = 1.4, line=1)
## legend
legend("bottomleft", legend=vls, pch=22, pt.bg=cols2, pt.cex=1, cex=0.75, bty="n",
       col="black", title="Number of checklists", x.intersp=1, y.intersp=0.5)
dev.off()

#------------ Write datasets to use later
write.table(effort_all, file = "cell_effort_new.txt", row.names=FALSE)


