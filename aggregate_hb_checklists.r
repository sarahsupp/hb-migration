# Cleans up the data needed for the hummingbird migration project. 
# Only needs to be done once with the raw files sent from FAL
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data"
setwd(wd)


# Aggregate the files by species or put them all together
files = list.files(path = getwd(), pattern = "hex-12_checklists_*", recursive=TRUE, full.names=TRUE)

for (f in 1:length(files)){
  data = read.table(files[f], header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  MONTH= as.numeric(months(data$DAY - 1))
  data = cbind(data, MONTH)
  
  if (f == 1) {
    agg_data = data
  }
  else{
    agg_data = rbind(agg_data, data)
  }
  print (paste("file", f, "is completed"))
}

# print the species names
spp = unique(agg_data$SCI_NAME)

#put migratory species in separate datafiles
bchu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Archilochus alexandri"),])
ruhu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Selasphorus rufus"),])
bthu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Selasphorus platycercus"),])
rthu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Archilochus colubris"),])
cahu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Selasphorus calliope"),])

#write the files to the folder for output
writewd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data/"

write.table(bchu, file = paste(writewd,"bchu.txt", sep=""), row.names=FALSE, sep=",")
write.table(ruhu, file = paste(writewd,"ruhu.txt", sep=""), row.names=FALSE, sep=",")
write.table(bthu, file = paste(writewd,"bthu.txt", sep=""), row.names=FALSE, sep=",")
write.table(rthu, file = paste(writewd,"rthu.txt", sep=""), row.names=FALSE, sep=",")
write.table(cahu, file = paste(writewd,"cahu.txt", sep=""), row.names=FALSE, sep=",")


#clean up the icosahedron file so it doesn't have to be brought in globally each time
hexgrid = readShapePoly(paste(getwd(), "/icosahedron_land_and_sea/cell_out2.shp", sep=""))
hexpoints = readShapePoints(paste(getwd(), "/icosahedron_land_and_sea/point_out2.shp", sep=""))
hexlonlat = data.frame(ID = hexpoints@data$global_id, LON = hexpoints@coords[,1], LAT = hexpoints@coords[,2])
hexlonlat$POLYFID = as.integer(as.character(hexlonlat$ID))

hexgrid$POLYFID = hexlonlat$POLYFID
hexgrid$LATITUDE = hexlonlat$LAT
hexgrid$LONGITUDE = hexlonlat$LON

# crop to just North America, where the migratory species occur
hexgrid = hexgrid[which(hexgrid$LATITUDE > 15 & hexgrid$LATITUDE <75 & 
                          hexgrid$LONGITUDE > -175 & hexgrid$LONGITUDE < -50),]

#write cropped file to the data folder
writePolyShape(hexgrid, "icosahedron.shp")
