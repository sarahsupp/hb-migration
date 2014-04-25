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

write.table(bchu, file = paste(writewd,"bchu.txt", sep=""), row.names=FALSE)
write.table(ruhu, file = paste(writewd,"ruhu.txt", sep=""), row.names=FALSE)
write.table(bthu, file = paste(writewd,"bthu.txt", sep=""), row.names=FALSE)
write.table(rthu, file = paste(writewd,"rthu.txt", sep=""), row.names=FALSE)
write.table(cahu, file = paste(writewd,"cahu.txt", sep=""), row.names=FALSE)


