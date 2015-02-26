# Cleans up the data needed for the hummingbird migration project. 
# Only needs to be done once with the raw files sent from FAL
# files in "eBird_checklists_2008-2014
# SRS 25 Feb 2015

library(tools)

#path to files
filepath = "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/"
writepath = "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"


#----------------------------------------------------FUNCTIONS
GroupDuplicates = function(humdat) { 
  #gets rid of duplicate records that are part of the same group. 
  #Takes only the first record to move to analysis, keeps all single records that are not part of a group.
  gid = sort(unique(humdat$GROUP_ID))
  
  #open start a new dataframe with same columns as the main dataframe
  keep = humdat[1,]
  out = 0
  
  for (g in 1:length(gid)){
    out = out + 1
    tmp = humdat[which(humdat$GROUP_ID == gid[g]),]
    #record the first line of the data (assume the group has the same information)
    if (nrow(tmp) == 1) { keep[out,] = tmp }
    else{ keep[out,] = tmp[1,] }
  }
  
  keepnongroup = humdat[which(is.na(humdat$GROUP_ID)),]
  keep = rbind(keep, keepnongroup)
  return(keep)
}


#------------------------------------------------ AGGREGATE THE FILES
files = list.files(path = filepath, pattern = "eBird_checklists_*", recursive=TRUE, full.names=TRUE)

for (f in 1:length(files)){
  data = read.table(files[f], header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  #really ugly regex to pull month from filename
  MONTH = as.numeric(sub("_20[0-9][0-9]", "", sub("eBird_checklists_", "", basename(file_path_sans_ext(files[f])))))
  MONTH = rep(MONTH, nrow(data))
  data = cbind(data, MONTH)
  
  if (f == 1) { agg_data = data }
  else{ agg_data = rbind(agg_data, data) }
  print (paste("file", f, "is completed:", basename(file_path_sans_ext(files[f]))))
}


# print the species names
unique(agg_data$SCI_NAME)

#put migratory species in separate datafiles
bchu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Archilochus alexandri"),])
ruhu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Selasphorus rufus"),])
bthu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Selasphorus platycercus"),])
rthu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Archilochus colubris"),])
cahu = GroupDuplicates(agg_data[which(agg_data$SCI_NAME == "Selasphorus calliope"),])

print(paste0("There are ", nrow(cahu), " CAHU, ", nrow(rthu), " RTHU, ", nrow(bthu), " BTHU, ", nrow(ruhu), " RUHU, and ",
             nrow(bchu), " BCHU, and ", nrow(cahu) + nrow(rthu) + nrow(bthu) + nrow(ruhu) + nrow(bchu), " total observations"))

#write the files to the folder for output
write.table(bchu, file = paste(writepath,"bchu08-14.txt", sep=""), row.names=FALSE, sep=",")
write.table(ruhu, file = paste(writepath,"ruhu08-14.txt", sep=""), row.names=FALSE, sep=",")
write.table(bthu, file = paste(writepath,"bthu08-14.txt", sep=""), row.names=FALSE, sep=",")
write.table(rthu, file = paste(writepath,"rthu08-14.txt", sep=""), row.names=FALSE, sep=",")
write.table(cahu, file = paste(writepath,"cahu08-14.txt", sep=""), row.names=FALSE, sep=",")
