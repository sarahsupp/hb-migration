# Aggregate the files by species or put them all together
files = list.files(path = getwd(), pattern = "hex-12_checklists_*", recursive=TRUE, full.names=TRUE)

for (f in 1:length(files)){
  data = read.table(files[f], header=TRUE, sep="\t", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  MONTH= as.numeric(months(data$DAY - 1))
  
  data = cbind(data, month)
  
  if (f == 1) {
    agg_data = data
  }
  
  else{
    agg_data = rbind(agg_data, data)
  }
  
  
  d = GroupDuplicates(data) #account for duplicate checklist data AFTER subset by species