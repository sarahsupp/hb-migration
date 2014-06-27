# This script will combine all movebank annotations for a species (all years, all variables) into one 
# file for analysis.
###########################################################################################
#directory containing annotations
ann_dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/downloaded_annotations/ruhu/"
###########################################################################################
#List files in ann_dir
ann_files <- list.files(ann_dir, pattern="*csv", full.names=T)

#find unique variables
unq <- unlist(strsplit(ann_files,"_"))
unq2 <- unlist(strsplit(unique(grep("*.csv", unq, value=T)),"\\."))
unq3 <- unq2[unq2 != "csv"]

#group tables by unique variables and rbind
for (i in unq3) {
  vars <- grep(i, ann_files, value=T)
  #name output tempfile
  nm <- unlist(strsplit(basename(ann_files[1]),"_"))
  out_tmp <- paste0(ann_dir,"/tmp/", nm[1], "_", nm[4],"_",i,".csv")
  
  #deal with presence and absence
  for (j in c(1:length(vars))) {
    
    if ((j == 1) && (grepl("_abs_", vars[j]))) {
      #For first table, pres/abs, then just write it out to csv
      abs <- read.csv(vars[j])
      abs$present <- 0
      #switch columns around so env vars are at the end
      abs <- abs[,c(1:4,ncol(abs),5:c(ncol(abs)-1))]
      write.csv(abs, out_tmp, quote=F, row.names=F)   
      
    } else if ((j == 1) && (grepl("_pres_", vars[j]))) {
        #For first table, pres/abs, then just write it out to csv
        abs <- read.csv(vars[j])
        abs$present <- 1
        #switch columns around so env vars are at the end
        abs <- abs[,c(1:4,ncol(abs),5:c(ncol(abs)-1))]
        write.csv(abs, out_tmp, quote=F, row.names=F)   
    
    } else if ((j > 1) && (grepl("_abs_", vars[j]))) {
        abs <- read.csv(vars[j])
        abs$present <- 0
        #switch columns around so env vars are at the end
        abs <- abs[,c(1:4,ncol(abs),5:c(ncol(abs)-1))]
        write.table(abs, out_tmp, sep=',', quote=F, row.names=F, append=T, col.names=F)
        
    } else if ((j > 1) && (grepl("_pres_", vars[j]))) {
        abs <- read.csv(vars[j])
        abs$present <- 1
        #switch columns around so env vars are at the end
        abs <- abs[,c(1:4,ncol(abs),5:c(ncol(abs)-1))]
        write.table(abs, out_tmp, sep=',', quote=F, row.names=F, append=T, col.names=F)
      
    } else {
        print("Cannot determine present or absent from file name. Check into this!")
    }#end j if
    }#end j loop
  }#end unq loop
 
#now loop over tmp files and cbind them
tmp_files <- list.files(paste0(ann_dir, "/tmp/"), pattern="*csv", full.names=T)
f_name <- paste0(ann_dir, "combined/", nm[1], "_", nm[4], "_allYears.csv")

for (tmp in c(1:length(tmp_files))) {
  if (tmp == 1) {
    #begin aggregation by writing out first file in entirety
    agg <- read.csv(tmp_files[tmp])
    #write.table(agg, f_name, sep=",", quote=F, row.names=F)
    
  } else if (tmp >1) {
    #agg <- read.csv(f_name)
    nxt <- read.csv(tmp_files[tmp], header=T)
    
    if (nrow(agg) == nrow(nxt)) {
      #combine variables into one table
      nxt.sub <- as.data.frame(nxt[,c(6:ncol(nxt))])
      names(nxt.sub) <- names(nxt[c(6:ncol(nxt))])
      agg <- cbind(agg, nxt.sub)
      #write.table(agg, f_name, sep=",", quote=F, row.names=F, append=T)

    } else {
      #skip this one until the rest of the data come in
      print(paste0("Number of rows differ between ", f_name, " and ", tmp_files[tmp]))
      print(paste0("Try again after all data have been downloaded."))
      
  }#end dimensions test if
  }#End tmp=1 if
}#end tmp file loop

write.table(agg, f_name, sep=",", quote=F, row.names=F)
