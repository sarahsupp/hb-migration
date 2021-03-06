# Cleans up the data needed for the hummingbird migration project. 
# Only needs to be done once with the raw files sent from FAL
# files in "eBird_checklists_2008-2014
# SRS 25 Feb 2015

library(tools)
library(maptools)

#path to files
filepath = "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/"
writepath = "/home/sarah/Dropbox/Hummingbirds/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"
hexpath = "/home/sarah/Dropbox/CompletedResearchProjects/Hummingbird_eBirdMigration/"
gitpath = "/home/sarah/Documents/GitHub/hb-migration/"

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


ID_windows = function(yeardat, spring, peak, fall, timewindow){
  # This function separates each year into a given number of time windows, where window 2 (w2) begins 
  # on the date for the onset of spring migration and the last window (wn) ends in the time frame that 
  # follows the time frame containing the date for the end of fall migration (wn-1)
    # yeardat: a data.frame subsetted for a single year
    # spring: date for onset of spring migration
    # peak: date for peak latitude for the population
    # fall: date for end of fall migration
    # timewindow: numeric value indicating the time window used to aggregate observation data (e.g. 7 would define a week)
  
  yeardat$increment = 0
  yeardat$window = 0
  yeardat$season="winter"
  
  if(fall < spring | fall < peak) {
    print(paste0("ERROR: migration dates do not make sense, check data ", yeardat$SCI_NAME[1], yeardat$YEAR[1]))
    return(yeardat)
  }
  
  start <-  spring - timewindow
  n.windows <- floor(((fall - start)/timewindow) + 2)
  mid.window <- ceiling((peak - start)/timewindow)
  end <- start + n.windows * timewindow

  yeardat[yeardat$DAY >= start & yeardat$DAY < end,]$increment <- yeardat[yeardat$DAY >= start & yeardat$DAY < end,]$DAY-start + 1
  yeardat[yeardat$DAY >= start & yeardat$DAY < end,]$window <- ceiling(yeardat[yeardat$DAY >= start & yeardat$DAY < end,]$increment/timewindow)
  yeardat[yeardat$DAY >= start & yeardat$DAY < end & yeardat$window < mid.window,]$season <- "spring"
  yeardat[yeardat$DAY >= start & yeardat$DAY < end & yeardat$window > mid.window,]$season <- "fall"
  yeardat[yeardat$DAY >= start & yeardat$DAY < end & yeardat$window == mid.window,]$season <- "breeding"
  
  return(yeardat)
}

#--------------------------------------------------------------------------------------------------------
#------------------------------------------------ AGGREGATE THE FILES -----------------------------------
#--------------------------------------------------------------------------------------------------------

files = list.files(path = filepath, pattern = "eBird_checklists_*", recursive=FALSE, full.names=TRUE)

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


#----------------------------------- GET DESCRIPTIVE DATA AND MIGRATION DATES FOR EACH SPECIES
# read in summary of effort data (Number of eBird checklists submitted per day per year), and hexgrid for mapping
effort = read.table(paste0(writepath, "eBird_checklists_effort.txt"), header=TRUE, as.is=TRUE)
hexgrid = readShapePoly(paste(hexpath, "/data/icosahedron_land_and_sea/icosahedron.shp", sep="")) #hex with land and sea, cropped to North America

#make a list of the file names to loop through
files = list.files(path=writepath, pattern = "*hu08-14.txt")

for (f in 1:length(files)){
  
  require(ggmap)
  require(ggplot2)
  require(plyr)
  require(reshape2)
  require(Rmisc)
  require(sp)
  require(raster)
  require(segmented)
  source(paste(gitpath, "/Project-1/migration-fxns.r", sep=""))
  
  humdat = read.table(paste0(writepath, files[f]), header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  
  #make the column names look nicer
  names(humdat) = gsub('.{1}$', '', substring(names(humdat),3))
  #make sure month is read as an ordered factor
  humdat$MONTH = factor(humdat$MONTH, levels=c(1:12), ordered=TRUE)
  
  #grab species name for setting directory paths and naming figures, identify years
  species = gsub("\"", "", humdat$SCI_NAME[1], fixed=TRUE) 
  spcode = gsub("08-14.txt", "", files[f])
  years = sort(unique(humdat$YEAR))
  
  #plot number of records by week and year, save frequency of obs to txt file
  freqobs = ggplot(humdat, aes(DAY)) + geom_histogram(binwidth=1) + theme_bw() + facet_wrap(~YEAR)
  ggsave(file=paste(writepath, spcode, "_obs_by_year.pdf", sep=""))
  
  yeartable = PlotRecords(humdat$YEAR, species)
  write.table(yeartable, file = paste(writepath, spcode,"_numobs.txt",sep=""), row.names=FALSE)
  
  for (y in 1:length(years)){
    yrdat = humdat[which(humdat$YEAR == years[y]),]
    yreffort = effort[which(effort$YEAR == years[y]),]
    
    meanlocs = AlternateMeanLocs(yrdat, species, hexgrid, yreffort)  #daily weighted mean location
    preds = EstimateDailyLocs(meanlocs) #predict daily location along a smoothing line
      preds$year = years[y]
    migration = round(Est3MigrationDates(meanlocs)) #gam to estimate starting points for segmentation from the mean loc latitude data
    dist = DailyTravel(preds, 4, 5, species, years[y], migration) #Great Circle distances traveled between predicted daily locations
    speed = MigrationSpeed(dist, migration) #estimate migration speed for spring and fall
    
    if (y == 1){
      pred_data = preds
      mig_data = data.frame("spr_begin" = migration[[1]], "peak_lat" = migration[[2]], 
                            "fal_end" = migration[[3]], "spr_spd" = speed[1], "fal_spd" = speed[2],
                            "species" = species, "year" = years[y])      
    }
    else{
      pred_data = rbind(pred_data, preds)
      dates = c(migration[[1]], migration[[2]], migration[[3]], speed[1], speed[2], species, years[y])
      mig_data = rbind(mig_data, dates)
      
      if (y == length(years)){
        write.table(pred_data, file = paste0(writepath, spcode, "_preds.txt"), append=FALSE,row.names=FALSE, sep=",")
        write.table(mig_data, file = paste0(writepath, spcode, "_migration.txt"), append=FALSE, row.names=FALSE, sep=",")
      }
    }

    
    #Collect data for WEST of 103 W Longitude
    if (spcode %in% c("bchu", "bthu", "ruhu", "cahu")) {
      #only use data west of the 103rd meridian (western flyway), where species had > 5000 obs 
      west_yrdat = yrdat[which(yrdat$LONGITUDE <= -103),]
      west_meanlocs = AlternateMeanLocs(west_yrdat, spcode, hexgrid, yreffort)
      west_preds = EstimateDailyLocs(west_meanlocs)
        west_preds$year = years[y]
      west_migration = round(Est3MigrationDates(west_meanlocs))
      west_dist = DailyTravel(west_preds, 4, 5, spcode, years[y], west_migration)
      west_speed = MigrationSpeed(west_dist, west_migration)
      #append migration dates (repated across the year for later use)
      west_yrdat$spr = west_migration[[1]]
      west_yrdat$peak = west_migration[[2]]
      west_yrdat$fal = west_migration[[3]]
      
      if (y == 1){
        west_pred_data = west_preds
        west_mig_data = data.frame("spr_begin" = west_migration[[1]], "peak_lat" = west_migration[[2]], 
                              "fal_end" = west_migration[[3]], "spr_spd" = west_speed[1], "fal_spd" = west_speed[2],
                              "species" = species, "year" = years[y])    
        west_humdat = west_yrdat
        
      }
      else{
        west_pred_data = rbind(west_pred_data, west_preds)
        west_dates = c(west_migration[[1]], west_migration[[2]], west_migration[[3]], west_speed[1], west_speed[2], species, years[y])
        west_mig_data = rbind(west_mig_data, west_dates)
        west_humdat = rbind(west_humdat, west_yrdat)
        
        if (y == length(years)){
          write.table(west_pred_data, file = paste0(writepath, spcode, "_preds_west.txt"), append=FALSE,row.names=FALSE, sep=",")
          write.table(west_mig_data, file = paste0(writepath, spcode, "_migration_west.txt"), append=FALSE, row.names=FALSE, sep=",")
          write.table(west_humdat, file = paste0(writepath, spcode, "_humdat_west.txt"), append=FALSE, row.names=FALSE, sep=",")
          
        }
      }
    }
    
    #Collect data for EAST of 103 W Longitude (eastern + central flyway), where species had > 5000 obs 
    if (spcode %in% c("bchu", "rthu")) {
      #only use data west of the 103rd meridian (western flyway)
      east_yrdat = yrdat[which(yrdat$LONGITUDE > -103),]
      east_meanlocs = AlternateMeanLocs(east_yrdat, spcode, hexgrid, yreffort)
      east_preds = EstimateDailyLocs(east_meanlocs)
        east_preds$year = years[y]
      east_migration = round(Est3MigrationDates(east_meanlocs))
      east_dist = DailyTravel(east_preds, 4, 5, spcode, years[y], east_migration)
      east_speed = MigrationSpeed(east_dist, east_migration)
      #append migration dates (repated across the year for later use)
      east_yrdat$spr = east_migration[[1]]
      east_yrdat$peak = east_migration[[2]]
      east_yrdat$fal = east_migration[[3]]
     
      if (y == 1){
        east_pred_data = east_preds
        east_mig_data = data.frame("spr_begin" = east_migration[[1]], "peak_lat" = east_migration[[2]], 
                                   "fal_end" = east_migration[[3]], "spr_spd" = east_speed[1], "fal_spd" = east_speed[2],
                                   "species" = species, "year" = years[y])    
        east_humdat = east_yrdat
      }
      else{
        east_pred_data = rbind(east_pred_data, east_preds)
        east_dates = c(east_migration[[1]], east_migration[[2]], east_migration[[3]], east_speed[1], east_speed[2], species, years[y])
        east_mig_data = rbind(east_mig_data, east_dates)
        east_humdat = rbind(east_humdat, east_yrdat)
        
        
        if (y == length(years)){
          write.table(east_pred_data, file = paste0(writepath, spcode, "_preds_east.txt"), append=FALSE,row.names=FALSE, sep=",")
          write.table(east_mig_data, file = paste0(writepath, spcode, "_migration_east.txt"), append=FALSE, row.names=FALSE, sep=",")
          write.table(east_humdat, file = paste0(writepath, spcode, "_humdat_east.txt"), append=FALSE, row.names=FALSE, sep=",")
        }
      }
    }
  }
}



#----------------------------------- ASSIGN TIME WINDOWS FOR EACH SPECIES
# Use the humdat files made above and annoate with time windows. 
# Keeping it separate allows you to skip repeating the time-consuming step of the gamms and estimating migration dates

#make a list of the file names to loop through
files = list.files(path=writepath, pattern = "_humdat_.*\\.*st.txt$")

#specify time frame (number of days) to group observations
timeframe = 7

for (f in 1:length(files)){
  
  #read in humdat files created during migration paths
  humdat = read.table(paste0(writepath, files[f]), header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")
  names(humdat) = gsub('.{1}$', '', substring(names(humdat),3))  #make the column names look nicer
 
  years = unique(humdat$YEAR)
  
  for (y in 1:length(years)){
    yrdat = humdat[humdat$YEAR == years[y],]
    #assign time frames for making alpha hulls (the last value is the desired time frame)
    humdat_win = ID_windows(yrdat, yrdat[1,]$spr, yrdat[1,]$peak, yrdat[1,]$fal, timeframe)

    if (y == 1){
      humwindows = humdat_win
    }
    else{
      humwindows = rbind(humwindows, humdat_win)
      
      if (y == length(years)){
        write.table(humwindows, file = paste0(writepath, "t", timeframe, "_", files[f]), append=FALSE, row.names=FALSE, sep=",")
      }
    }
  }
}

