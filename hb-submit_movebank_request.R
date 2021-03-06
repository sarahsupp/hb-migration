# Purpose: This R script uses curl to submit an annotation request to movebank.
# The arguments required include: url=url to which to submit request; xy=a csv file containing a list
# of xy coordinates for annotation (see "/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/test_calliope.csv")
# for an example of how it must be structured); xml=an xml file containing the request parameters (i.e. sensor, product, 
# interpolation etc. See /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/NDVI_EVI.xml for
# example); un=movebank username; pw=movebank password; req.outdir=output directory for returned xml document with access key;
# ann.outdir=annotation output directory - will be written to db.
#
# Author: Tina Cormier
# Date: April 21, 2014
# Status: In Progress

library(RPostgreSQL)
library(XML)
source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/hb_RS_functions.R")

# options(echo=TRUE)
# args <- commandArgs(trailingOnly=TRUE)
# 
# url <- args[1]
# xy <- args[2]
# xml <- args[3]
# un <- args[4]
# pw <- args[5]
# outdir <- args[6]

###########################################################################################
# #hardcoded for testing:
#
spp <- "ruhu"
#address to submit request
#url <- ""

#CSV file list of track csvs that contain xy locations of points to be annotated. Each track csv must 
#have the following columns with no spaces: Timestamp, location-long, location-lat, height-above-ellipsoid (optional)
#xy.csv <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/track_csvs/submit_lists/bchu_lag0_allyrs.csv"
xy.csv <- paste0("/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/track_csvs/submit_lists/",spp,"_lag0_allyrs.csv")

#CSV file list of XMLs containing movebank request details (see hb_movebank_createXML.py)
#xml.csv <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/env_request_xmls/submit_lists/submit_vars_20140626.csv"
xml.csv <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/env_request_xmls/submit_lists/veg_pixel_reliability_20140814.csv"
  
#user name
#un <- ""

#password
#pw <- ""

#request outdir
#req.outdir <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/submitted_requests/"
req.outdir <- paste0("/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/submitted_requests/",spp,"/")

#annotation outdir
ann.outdir <- paste0("/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/downloaded_annotations/",spp,"/")
##########################################################################################

#SKIP DB stuff on arctic for now - can't get it to work on arctic - works on the mac though :)
#connect to db
# Connect to database
con <- dbConnect(drv="PostgreSQL", port="5433",host="fusion", dbname="tcormier")
#dbListTables(con)

#read in submit lists:
xy.list <- as.vector(read.csv(xy.csv, header=F, as.is=T)[,1])
xml.list <- as.vector(read.csv(xml.csv, header=F, as.is=T)[,1])

#FOR NOW - until I can get database working from arctic - do this on mac:

#Now, loop over each xy list and submit a request for each xml in xml.list
for (xy in xy.list) {
  tbl <- read.csv(xy)
  
  for (xml in xml.list) {
    #submit request - returns two objects in a list
    mb.return <- submitMB(url,un,pw,xy,xml,req.outdir)
    ret.xml <- mb.return[[1]]
    db.now <- format(mb.return[[2]], format='%Y-%m-%d %H:%M:%S')
    
    #check to see if ret.xml exists (it should if the curl system call did not produce an error).
    #If it exists, write access key to db
    if (file.exists(ret.xml)) {
      #parse returned xml file
      doc <- xmlTreeParse(ret.xml)
      r <- xmlRoot(doc)
      ak <- xmlAttrs(r)[names(xmlAttrs(r))=="accessKey"]
      ak.status <- xmlAttrs(r)[names(xmlAttrs(r))=="status"]
      ann_name <- paste0(ann.outdir,unlist(strsplit(basename(xy), "\\."))[1], "_", unlist(strsplit(basename(xml), "\\."))[1], ".csv")
      spp <- unlist(strsplit(basename(ann_name), "_"))[1]
      
      #talk to Jesse about making this a real date-time (DONE)
      status.query <- sprintf("INSERT INTO access_key_lut (access_key, date_time, tracks, xml, status, ann_name,spp) VALUES ('%s','%s','%s','%s','%s','%s','%s')", ak,db.now,xy,xml,ak.status,ann_name,spp)
      insert.ak <- dbGetQuery(con, status.query)
    } else {
        print(paste0(ret.xml, " does not exist - check to see if request was submitted.")) 
    }# end ret.xml if
  }#end xml for loop
}#end xy loop




