# Purpose: This R script uses curl to submit an annotation request to movebank.
# The arguments required include: url=url to which to submit request; xy=a csv file containing a list
# of xy coordinates for annotation (see "/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/test_calliope.csv")
# for an example of how it must be structured); xml=an xml file containing the request parameters (i.e. sensor, product, 
# interpolation etc. See /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/NDVI_EVI.xml for
# example); un=movebank username; pw=movebank password; outdir=output directory for returned xml document with access key.
#
# Author: Tina Cormier
# Date: April 21, 2014
# Status: In Progress



library(RPostgreSQL)
library(XML)


# options(echo=TRUE)
# args <- commandArgs(trailingOnly=TRUE)
# 
# url <- args[1]
# xy <- args[2]
# xml <- args[3]
# un <- args[4]
# pw <- args[5]
# outdir <- args[6]

###############################
# #hardcoded for testing:
#
#address to submit request
url <- "http://www.bioinfo.mpg.de/orn-gateway/request-annotation-xml.jsp"

#CSV file list of track csvs that contain xy locations of points to be annotated. Each track csv must 
#have the following columns with no spaces: Timestamp, location-long, location-lat, height-above-ellipsoid (optional)
xy.csv <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/track_csvs/submit_lists/bchu_lag0_allyrs.csv"


#CSV file list of XMLs containing movebank request details (see hb_movebank_createXML.py)
#xml <- "/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/env_request_xmls/Quality_NDVI_EVI.xml"
xml.csv <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/env_request_xmls/submit_lists/submit_vars_20140626.csv"

#user name
un <- "tcormier"

#password
pw <- "nohg3ITh"

#outdir
outdir <- "C:/Share/tcormier/hummingbirds/migration_study/movebank/submitted_requests/"
#outfile <- "/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests/testR.xml"
##############################

#SKIP DB stuff for now - can't get it to work on arctic.
#connect to db
# Connect to database
# con <- dbConnect(drv="PostgreSQL", port="5433",host="fusion", dbname="tcormier")
#dbListTables(con)

#read in submit lists:
xy.list <- as.vector(read.csv(xy.csv, header=F, as.is=T)[,1])
xml.list <- as.vector(read.csv(xml.csv, header=F, as.is=T)[,1])

#FOR NOW - until I can get database working from arctic:

#Now, loop over each xy list and submit a request for each xml in xml.list
for (xy in xy.list) {
  tbl <- read.csv(xy)
  
  for (xml in xml.list) {
    #Curl command will return an xml with the access key in it - filepath of that returned file:
    now <- format(Sys.time(), format="%Y%m%d_%H%M%S")
    ret.xml <- paste(outdir, "/", unlist(strsplit(basename(xy), "\\."))[1], "_", unlist(strsplit(basename(xml), "\\."))[1], "_", now, ".xml", sep="")
    
    #system call to curl (was able to get this to work somehow when I couldn't get Rcurl to work! - just need to plug ahead for now.)
    curl.cmd <- paste("curl -o ", ret.xml, " -F \"request=@", xml, "\" -F \"tracks=@", xy, "\" -F \"login=", un, "\" -F \"password=", pw, "\" ", url, sep="")
    system(curl.cmd)
    
    #check to see if ret.xml exists (it should if the curl system call did not produce an error).
    #If it exists, write access key to db
    if (file.exists(ret.xml)) {
      #parse returned xml file
      doc <- xmlTreeParse(ret.xml)
      r <- xmlRoot(doc)
      ak <- xmlAttrs(r)[names(xmlAttrs(r))=="accessKey"]
      ak.status <- xmlAttrs(r)[names(xmlAttrs(r))=="status"]
      
      #FIX this to non-test table when script is stable and production-ready
      akl.query <- sprintf("INSERT INTO test_access_key_list (access_key) VALUES (%s)", ak)
      insert.ak <- dbGetQuery(con, akl.query)
      
      #talk to Jesse about making this a real date-time; also look up how to more easily write a query to avoid
      #tedious issues with quotes fixed by using sprintf (can use this to convert python xml generator code over to R?)
      status.query <- sprintf("INSERT INTO test_access_key_lut (access_key, date_time, tracks, xml, status) VALUES (%s,'%s','%s','%s','%s')", ak,now,xy,xml,ak.status)
      insert.ak <- dbGetQuery(con, status.query)
    }# end ret.xml if
  }
  
}#end xy loop




