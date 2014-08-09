library(RPostgreSQL)
library(RCurl)
library(XML)

#connect to db
# Connect to database
con <- dbConnect(drv="PostgreSQL", port="5433",host="fusion", dbname="tcormier")

#select requests that are not "downloaded"
#FIX this to non-test table when script is stable and production-ready
#options(digits=22)
status.query <- "SELECT * FROM access_key_lut WHERE status != 'delivered' AND status != 'failed';"
check.list <- dbGetQuery(con,status.query)

#Now check status of these acess keys. If available, download it and record to db.
for (i in check.list$access_key) {
  #get current status from db
  st.query <- sprintf("SELECT status FROM access_key_lut WHERE access_key = '%s';",i)
  st.db <- dbGetQuery(con,st.query)
  
  #URL where you can check the status of your request
  url <- paste0("http://www.bioinfo.mpg.de/orn-gateway/request-info-xml.jsp?access-key=",as.character(i))
  
  #query name of output file from db
  ann_file.query <- sprintf("SELECT ann_name FROM access_key_lut WHERE access_key='%s';",i)
  ann_file <- as.vector(dbGetQuery(con,ann_file.query)[,1])
  ann_name <- unlist(strsplit(basename(ann_file), "\\."))[1]
  
  #This is way more difficult than it needs to be - don't think I need to write out a file to get at the status, 
  #but it's currently the only way I can get it to work. Ugh, web stuff! :)
  st.xmlfile <- paste0("/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/status_xmls/",ann_name,"_", i, ".xml")
  status.cmd <- sprintf("curl -o %s %s", st.xmlfile,url)
  st.xml <- system(status.cmd)

  #parse resulting xml to get status
  #parse returned xml file
  doc <- xmlTreeParse(st.xmlfile)
  r <- xmlRoot(doc)
  ak.status <- xmlAttrs(r)[names(xmlAttrs(r))=="status"]
  
  #if movebank status is "available", download annotation.
  if (ak.status == "available") {
    #update status in DB
    up.query <- sprintf("UPDATE access_key_lut SET status = '%s' WHERE access_key = '%s'", ak.status, i)
    up.ak <- dbGetQuery(con, up.query)
    
    #downlod file
    dl_url <- sprintf("http://www.bioinfo.mpg.de/orn-gateway/download.jsp?access-key=%s",i)
    dl.cmd <- sprintf("curl -o %s %s", ann_file, dl_url)
    dl.get <- system(dl.cmd)
    
    #check movebank status again (should be "delivered")
    #This is way more difficult than it needs to be - don't think I need to write out a file to get at the status, 
    #but it's currently the only way I can get it to work. Ugh, web stuff! :)
    st.xmlfile <- paste0("/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/status_xmls/",ann_name,"_", i, ".xml")
    status.cmd <- sprintf("curl -o %s %s", st.xmlfile,url)
    st.xml <- system(status.cmd)
    
    #parse resulting xml to get status
    #parse returned xml file
    doc <- xmlTreeParse(st.xmlfile)
    r <- xmlRoot(doc)
    ak.status <- xmlAttrs(r)[names(xmlAttrs(r))=="status"]
    
    if (file.exists(ann_file) && ak.status == "delivered"){
      #update status in DB
      up.query <- sprintf("UPDATE access_key_lut SET status = '%s' WHERE access_key = '%s'", ak.status, i)
      up.ak <- dbGetQuery(con, up.query)
    } else {
      print(paste0("conflict with access key ", i, ". Check if ", ann_file, " was downloaded!"))
      up.query <- sprintf("UPDATE test_access_key_lut SET status = 'download error' WHERE access_key = '%s'", i)
      up.ak <- dbGetQuery(con, up.query)
    }#end db update to delivered "if"
    
  #Is db status different from movebank status (i.e., has status changed since last check?)  
  } else if (ak.status != st.db){
      up.query <- sprintf("UPDATE access_key_lut SET status = '%s' WHERE access_key = '%s'", ak.status, i)
      up.ak <- dbGetQuery(con, up.query)
  }#else do nothing - end if status is available if statement
}# end access key loop


