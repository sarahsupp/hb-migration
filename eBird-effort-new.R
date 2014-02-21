# This code is for the eBird observer effort tables for the migration project and to plot the results
# In this case, effort represents the total number of submitted checklists per hex cell
# The data can be aggregated by day, year, or hex cell. 
# Data was obtained from FAL, 2014/02/21
#(c) 2014 Sarah Supp


#set working directory
wd = "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data/"
setwd(wd)

#read in checklist effort data from FAL
eft = read.table("icosohedron_ebd_effort_2004-2013.txt", sep=" ", header=T, as.is=TRUE)

#convert to year and julian date
YEAR = sapply(eft$OBSERVATION.DATE,function(x){
  as.numeric(substring(x, 1, 4))
})

DAY = sapply(eft$OBSERVATION.DATE, function(x){
  julian(as.numeric(substring(x, 6, 7)), as.numeric(substring(x, 9, 10)), as.numeric(substring(x, 1, 4)), 
         origin. = c(1, 1, as.numeric(substring(x, 1, 4)))) + 1
})


posixdates = strptime(eft$OBSERVATION.DATE, format = "%Y-%D-%M")

posix.dates <- strptime( x=paste(2009,"-",jdate.ttt.seq, sep=""),"%Y-%j")
date.string <- paste(months(posix.dates, abbreviate=T), posix.dates$mday,sep="_")