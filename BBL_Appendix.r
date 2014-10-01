# This code is for mapping North American hummingbird banding data from the Bird Banding 
# Laboratory (BBL). Recapture data is available through 2013 and was last updated 11 April 2014. 
# BBL Data is shared with permission from the USGS Bird Banding Laboratory
# BBL 1991-2012 summaries of bander data is shared with permission of USGS, 
#     emailed by Danny Bystrak, 29 Sep 2014
#     bander summary data uses the 39th parallel as the northern boundary and Louisiana, 
#     Arkansas, Missouri as the western end to define the "Southeast". 
#     The first table are all bandings all year (not just wintering birds). 
#     The second is wintering only (November through March).
#(c) 2014 Sarah Supp, modified from Marisa Lim, Stony Brook University
# Makes the figures for Appendix A for the manuscript 
# "Citizen-science data provides new insight into annual and seasonal variation in migration patterns"

# required packages
library(maps)
library(geosphere)
library(raster)
library(ggplot2)
library(maptools)
library(mapdata)
        
# set working directory
wd <- "C:/Users/sarah/Dropbox/Hummingbirds/NASA_Anusha/MarisaRAstuff/BBLrecapproj"
figpath <- "C:/Users/sarah/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/Figures"
datpath <- "/Users/sarah/Desktop/Dropbox/ActiveResearchProjects/Hummingbird_eBirdMigration/data/BBL_data/"
setwd(wd)

# import the data
BBLdat <- read.csv("11April14_BBLdata.csv", sep=",", header=T)
year = read.csv(paste(datpath, "fullyear-banding-1991-2012.csv", sep=""), header=T)
winter = read.csv(paste(datpath, "winterbanding-1991-2012.csv", sep=""), header=T)


#Plot the number of banders for the winter season
num_banders = ggplot(winter, aes(yr, num_banders)) + geom_line() + theme_classic() + 
  xlab("year") + ylab("number of registered banders")
ggsave(filename=paste(figpath,"figA4_winterbanders.png", sep=""), plot=num_banders, dpi=600, height=3, width=4)
ggplot(year, aes(yr, num_banders)) + geom_line() + theme_bw()


# First, let's clean up the data. Information that needs to be cleaned: (quick filter check of excel file shows that these need cleaning)
# 1. incorrect months of year (in ENCOUNTER_MONTH)
# 2. incorrect days of month (in ENCOUNTER_DAY)
# 3. unidentified species (in B_SPECIES_NAME)
# 4. missing lat/long information (in E_LAT_DECIMAL_DEGREES, E_LON_DECIMAL_DEGREES)
BBLdat1 <- subset(BBLdat, BBLdat$ENCOUNTER_MONTH <= 12)
BBLdat2 <- subset (BBLdat1, BBLdat1$ENCOUNTER_DAY <= 31)
BBLdat3 <- subset(BBLdat2, BBLdat2$B_SPECIES_NAME != "Unidentified Hummingbird")
BBLdat4 <- subset(BBLdat3, BBLdat3$E_LAT_DECIMAL_DEGREES != "NA")

# use only the North American migratory species
species <- c("Ruby-throated Hummingbird", "Black-chinned Hummingbird", "Broad-tailed Hummingbird", "Calliope Hummingbird", "Rufous Hummingbird")

for (s in 1:length(species)){
  
  spdat <- BBLdat4[BBLdat4$B_SPECIES_NAME == species[s],]
  bands <- unique(spdat$BAND_NUM)
  eastern <- spdat[spdat$B_LON_DECIMAL_DEGREES > -103,] #banded in the east
  western <- spdat[spdat$B_LON_DECIMAL_DEGREES <= -103,] #banded in the west
  alleast <- spdat[spdat$B_LON_DECIMAL_DEGREES > -103 & spdat$E_LON_DECIMAL_DEGREES > -103,] #banded and recaptured in the east
  
  #-------------- plot histogram of months captured outside of Western Flyway
  ggplot(alleast, aes(ENCOUNTER_MONTH)) + xlab("encounter month") + 
    geom_histogram(fill = "gray", aes = 0.5, binwidth=1) +
    theme_classic() + theme(text=element_text(size=10)) + 
    scale_x_continuous(breaks = seq(1, 12, 2), limits = c(1,12))
  ggsave(filename = paste(figpath, "/figA1_", species[s], ".pdf", sep=""), width=3, height=3, dpi=600, units="in")
                                                                                     
  #-------------- plot histogram of hb captured outside of Western Flyway
  ggplot(spdat, aes(ENCOUNTER_YEAR)) + ylab("Number birds banded") +
    xlab("Banding Year") + geom_histogram(fill = "black", binwidth=1, alpha=0.25) + 
    geom_histogram(data=eastern, aes(ENCOUNTER_YEAR), fill="red", alpha = 0.5, binwidth=1) +
    theme_classic() + theme(text=element_text(size=10), axis.text.x = element_text(angle = 60, hjust=1)) + 
    scale_x_continuous(breaks = seq(1980, 2013, 5), limits = c(1980,2013))
  ggsave(filename = paste(figpath, "/figA2_", species[s], ".pdf", sep=""), width=3, height=3, dpi=600, units="in")
  
  #-------------- plot where these eastern captures went
  pdf(file = paste(figpath, "/figA3_", species[s], ".pdf", sep=""), width = 3, height = 3)
  
  plot(NA, NA, xlim=c(-170,-50), ylim=c(15,75), xlab="", ylab="", main = "", axes = FALSE,
       bty="n")
  for(k in 1:nrow(western)){
    if(nrow(western) > 0){
    if(western[k,]$E_LON_DECIMAL_DEGREES <= -103) { color = "black"}
    else { color = "indianred" }
    inter <- gcIntermediate(c(western[k,]$B_LON_DECIMAL_DEGREES, western[k,]$B_LAT_DECIMAL_DEGREES),
                            c(western[k,]$E_LON_DECIMAL_DEGREES, western[k,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col=color, lwd=0.5)
  }}
  for(j in 1:nrow(eastern)){
    if(eastern[j,]$E_LON_DECIMAL_DEGREES > -103) { color = "cadetblue"}
    else { color = "indianred" }
    inter <- gcIntermediate(c(eastern[j,]$B_LON_DECIMAL_DEGREES, eastern[j,]$B_LAT_DECIMAL_DEGREES),
                            c(eastern[j,]$E_LON_DECIMAL_DEGREES, eastern[j,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col=color, lwd=0.5)
  }
  map("worldHires", c("usa", "canada", "mexico"), add=TRUE, cex = 0.5)
  dev.off()
}




