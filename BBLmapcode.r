# This code is for mapping North American hummingbird banding data from the Bird Banding 
# Laboratory (BBL). Recapture data is available through 2013 and was last updated 11 April 2014. 
# Miscellaneous old code at end (to be cleaned/deleted)
#(c) 2013 Marisa Lim, Stony Brook University

# required packages
require(maps)
require(geosphere)
require(raster)

--------------------------------------------------------------------------------------------
# To Do:
# Clean up code
--------------------------------------------------------------------------------------------

# set working directory

wd = "C:/Users/mcwlim/Dropbox/NASA_Anusha/MarisaRAstuff/BBLrecapproj"
setwd(wd)

################# Updated BBL 11April2014 #################
BBLdat <- read.csv("11April14_BBLdata.csv", sep=",", header=T)
head(BBLdat)
names(BBLdat)
dim(BBLdat)

# First, let's clean up the data. Information that needs to be cleaned: (quick filter check of excel file shows that these need cleaning)
# 1. incorrect months of year (in ENCOUNTER_MONTH)
# 2. incorrect days of month (in ENCOUNTER_DAY)
# 3. unidentified species (in B_SPECIES_NAME)
# 4. missing lat/long information (in E_LAT_DECIMAL_DEGREES, E_LON_DECIMAL_DEGREES)

BBLdat1 <- subset(BBLdat, BBLdat$ENCOUNTER_MONTH <= 12)
BBLdat2 <- subset (BBLdat1, BBLdat1$ENCOUNTER_DAY <= 31)
BBLdat3 <- subset(BBLdat2, BBLdat2$B_SPECIES_NAME != "Unidentified Hummingbird")
BBLdat4 <- subset(BBLdat3, BBLdat3$E_LAT_DECIMAL_DEGREES != "NA")
dim(BBLdat4)

# based on 11April14 data, there were 28 entries removed, 1293 records used 

# set extent
xlim <- c(-171.738281, -56.601563)
ylim <- c(8, 71.856229)
# set color for connecting lines
colors <- "#FF3300"

# maps of each species, banding and recapture points connected with a line
species <- unique(BBLdat4$B_SPECIES_NAME)
for(i in 1:length(species)){
  jpeg(paste("species_", species[i], ".jpg"), height=5, width=5, units="in", res=500)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- BBLdat4[BBLdat4$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  bandnums <- unique(birdies$BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.35)
  title(sub=paste("# unique birds = ", length(bandnums)), cex.sub=1.5, line=0.5)
  for(j in 1:nrow(birdies)){
    inter <- gcIntermediate(c(birdies[j,]$B_LON_DECIMAL_DEGREES, birdies[j,]$B_LAT_DECIMAL_DEGREES),
                            c(birdies[j,]$E_LON_DECIMAL_DEGREES, birdies[j,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col=colors, lwd=2)
    points(birdies$B_LAT_DECIMAL_DEGREES ~ birdies$B_LON_DECIMAL_DEGREES, cex=1, pch=20, col="black")  
  }
  dev.off()
}


################# Test data with time, connections, sex, and age #################

humtest <- read.csv("humtest27aug.csv", sep=",", header=T)
humtest <- humtest[order(humtest$B_SPECIES_NAME, humtest$B_BAND_NUM, humtest$BANDING_YEAR, humtest$BANDING_MONTH, humtest$BANDING_DATE),]
head(humtest)
names(humtest)
dim(humtest)

# set extent
xlim <- c(-171.738281, -56.601563)
ylim <- c(8, 71.856229)
# set color for connecting lines
colors <- "#FF3300"

# maps of each species, banding and recapture points connected with a line
species <- unique(humtest$B_SPECIES_NAME)
for(i in 1:length(species)){
  pdf(paste("species_", species[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  bandnums <- unique(birdies$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  for(j in 1:nrow(birdies)){
    inter <- gcIntermediate(c(birdies[j,]$B_LON_DECIMAL_DEGREES, birdies[j,]$B_LAT_DECIMAL_DEGREES),
                            c(birdies[j,]$E_LON_DECIMAL_DEGREES, birdies[j,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col=colors, lwd=2)
    points(birdies$B_LAT_DECIMAL_DEGREES ~ birdies$B_LON_DECIMAL_DEGREES, cex=1, pch=20, col="black")  }
  dev.off()
}

# maps of each species, points colored by sex (classes are defined as female or male) 
species <- unique(humtest$B_SPECIES_NAME)
for(i in 1:length(species)){
  pdf(paste("Sp_", species[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ]
  females <- birdies[birdies$B_SEX_CODE == "5" | birdies$B_SEX_CODE == "7", ]
  males <- birdies[birdies$B_SEX_CODE == "4" | birdies$B_SEX_CODE == "6", ]
  points(females$B_LAT_DECIMAL_DEGREES ~ females$B_LON_DECIMAL_DEGREES, cex=1, pch=1, col="red")
  points(males$B_LAT_DECIMAL_DEGREES ~ males$B_LON_DECIMAL_DEGREES, cex=1, pch=1, col="blue")
  legend("left", ncol=2, y.intersp=1, x.intersp=0.5, xjust=0, bg="whitesmoke", pch=20, cex=1, col=c("red", "blue"), c("Female", "Male"))
  title(main=species[i])
  dev.off()
}

# maps connecting recaptures, with 1st banding and recapture in same row and,
# unique individuals (male or female) are marked as points, they are shown for the original banding location only
species <- unique(humtest$B_SPECIES_NAME)
for(i in 1:length(species)){
  pdf(paste("species_", species[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  bandnums <- unique(birdies$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  for(j in 1:nrow(birdies)){
    inter <- gcIntermediate(c(birdies[j,]$B_LON_DECIMAL_DEGREES, birdies[j,]$B_LAT_DECIMAL_DEGREES),
                            c(birdies[j,]$E_LON_DECIMAL_DEGREES, birdies[j,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col=colors, lwd=2)
    females <- birdies[birdies$B_SEX_CODE == "5" | birdies$B_SEX_CODE == "7", ]
    points(females$B_LAT_DECIMAL_DEGREES ~ females$B_LON_DECIMAL_DEGREES, cex=1, pch=1, col="red")
    males <- birdies[birdies$B_SEX_CODE == "4" | birdies$B_SEX_CODE == "6", ]
    points(males$B_LAT_DECIMAL_DEGREES ~ males$B_LON_DECIMAL_DEGREES, cex=1, pch=1, col="blue")
  }
  dev.off()
}

# mapping males and females with the lines, separately
# Males only.
species <- unique(humtest$B_SPECIES_NAME)
for(i in 1:length(species)){
  pdf(paste("species_", species[i], "_m", ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  males <- birdies[birdies$B_SEX_CODE == "4" | birdies$B_SEX_CODE == "6", ] #all rows for a specificed sex
  bandnums <- unique(males$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  points(males$B_LAT_DECIMAL_DEGREES ~ males$B_LON_DECIMAL_DEGREES, cex=1, pch=1, col="blue")
  for(j in 1:nrow(males)){
    if(nrow(males) != 0){
      inter <- gcIntermediate(c(males[j,]$B_LON_DECIMAL_DEGREES, males[j,]$B_LAT_DECIMAL_DEGREES),
                              c(males[j,]$E_LON_DECIMAL_DEGREES, males[j,]$E_LAT_DECIMAL_DEGREES),
                              n=50, addStartEnd=TRUE,)
      lines(inter, col=colors, lwd=0.5) 
    }
  }
  dev.off()
}
# Females only. 
# note: had to add the if statement because some species have no female records
species <- unique(humtest$B_SPECIES_NAME)
for(i in 1:length(species)){
  pdf(paste("species_", species[i], "_f", ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  females <- birdies[birdies$B_SEX_CODE == "5" | birdies$B_SEX_CODE == "7", ] #all rows for a specificed sex
  bandnums <- unique(females$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  points(females$B_LAT_DECIMAL_DEGREES ~ females$B_LON_DECIMAL_DEGREES, cex=1, pch=1, col="red")
  for(j in 1:nrow(females)){
    if(nrow(females) != 0){
      inter <- gcIntermediate(c(females[j,]$B_LON_DECIMAL_DEGREES, females[j,]$B_LAT_DECIMAL_DEGREES),
                              c(females[j,]$E_LON_DECIMAL_DEGREES, females[j,]$E_LAT_DECIMAL_DEGREES),
                              n=50, addStartEnd=TRUE,)
      lines(inter, col=colors, lwd=0.5)
    }
  }
  dev.off()
}

# age data
# Males and age
species <- unique(humtest$B_SPECIES_NAME)

for(i in 1:length(species)){
  pdf(paste("species_", species[i], "_m", ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  males <- birdies[birdies$B_SEX_CODE == "4" | birdies$B_SEX_CODE == "6", ] #all rows for a specificed sex
  males <- males[order(males$B_AGE_CODE), ]
  bandnums <- unique(males$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  for(j in 1:nrow(males)){
    if(nrow(males) != 0){
      inter <- gcIntermediate(c(males[j,]$B_LON_DECIMAL_DEGREES, males[j,]$B_LAT_DECIMAL_DEGREES),
                              c(males[j,]$E_LON_DECIMAL_DEGREES, males[j,]$E_LAT_DECIMAL_DEGREES),
                              n=50, addStartEnd=TRUE,)
      lines(inter, col=colors, lwd=0.5) 
      ages <- unique(males$B_AGE_CODE)
      points(males$B_LAT_DECIMAL_DEGREES ~ males$B_LON_DECIMAL_DEGREES, cex=1, pch=ages, col="blue")
      legend("left", col="blue", pch=ages, legend=ages, bg="whitesmoke")
    }
  }
  dev.off()
}
# Females & age
# note: had to add the if statement because some species have no female records
species <- unique(humtest$B_SPECIES_NAME)
for(i in 1:length(species)){
  pdf(paste("species_", species[i], "_f", ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=T, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  females <- birdies[birdies$B_SEX_CODE == "5" | birdies$B_SEX_CODE == "7", ] #all rows for a specificed sex
  females <- females[order(females$B_AGE_CODE), ]
  bandnums <- unique(females$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  for(j in 1:nrow(females)){
    if(nrow(females) != 0){
      inter <- gcIntermediate(c(females[j,]$B_LON_DECIMAL_DEGREES, females[j,]$B_LAT_DECIMAL_DEGREES),
                              c(females[j,]$E_LON_DECIMAL_DEGREES, females[j,]$E_LAT_DECIMAL_DEGREES),
                              n=50, addStartEnd=TRUE,)
      lines(inter, col=colors, lwd=0.5)
      ages <- unique(females$B_AGE_CODE)
      points(females$B_LAT_DECIMAL_DEGREES ~ females$B_LON_DECIMAL_DEGREES, cex=1, pch=ages, col="red")
      legend("left", col="red", pch=ages, legend=ages, bg="whitesmoke")
    }
  }
  dev.off()
}

# For RUHU: showing how points change within year, can track the migration pattern, also showing ages 
RUHU <- humtest[humtest$B_SPECIES_NAME == "RUFOUS HUMMINGBIRD", ]
RUHUorder <- RUHU[order(RUHU$BANDING_MONTH), ]
RUHUmonth <- unique(RUHUorder$BANDING_MONTH)
for(i in 1:length(RUHUmonth)){
  pdf(paste("Month_", RUHUmonth[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  title(main=RUHUmonth[i])
  mopts <- RUHU[RUHU$BANDING_MONTH == RUHUmonth[i], ]
  mopts <- mopts[order(mopts$B_AGE_CODE), ]
  ages <- unique(mopts$B_AGE_CODE)
  points(mopts$B_LAT_DECIMAL_DEGREES ~ mopts$B_LON_DECIMAL_DEGREES, col="black", cex=1, pch=ages)
  legend("left", col="black", pch=ages, legend=ages, bg="whitesmoke")
  dev.off()
}

# For RUHU: same as above, but adding points to show male vs. female - is there a difference in timing of migration between the sexes?
RUHU <- humtest[humtest$B_SPECIES_NAME == "RUFOUS HUMMINGBIRD", ]
RUHUorder <- RUHU[order(RUHU$BANDING_MONTH), ]
RUHUmonth <- unique(RUHUorder$BANDING_MONTH)
for(i in 1:length(RUHUmonth)){
  pdf(paste("Month_", RUHUmonth[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  title(main=RUHUmonth[i])
  mopts <- RUHU[RUHU$BANDING_MONTH == RUHUmonth[i], ]
  mopts <- mopts[order(mopts$B_AGE_CODE), ]
  ages <- unique(mopts$B_AGE_CODE)
  females <- mopts[mopts$B_SEX_CODE == "5" | mopts$B_SEX_CODE == "7", ]
  males <- mopts[mopts$B_SEX_CODE == "4" | mopts$B_SEX_CODE == "6", ]
  points(mopts$B_LAT_DECIMAL_DEGREES ~ mopts$B_LON_DECIMAL_DEGREES, col="black", cex=1, pch=ages)
  points(females$B_LAT_DECIMAL_DEGREES ~ females$B_LON_DECIMAL_DEGREES, cex=1, pch=20, col="red")
  points(males$B_LAT_DECIMAL_DEGREES ~ males$B_LON_DECIMAL_DEGREES, cex=1, pch=20, col="blue")
  legend("left", col="black", pch=ages, legend=ages, bg="whitesmoke")
  dev.off()
}

# maps for each species by month (pts are colored by sex and shaped by age)
species <- unique(humtest$B_SPECIES_NAME)
hummonth <- unique(humtest$BANDING_MONTH)
for(i in 1:length(species)){
  birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ]
  for(j in 1:length(hummonth)){
    pdf(paste("Sp_", species[i], "_Month_", hummonth[j], ".pdf", sep=""), width=11, height=7)
    map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
    title(main=hummonth[j])
    mopts <- birdies[birdies$BANDING_MONTH == hummonth[j], ]
    mopts <- mopts[order(mopts$B_AGE_CODE), ]
    females <- mopts[mopts$B_SEX_CODE == "5" | mopts$B_SEX_CODE == "7", ]
    males <- mopts[mopts$B_SEX_CODE == "4" | mopts$B_SEX_CODE == "6", ]
    points(females$B_LAT_DECIMAL_DEGREES ~ females$B_LON_DECIMAL_DEGREES, cex=1, pch=20, col="red")
    points(males$B_LAT_DECIMAL_DEGREES ~ males$B_LON_DECIMAL_DEGREES, cex=1, pch=20, col="blue")
    ages <- unique(mopts$B_AGE_CODE)
    if(length(ages) != 0){
      points(mopts$B_LAT_DECIMAL_DEGREES ~ mopts$B_LON_DECIMAL_DEGREES, col="black", cex=1, pch=ages)
      colors <- c(rep("black", length(ages)), "red", "blue")
      legend("left", col=colors, pch=c(ages, 20, 20), legend=c(ages, "Female", "Male"), bg="whitesmoke") 
      }
    dev.off()
  }
}

# by year for RUHU, RTHU, BCHU - where banded and recap'ed in same year
# pattern with el nino years? are hummers doing something different in each year? 
colors <- rainbow(12, start=0.4, end=1, alpha=0.8)
BCHUpts <- humtest[humtest$B_SPECIES_NAME == "BLACK-CHINNED HUMMINGBIRD", ]
BCHUpts <- BCHUpts[order(BCHUpts$BANDING_YEAR), ]
# new df of all birds that were recaptured in the same year that they were banded
BCHUpts <- BCHUpts[BCHUpts$BANDING_YEAR == BCHUpts$ENCOUNTER_YEAR, ] 
BCHUyr <- unique(BCHUpts$BANDING_YEAR)
for(i in 1:length(BCHUyr)){
  yrpts <- BCHUpts[BCHUpts$BANDING_YEAR == BCHUyr[i], ]
  pdf(paste("BCHU_", BCHUyr[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  title(main=BCHUyr[i])
  for(j in 1:nrow(yrpts)){
    inter <- gcIntermediate(c(yrpts[j,]$B_LON_DECIMAL_DEGREES, yrpts[j,]$B_LAT_DECIMAL_DEGREES),
                          c(yrpts[j,]$E_LON_DECIMAL_DEGREES, yrpts[j,]$E_LAT_DECIMAL_DEGREES),
                          n=50, addStartEnd=TRUE,)
  lines(inter, col="#FF3300", lwd=2)
  points(yrpts$B_LAT_DECIMAL_DEGREES ~ yrpts$B_LON_DECIMAL_DEGREES, col="black", bg=colors, cex=1.5, pch=21)
  legend("left", ncol=4, y.intersp=1, x.intersp=0.5, xjust=0, bg="whitesmoke", pch=21, cex=1, col="black", pt.bg=colors, c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
  }
  dev.off()
}
RUHUpts <- humtest[humtest$B_SPECIES_NAME == "RUFOUS HUMMINGBIRD", ]
RUHUpts <- RUHUpts[order(RUHUpts$BANDING_YEAR), ]
RUHUpts <- RUHUpts[RUHUpts$BANDING_YEAR == RUHUpts$ENCOUNTER_YEAR, ] 
RUHUyr <- unique(RUHUpts$BANDING_YEAR) 
for(i in 1:length(RUHUyr)){
  yrpts <- RUHUpts[RUHUpts$BANDING_YEAR == RUHUyr[i], ]
  pdf(paste("RUHU_", RUHUyr[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  title(main=RUHUyr[i])
  for(j in 1:nrow(yrpts)){
    inter <- gcIntermediate(c(yrpts[j,]$B_LON_DECIMAL_DEGREES, yrpts[j,]$B_LAT_DECIMAL_DEGREES),
                            c(yrpts[j,]$E_LON_DECIMAL_DEGREES, yrpts[j,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col="#FF3300", lwd=2)
    points(yrpts$B_LAT_DECIMAL_DEGREES ~ yrpts$B_LON_DECIMAL_DEGREES, col="black", bg=colors, cex=1.5, pch=21)
    legend("left", ncol=4, y.intersp=1, x.intersp=0.5, xjust=0, bg="whitesmoke", pch=21, cex=1, col="black", pt.bg=colors, c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
  }
  dev.off()
}  
RTHUpts <- humtest[humtest$B_SPECIES_NAME == "RUBY-THROATED HUMMINGBIRD", ]
RTHUpts <- RTHUpts[order(RTHUpts$BANDING_YEAR), ]
RTHUpts <- RTHUpts[RTHUpts$BANDING_YEAR == RTHUpts$ENCOUNTER_YEAR, ] 
RTHUyr <- unique(RTHUpts$BANDING_YEAR)
for(i in 1:length(RTHUyr)){
  yrpts <- RTHUpts[RTHUpts$BANDING_YEAR == RTHUyr[i], ]
  pdf(paste("RTHU_", RTHUyr[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  title(main=RTHUyr[i])
  for(j in 1:nrow(yrpts)){
    inter <- gcIntermediate(c(yrpts[j,]$B_LON_DECIMAL_DEGREES, yrpts[j,]$B_LAT_DECIMAL_DEGREES),
                            c(yrpts[j,]$E_LON_DECIMAL_DEGREES, yrpts[j,]$E_LAT_DECIMAL_DEGREES),
                            n=50, addStartEnd=TRUE,)
    lines(inter, col="#FF3300", lwd=2)
    points(yrpts$B_LAT_DECIMAL_DEGREES ~ yrpts$B_LON_DECIMAL_DEGREES, col="black", bg=colors, cex=1.5, pch=21)
    legend("left", ncol=4, y.intersp=1, x.intersp=0.5, xjust=0, bg="whitesmoke", pch=21, cex=1, col="black", pt.bg=colors, c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
  }
  dev.off()
}

################# BBL basic metrics of the data #################
# counting number of individuals, females, males, and juveniles per species
numbers <- function(){
  indiv.v <- c()
  species <- unique(humtest$B_SPECIES_NAME)
  for(i in 1:length(species)){
    birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ]
    indivs <- length(unique(birdies$B_BAND_NUM))
    indiv.v[i] <- indivs
  }
  fem.v <- c()
  for(i in 1:length(species)){
    birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ]
    females <- birdies[birdies$B_SEX_CODE == "5" | birdies$B_SEX_CODE == "7", ]
    females <- length(unique(females$B_BAND_NUM))
    fem.v[i] <- females
  }
  male.v <- c()
  for(i in 1:length(species)){
    birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ]
    males <- birdies[birdies$B_SEX_CODE == "4" | birdies$B_SEX_CODE == "6", ]
    males <- length(unique(males$B_BAND_NUM))
    male.v[i] <- males
  }
  #assuming juveniles are AHY, HY, or J
  juv.v <- c()
  for(i in 1:length(species)){
    birdies <- humtest[humtest$B_SPECIES_NAME == species[i], ]
    juvs <- birdies[birdies$B_AGE_CODE == "1" | birdies$B_AGE_CODE == "2" | birdies$B_AGE_CODE == "3", ]
    juvs <- length(unique(juvs$B_BAND_NUM))
    juv.v[i] <- juvs
  }
  allnums <- data.frame(species, indiv.v, fem.v, male.v, juv.v)
  write.table(allnums, file="24sep13_BBLdata.csv", sep=",")
}
numbers()

# calculating the distance between points in km, Euclidian distance
A <- SpatialPointsDataFrame(coords=cbind(humtest$B_LON_DECIMAL_DEGREES,humtest$B_LAT_DECIMAL_DEGREES),humtest)
B <- SpatialPointsDataFrame(coords=cbind(humtest$E_LON_DECIMAL_DEGREES,humtest$E_LAT_DECIMAL_DEGREES),humtest)
distances <- pointDistance(A, B, longlat=TRUE)/1000 # /1000 to get km instead of meters
hist(distances, col="plum4")
newhumtest <- data.frame(A, distances)
newhumtest2 <- data.frame(humtest, distances)

distcalcs <- function(){
  species <- unique(newhumtest2$B_SPECIES_NAME)
  meandist.v <- c()
  sddist.v <- c()
  meddist.v <- c()
  rangedist.m <- matrix(nrow=14, ncol=2)
  for(i in 1:length(species)){
    birdies <- newhumtest2[newhumtest2$B_SPECIES_NAME == species[i], ]
    #for each species, i want to calculate the mean distance, sd, median + include the range of values
    meandist <- mean(birdies$distances)
    sddist <- sd(birdies$distances)
    meddist <- median(birdies$distances)
    rangedist <- range(birdies$distances)
    meandist.v[i] <- meandist
    sddist.v[i] <- sddist
    meddist.v[i] <- meddist  
    rangedist.m[i, ] <- rangedist
  }
  allcalcs <- data.frame(meandist.v, sddist.v, meddist.v, rangedist.m)
  write.table(allcalcs, file="25sep13_distcalcs.csv", sep=",")
}

distcalcs()

################# SDM test code #################
# not working yet
hum.sp<-SpatialPointsDataFrame(coords=cbind(humtest$B_LON_DECIMAL_DEGREES,humtest$B_LAT_DECIMAL_DEGREES),humtest)

require(dismo)
USA<-getData('GADM', country='USA', level=1)

plot(USA)

require(raster)
u<-raster()
r<-rasterize(USA,u)
plot(r)
r[r>0]<-1
plot(r)
plot(r,ext=extent(hum.sp))
points(hum.sp)

predictors <- stack(list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
                               pattern='grd', full.names=TRUE ))
me<-maxent(predictors,hum.sp)

################# OLD/test code for recap connections by species #################
# need to clean up/delete? 

library(maps)
library(geosphere)
library(adegenet)

xlim <- c(-171.738281, -56.601563)
ylim <- c(8, 71.856229)

allhum <- read.csv("testhummerdata.csv", sep=",", header=TRUE)
names(allhum)
head(allhum)
allhum2 <- data.frame(allhum[,1:9]) #allhum has 6 empty columns for some reason, got rid of here
head(allhum2)
#putting data in order by species, bandnumber, year, month, and then day
allhum2 <- allhum2[order(allhum2$B_SPECIES_NAME, allhum2$B_BAND_NUM, allhum2$BANDING_YEAR, allhum2$BANDING_MONTH, allhum2$BANDING_DATE), ]
head(allhum2)
species <- unique(allhum2$B_SPECIES_NAME)
colors <- "#FF3300"
# pal <- colorRampPalette(c("orange", "green"))
# colors <- pal(50)

# ##test plotting of datapoints by species
# plot.new()
# map("world", xlim=xlim, ylim=ylim, col="whitesmoke", fill=T, bg="azure3", lwd=0.05)
# birdies <- allhum2[allhum2$B_SPECIES_NAME == species[11], ]
# points(birdies$B_LAT_DECIMAL_DEGREES ~ birdies$B_LON_DECIMAL_DEGREES, cex=2, col="red", pch=20)


for(i in 1:length(species)){
  pdf(paste("species_", species[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  birdies <- allhum2[allhum2$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
  bandnums <- unique(birdies$B_BAND_NUM)
  par(adj=0.5)
  title(main=species[i], col.main="black")
  par(adj=0.55)
  title(sub=length(bandnums), cex.sub=1.5, line=0.5)
  par(adj=0.35)
  title(sub="# unique birds =", line=0.5, cex.sub=1.5)
  for(j in 1:length(bandnums)){
    bird <- subset(birdies, B_BAND_NUM == bandnums[j])
    #     bird <- bird[order(bird$BANDING_DATE), ] #how does R read dates?
    for(k in 1:nrow(bird)){
      if((k+1) <= nrow(bird)){
        temp <- bird[k:(k+1), ]
        inter <- gcIntermediate(c(temp[1,]$B_LON_DECIMAL_DEGREES, temp[1,]$B_LAT_DECIMAL_DEGREES), 
                                c(temp[2,]$B_LON_DECIMAL_DEGREES, temp[2,]$B_LAT_DECIMAL_DEGREES), 
                                n=50, addStartEnd=TRUE, )
        lines(inter, col=colors, lwd=2)
      }
    }  
  }
  dev.off()
}

#new code gets rid of this: Error in if (antipodal(p1, p2)) { : missing value where TRUE/FALSE needed

#code below goes through individual species at a time for recap maps
ALHU <- allhum2[allhum2$B_SPECIES_NAME == "ALLEN'S HUMMINGBIRD", ] #df with just data for Allen's hbird
ANHU <- allhum2[allhum2$B_SPECIES_NAME == "ANNA'S HUMMINGBIRD", ]
BCHU <- allhum2[allhum2$B_SPECIES_NAME == "BLACK-CHINNED HUMMINGBIRD", ]
BLUH <- allhum2[allhum2$B_SPECIES_NAME == "BLUE-THROATED HUMMINGBIRD", ]
BBLH <- allhum2[allhum2$B_SPECIES_NAME == "BROAD-BILLED HUMMINGBIRD", ]
BTLH <- allhum2[allhum2$B_SPECIES_NAME == "BROAD-TAILED HUMMINGBIRD", ] #problematic
BUFH <- allhum2[allhum2$B_SPECIES_NAME == "BUFF-BELLIED HUMMINGBIRD", ]
CAHU <- allhum2[allhum2$B_SPECIES_NAME == "CALLIOPE HUMMINGBIRD", ]
COHU <- allhum2[allhum2$B_SPECIES_NAME == "COSTA'S HUMMINGBIRD", ]
MAHU <- allhum2[allhum2$B_SPECIES_NAME == "MAGNIFICENT HUMMINGBIRD", ]
RTHU <- allhum2[allhum2$B_SPECIES_NAME == "RUBY-THROATED HUMMINGBIRD", ] #problematic
RUHU <- allhum2[allhum2$B_SPECIES_NAME == "RUFOUS HUMMINGBIRD", ]
VCHU <- allhum2[allhum2$B_SPECIES_NAME == "VIOLET-CROWNED HUMMINGBIRD", ]
WEHU <- allhum2[allhum2$B_SPECIES_NAME == "WHITE-EARED HUMMINGBIRD", ]

# prints out the number of samples/species
for(i in 1:length(species)){
  birdies <- allhum2[allhum2$B_SPECIES_NAME == species[i], ]
  bands <- unique(birdies$B_BAND_NUM)
  print(length(bands))
}

#returns rows with same band number -> so it shows the original and recapture
# bird1 <- allhum2[allhum2$B_BAND_NUM == allens[20, ]$B_BAND_NUM, ] 

#tutorial testdrive data
# plot.new()
# map("world", col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05, xlim=xlim, ylim=ylim)
# 
# allenslato <- 37.91667
# allenslongo <- -122.75
# allenslatr <- 29.00
# allenslongr <- -120.00
# inter <- gcIntermediate(c(allenslongo, allenslato), c(allenslongr, allenslatr), n=50, addStartEnd=TRUE)
# lines(inter, col="green")
# lat_tx <- 29.954935
# lon_tx <- -98.701172
# inter2 <- gcIntermediate(c(allenslongo, allenslato), c(lon_tx, lat_tx), n=50, addStartEnd=TRUE)
# lines(inter2, col="red")
# inter3 <- gcIntermediate(c(allenslongr, allenslatr), c(lon_tx, lat_tx), n=50, addStartEnd=TRUE)
# lines(inter3, col="blue")

############################### Old temporal data by species code
# all points
plot.new()
map("world", col="#f2f2f2", fill=TRUE, bg="steelblue", lwd=0.05, xlim=xlim, ylim=ylim)
points(allhum2$B_LAT_DECIMAL_DEGREES ~ allhum2$B_LON_DECIMAL_DEGREES, pch=20, cex=0.5)

# pal <- colorRampPalette(c("orange", "blue"))
# colors <- pal(14)
# plot.new()
# map("world", col="#f2f2f2", fill=TRUE, bg="steelblue", lwd=0.05, xlim=xlim, ylim=ylim)
# points(RUHU$B_LAT_DECIMAL_DEGREES ~ RUHU$B_LON_DECIMAL_DEGREES, pch=20, cex=1, col=transp(colors))

#maps out points by species, colors by month
# birdyears <- unique(birdies$BANDING_YEAR)
# for(i in 1:length(species)){
#   pdf(paste("species_", species[i], ".pdf", sep=""), width=11, height=7)
#   map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
#   birdies <- allhum2[allhum2$B_SPECIES_NAME == species[i], ] #all rows for a specified species i
#   par(adj=0.5)
#   title(main=species[i], col.main="black")
# #   legend("bottomleft", ncol=2, y.intersp=0.8, x.intersp=0.4, xjust=0, bg="whitesmoke", col=birdies$BANDING_MONTH, pch=20, cex=2, c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
#   for(j in 1:length(birdyears)){
#     bird <- subset(birdies, BANDING_YEAR == birdyears[j])
#     points(birdies$B_LAT_DECIMAL_DEGREES ~ birdies$B_LON_DECIMAL_DEGREES, col=birdies$BANDING_MONTH, cex=1.5, pch=21)
#   } 
#   dev.off()
# }

# pal <- colorRampPalette(c("orange", "blue"))
# colors <- pal(12)

#loop goes through each species and plots all points for unique years of data of that species
plot.new()
species <- unique(allhum2$B_SPECIES_NAME)
colors <- rainbow(12, start=0.4, end=1, alpha=0.8)
for(i in 1:length(species)){
  birdies <- allhum2[allhum2$B_SPECIES_NAME == species[i], ]
  birdyears <- unique(birdies$BANDING_YEAR)
  for(j in 1:length(birdyears)){
    yrpts <- birdies[birdies$BANDING_YEAR == birdyears[j], ]
    pdf(paste("Sp_", species[i], "_Year_", birdyears[j], ".pdf", sep=""), width=11, height=7)
    map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
    points(yrpts$B_LAT_DECIMAL_DEGREES ~ yrpts$B_LON_DECIMAL_DEGREES, cex=1.5, pch=20, col=colors)
    title(main=birdyears[j])
    legend("left", ncol=4, y.intersp=1, x.intersp=0.5, xjust=0, bg="whitesmoke", pch=21, cex=1, col="black", pt.bg=colors, c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
    dev.off()
  }
}

#for one species at a time
plot.new()
RUHUorder <- RUHU[order(RUHU$BANDING_YEAR), ]
RUHUyears <- unique(RUHUorder$BANDING_YEAR)
for(i in 1:length(RUHUyears)){
  yrpts <- RUHU[RUHU$BANDING_YEAR == RUHUyears[i], ]
  pdf(paste("Year_", RUHUyears[i], ".pdf", sep=""), width=11, height=7)
  map("world", col="whitesmoke", fill=TRUE, bg="azure3", lwd=0.05, xlim=xlim, ylim=ylim)
  points(yrpts$B_LAT_DECIMAL_DEGREES ~ yrpts$B_LON_DECIMAL_DEGREES, col=colors, cex=1.5, pch=20)
  title(main=RUHUyears[i])
  legend("left", ncol=4, y.intersp=1, x.intersp=0.5, xjust=0, bg="whitesmoke", pch=21, cex=1, col="black", pt.bg=colors, c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
  dev.off()
}

#manual way of plotting different years, have to change ALHUyears number
plot.new()
map("world", xlim=xlim, ylim=ylim)
yrpts <- ALHU[ALHU$BANDING_YEAR == ALHUyears[4], ]
points(yrpts$B_LAT_DECIMAL_DEGREES ~ yrpts$B_LON_DECIMAL_DEGREES, col="pink", cex=1.5, pch=20)
