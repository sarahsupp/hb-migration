# runs through all hb species to analyze impacts of environmental data on migration route
# Sarah Supp and Tina Cormier

# import libraries
library(randomForest)
library(sm)
library(ggplot2)
library(ltm)
library(Rarity)
library(lme4)

# define pathnames
function.dir <- "C:/Users/sarah/Documents/github/hb-migration/hb_RS_functions.R"
agan.dir <- "C:/Users/sarah/Dropbox/ebird_annotated_raw/"
migtime.dir <- "C:/Users/sarah/Dropbox/ebird_annotated_raw/"
fig.dir <- "C:/Users/sarah/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/"

# species codes
spcodes <- c("rthu", "bchu", "bthu", "cahu", "ruhu")

for (spp in unique(spcodes)){

  #read in annotaed data
  ann <- read.csv(paste0(agan.dir, spp, "/", spp, "_lag0_allYears.csv"), as.is=T)
  
  # set column types and add new date columns
  ann$timestamp <- as.Date(ann$timestamp, format='%Y-%m-%d %H:%M:%S.000')
  ann$presence <- as.factor(ann$presence)
  ann$GlobCover <- as.factor(ann$GlobCover)
  ann$month <- factor(format(ann$timestamp, "%m"), ordered=T)
  ann$year <- factor(format(ann$timestamp, "%Y"), ordered=T)
  ann$doy <- as.numeric(format(ann$timestamp, "%j"))
  
  #keep only the columns we need for analysis
  ann <- ann[,names(ann) %in% c()]
  
  #Separate by season using the migration timing table from the previous hb study
  migtime <- read.table(paste0(migtime.dir, spp, "/", "west_migration_", spp, ".txt"), header=T, sep=" ")
  years <- unique(ann$year)
  
  #make a dataframe to hold data during spring and fall migrations
  # Note: data that falls outside of seasons is not included in the new dataframe (sf)
  sf <- as.data.frame(matrix(data=NA, nrow=0, ncol=ncol(ann)+1))
  names(sf) <- c(names(ann), "season")
  
  for (yr in years) {
    seasonal <- seasonalSub(ann, migtime, yr)
    sf <- rbind(sf, seasonal)
    }

  #convert K to celsius
  sf$t10m <- sf$t10m - 273.15
  sf$Temp_sfc <- sf$Temp_sfc - 
  
  sf$occurrence[sf$presence==1] <- "present"
  sf$occurrence[sf$presence==0] <- "absent"
  

#---------------------------------plot doy representation to check that it is similar for pres and abs
  pres <- ann[ann$presence == 1,]
  abs <- ann[ann$presence == 0,]

  ggplot(data=pres, aes(doy)) + geom_histogram(fill="red", alpha=0.5) + 
    geom_histogram(data=abs, aes(doy), fill="blue", alpha=0.5) + 
    theme_classic() + ggtitle("blue = absent, red = present")
  ggsave(filename=paste0(fig.dir, spp, "_doy_hist.png"), width=4, height=3, res=400)


#----------------------------------pairwise correlation plot
  keepcols <- c("Temp_sfc", "Total_precipitation_sfc", "EVI", "swrf", "lwrf", "u10m", "uplift", "SRTM_elev" )
  df.cor <- sf[,names(sf) %in% keepcols]
  df.cor <- df.cor[complete.cases(df.cor),]

  out <- paste0(fig.dir, spp, "_variable_correlations.png")
  png(out, width=12, height=12, units="in", res=600)
  pairs(df.cor, diag.panel=panel.hist, upper.panel=panel.cor)
  dev.off()


#----------------------------Random forest - for variable importance
  preds <- preds[,-1]
  response <- preds$presence

  #keep just the environmental vars here
  preds.rf <- preds[,c(2,3,7,8,9,10,11,12,13,15,16,18,21,23,24,35,36,37)]
  rf <- randomForest(preds.rf, response, na.rm=T, importance=T)

  # Use variable importance plot based on mean decrease in GINI. Calle & Urrea (2010) show its rankings are more stable.
  #   Explanation here: http://www.stat.berkeley.edu/~breiman/RandomForests/cc%5Fhome.htm#varimp
  varImp <- importance(rf,type = 2)

  # Plot results
  rf.graph <- paste0(fig.dir, spp, "/", spp, "_varImp_allvars.pdf")
  pdf(file = rf.graph, width=8, height=8)
  par(mar=c(8,6,6,4))
  varImpPlot(rf, type=2,pch=16, col="blue",main=paste(spp, "\nRandom Forest Variable Importance", sep=' '))
  dev.off()


#--------------------------------- Generalized linear mixed model with year and season as random effects
  # use only rows without NA (can't standardize NA values below)
  sf2 <- sf[complete.cases(sf),]

  #standardize the data
  zscore = apply(sf2[,names(sf2)%in% vars ], 2, function(x) {
    y = (x - mean(x))/sd(x)
    return(y)
  })

  #make into data frame and add ID and time vars back in
  zscore <- as.data.frame(zscore)
  zscore$presence <- sf2$presence
  zscore$year <- as.factor(sf2$year)
  zscore$season <- sf2$season

  fit <- glmer(presence ~ Temp_sfc + t10m + swrf + lwrf + EVI + Total_precipitation_sfc + uplift + SRTM_elev + 
               (1|year) + (1|season), data = sf, family = "binomial", 
               control = glmerControl(optimizer = "bobyqa"))
  summary(fit)
}

