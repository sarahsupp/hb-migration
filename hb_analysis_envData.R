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
function.dir <- "/Users/sarah/Documents/github/hb-migration/hb_RS_functions.R"
agan.dir <- "/Users/sarah/Desktop/Dropbox/Hummingbirds/hb_migration_data/ebird_annotated_fil/"
migtime.dir <- "/Users/sarah/Desktop/Dropbox/Hummingbirds/hb_migration_data/ebird_annotated_fil/"
fig.dir <- "/Users/sarah/Desktop/Dropbox/Hummingbirds/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/"

# source function script
source(function.dir)

# species codes
spcodes <- c("rthu", "bchu", "bthu", "cahu", "ruhu")

for (spp in unique(spcodes)){

  #read in annotaed data
  ann <- read.csv(paste0(agan.dir, spp, "/", spp, "_lag0_allYears_fil.csv"), as.is=T)
  
  # set column types and add new date columns
  ann$timestamp <- as.Date(ann$timestamp, format='%Y-%m-%d %H:%M:%S.000')
  ann$presence <- as.factor(ann$presence)
  ann$GlobCover <- as.factor(ann$GlobCover)
  ann$month <- factor(format(ann$timestamp, "%m"), ordered=T)
  ann$year <- factor(format(ann$timestamp, "%Y"), ordered=T)
  ann$doy <- as.numeric(format(ann$timestamp, "%j"))
  
  #keep only the columns we need for analysis
  ann <- ann[,names(ann) %in% c("location.long", "location.lat", "presence", "lwrf", "swrf", "t10m", "u10m", 
                                "lwrf_up", "swrf_up", "EVI", "SRTM_elev", "Total_precipitation_sfc", 
                                "month", "year", "doy")]
  
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
  
  #TODO: double check main vars to use, update here and elsewhere throughout code
#  vars <- c("t10m", "EVI", "swrf", "lwrf", "swrf_up", "lwrf_up", "u10m", "v10m", "uplift", "SRTM_elev", 
#            "pdj", "windspeed", "winddir", "Tes")
  vars <- c("t10m", "EVI", "SRTM_elev")

  # subset data by season  
  spr <- sf[sf$season == "spring",]
  fal <- sf[sf$season == "fall",]
  
  
#---------------------------------plot doy representation to check that it is similar for pres and abs across the dataset
  pres <- ann[ann$presence == 1,]
  abs <- ann[ann$presence == 0,]

  ggplot(data=pres, aes(doy)) + geom_histogram(fill="red", alpha=0.5) + 
    geom_histogram(data=abs, aes(doy), fill="blue", alpha=0.5) + 
    theme_classic() + ggtitle("blue = absent, red = present")
  ggsave(filename=paste0(fig.dir, spp, "/", spp, "_doy_hist.png"), width=4, height=3, dpi=600, units="in")


#----------------------------------pairwise correlation plot
#   df.cor <- ann[,names(ann) %in% vars]
# 
#   out <- paste0(fig.dir, spp, "/", spp, "_variable_correlations.png")
#   png(out, width=12, height=12, units="in", res=600)
#   pairs(df.cor, diag.panel=panel.hist, upper.panel=panel.cor)
#   dev.off()


#----------------------------Random forest - for variable importance
  response <- ann$presence

  #keep just the environmental vars here
  preds.rf <- ann[,c(4:ncol(ann))]
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


#--------------------------------- Generalized linear mixed model with year as random effect
  # using only migration + breeding season data
  #center the data
  zscore = apply(sf[,names(sf)%in% vars ], 2, function(x) {
    y = (x - mean(x))/sd(x)
    return(y)
  })

  #make into data frame and add ID and time vars back in
  zscore <- as.data.frame(zscore)
  zscore$presence <- sf$presence
  zscore$year <- as.factor(sf$year)
  zscore$season <- sf$season

  spring <- zscore[zscore$season == "spring",]
  fall <- zscore[zscore$season == "fall",]

  # run for each season and print results to screen
  fit <- glmer(presence ~ Temp_sfc + swrf_up + swrf + lwrf_up + lwrf + EVI + u10m + uplift + SRTM_elev + 
               (1|year), data = spring, family = "binomial", 
               control = glmerControl(optimizer = "bobyqa"))
  s<-capture.output(summary(fit))
  write(s, file=paste0(fig.dir, spp, "/", spp, "spring_glmer.txt"))


  fit <- glmer(presence ~ Temp_sfc + swrf_up + swrf + lwrf_up + lwrf + EVI + u10m + uplift + SRTM_elev + 
               (1|year), data = fall, family = "binomial", 
               control = glmerControl(optimizer = "bobyqa"))
  s<-capture.output(summary(fit))
  write(s, file=paste0(fig.dir, spp, "/", spp, "fall_glmer.txt"))


#--------------------------------- Generalized additive mixed model with year as random effect
  # using all data to plot how conditions change across flying during the seasons and across years
  #TODO: This part needs work
  daygam <- gamm(doy~t10m + swrf + lwrf + EVI + u10m + v10m + windspeed + pdj + (1|year), data=ann)


#--------------------------------- Plot the data during the migration season
  #subset data by season
  spr_pres <- spr[spr$presence==1,]
  fal_pres <- fal[fal$presence==1,]

  for (i in 1:length(vars)){
  
    p <- ggplot(sf, aes(doy, get(vars[i]))) + geom_point(alpha=0.01) + geom_smooth(col="black") + 
      geom_smooth(data=spr_pres, aes(doy, get(vars[i]), group=year), col="cadetblue", fill="cadetblue") + 
      geom_smooth(data=fal_pres, aes(doy, get(vars[i]), group=year), col="orange", fill="orange") + 
      geom_vline(data=migtime, aes(xintercept=spr_begin), col="cadetblue") +
      geom_vline(data=migtime, aes(xintercept=fal_end), col="orange") +
      theme_classic() + theme(text=element_text(size=14)) + ylab(vars[i]) + ggtitle(spp) +
      scale_x_continuous(breaks = seq(0, 365, by = 50), limits = c(0,365)) 
    ggsave(plot=p, filename=paste0(fig.dir, spp, "/", spp,"_", vars[i],".png"), dpi=600, height=4, width=5, units="in")
  }
}
