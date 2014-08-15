library(randomForest)
#library(party)
library(sm)
library(ggplot2)
library(ltm)
library(Rarity)

source("/Users/tcormier/Documents/scripts/git_repos/hb-migration/hb_RS_functions.R")

###########################################################################################
#spp abbreviation (same one used in file paths and names)
spp <- "ruhu"

#what lag do you want to explore?
lag <- 0

#Aggregated, Annotated File dir
agan.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/movebank/downloaded_annotations/"

#directory containing migration timing text files
migtime.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/data/supp_migration/"

#fig dir 
fig.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/figures/ruhu/"


###########################################################################################

#ann file
agan.file <- paste0(agan.dir, spp, "/combined/", spp, "_lag",lag,"_allYears.csv")

ann <- read.csv(agan.file, as.is=T)
ann$timestamp <- as.Date(ann$timestamp, format='%Y-%m-%d %H:%M:%S.000')
#ann$present <- as.factor(ann$present)
ann$GlobCover <- as.factor(ann$GlobCover)
ann$month <- as.numeric(format(ann$timestamp, "%m"))
ann$year <- as.numeric(format(ann$timestamp, "%Y"))
ann$doy <- as.numeric(format(ann$timestamp, "%j"))

#glob cover labels - legend for map comes from http://dup.esrin.esa.it/files/p68/GLOBCOVER2009_Product_Description_Manual_1.0.pdf
#glob.label <- unique(ann$GlobCover)
#oops, missing some here - go back and figure out which classes are missing.
#glob.names <- c("closed Broadleaved","Closed to Open Grassland","Mosaic Vegetation Cropland","Closed to Open Broadleaved Evergreen and/or Semi-Deciduous Forest",
#                 "Closed Needle-Leaved Evergreen Forest","Urban","Water","Closed to Open Mixed Broadleaved and Needleleaved Forest","Mosaic Forest or Shrubland and Grassland",
#                 "Closed to Open Shrubland","Mosaic Cropland/Vegetation", "Closed Broadleaved Forest Regularly Flooded, Fresh Water",
#                 "Closed Broadleaved Semi-Deciduous and/or Evergreen Forest Regularly Flooded, Saline Water", "Sparse Vegetation",
#                 "Close to Open Grassland or Shurbland or Woody Vegetation on Regularly Flooded or Waterlogged soil, fresh brakish, or saline water",
#                 "Open Broadleaved Deciduous Forest/Woodland","Bare","Permanent Snow and Ice","No Data")

#separate pres from abs
 pres <- ann[ann$present == 1,]
 abs <- ann[ann$present == 0,]

#plot doy representation to check that it is similar for pres and abs
ggplot(data=pres, aes(doy)) + geom_histogram(fill="red", alpha=0.5) + 
  geom_histogram(data=abs, aes(doy), fill="blue", alpha=0.5) + 
  theme_classic() + ggtitle("blue = absent, red = present")


#Separate by season
#first read in the migration timing table and format
migtime.file <- paste0(migtime.dir, "west_migration_", spp, ".txt")
migtime <- read.table(migtime.file, header=T, sep=" ")
years <- unique(ann$year)

#set up dfs that will hold data during spring and fall migrations 
# Note: data that falls outside of seasons is not included in the new dataframe (sf)
sf <- as.data.frame(matrix(data=NA, nrow=0, ncol=ncol(ann)+1))
names(sf) <- c(names(ann), "season")

for (yr in years) {
  seasonal <- seasonalSub(ann, migtime, yr)
  #append this year's data to sf table
  sf <- rbind(sf, seasonal)
}#end yr loop

###########################################################################################
# habitat utilization graphs (i.e. variable histograms)
#test to see if we can get any signal from swrf - this might be wrong to do, but helps visualization of 
#density plots - otherwise, the spike at 0 is so great, it's impossible to see what's happening in the rest
#of the graph.
#convert K to celsius
sf$t10m <- sf$t10m - 273.15

sf2 <- sf
sf <- sf[sf$swrf != 0,]

sf$occurrence[sf$present==1] <- "present"
sf$occurrence[sf$present==0] <- "absent"


#subset by presence, absence
pres <- sf[sf$occurrence=="present",]
abs <- sf[sf$occurrence=="absent",]

#subset by spring, fall
spr <- sf[sf$season == "spring",]
fall <- sf[sf$season == "fall",]

#END here on 7/9: refine this into a nice plot!
# Kernel density plots for spring and fall migrations

#Set up variables and titles for looping
seasons <- c("spr", "fall")
season.titles <- c("Spring", "Fall")

vars <- c("SRTM_elev", "t10m", "EVI", "swrf", "rugosity25", "u10m")
var.titles <- c("Elevation", "Temperature at 10 m", "EVI", "Downward Shortwave Radiation Flux", "Surface Roughness", "East-West Wind at 10 m")
xlab.titles <- c("Elevation (m)", expression("Temperature at 10 m" ~ (degree~C)), "EVI", expression("Downward Shortwave Radiation Flux" ~ (W ~ m^-2)), "Surface Roughness", "East-West Wind at 10 m (m/s)")

#for each variable, let's look at presence vs. absence in spring and fall
for (i in 1:length(vars)) {
  #print(paste(vars[i]))
  for (s in 1:length(seasons)) {
    #print(paste(vars[i], seasons[s]))
    #create graphs for spring presence/absence, spring first
    title <- paste0(season.titles[s], " Rufous Hummingbird Habitat Utilization - \n", var.titles[i])
    outfile <- paste0(fig.dir, season.titles[s], "_ruhu_habitat_utilization_", vars[i], ".pdf")
    spr.col <- c("gray30","cadetblue")
    fall.col <- c("gray30", "orange")
    
    if (seasons[s] == "spr") {
      col <- spr.col
    } else {
      col <- fall.col
    }#end color if
    
    pdf(outfile, width=9, height=8)
    p <- ggplot(get(seasons[s]), aes(x=get(vars[i]),fill=occurrence)) + geom_density(alpha=.4)
    p <- p + scale_fill_manual( values = col)
    p <- p + theme_classic() + theme(text=element_text(size=20))
    p <- p + ggtitle(title) + xlab(xlab.titles[i])
    print(p)
    
    dev.off()
   }#end seasons loop
} #end vars loop
#p + scale_fill_manual( values = c("red", "mediumaquamarine"))


# Now, look at presence only and the difference between seasons on the same plot
for (j in 1:length(vars)) {
  #create graphs of presence spring vs fall
  title <- paste0("Spring vs. Fall Rufous Hummingbird Habitat Utilization - \n", var.titles[j])
  outfile <- paste0(fig.dir, "ruhu_springVsFall_habitat_utilization_", vars[j], ".pdf")
  pdf(outfile, width=9, height=8)
  #FIX COLORS
  p <- ggplot(pres, aes(x=get(vars[j]), fill=season)) + geom_density(alpha=.4)   
  p <- p + scale_fill_manual( values = c("orange", "cadetblue"))
  p <- p + theme_classic() + theme(text=element_text(size=20))
  p <- p + ggtitle(title) + xlab(xlab.titles[j])
  print(p)
  dev.off()

} #end vars loop

###########################################################################################
#pairwise correlation plot
#first, we only want presence/absence and the environmental variables
#use sf2, which still have 0's in swrf
#df.cor <- sf2[,-c(1:4,7,10,12:14,21:24)]
df.cor <- sf2[,c(5,6,9,18,17,15)]
#df.cor$num_season[df.cor$season=="spring"] <- 1
#df.cor$num_season[df.cor$season=="fall"] <- 0
#df.cor <- df.cor[,-12]

#some labels - later, create a lookup table for these so you can add/remove columns and not 
#have to manually go through this.
lab <- c("presence","elev","EVI","t10m","swrf","uplift")
names(df.cor) <- lab
#set color of points with desired transparency
#cols <- makeTransparent("steelblue3", alpha=50)
cols <- addTrans("cadetblue", 50)
#pairwise correlation plot - par(mar=) does not work here bc it's set in
#the corPlot function. Need to adjust at some point.
outcp <- paste0(fig.dir, "Rufous_variable_correlations.png")
png(outcp, width=12, height=12, units="in", res=300)
par(mar=c(0.5, 1, 0.5, 0.5),pch=20, cex.lab=1.5)
corPlot(df.cor, method="pearson", xlab=lab, ylab=lab, col=cols)
dev.off()


#corr.pa <- biserial.cor(sf2$swrf, sf2$present, use="complete.obs", level=2)


###########################################################################################
# Correlate graphs, omitting NA values
ann.sub <- ann[,c(6,8,9,11,15:ncol(ann))]
corr.names <- c("elev","NDVI","EVI","Rugosity","uplift","lwrf","swrf","temp10m","uwind10m","vwind10m")
yy <- cor(ann.sub, use = "pairwise.complete.obs")

pres.sub <- pres[,c(6,8,9,11,15:ncol(pres))]
pp <- cor(pres.sub, use = "pairwise.complete.obs")

abs.sub <- abs[,c(6,8,9,11,15:ncol(abs))]
aa <- cor(abs.sub, use = "pairwise.complete.obs")

rownames(yy) <- corr.names
colnames(yy) <- corr.names
rownames(pp) <- corr.names
colnames(pp) <- corr.names
rownames(aa) <- corr.names
colnames(aa) <- corr.names

# Matrix plot of variable correlations ordered by principle components
circle.corr( yy, order = TRUE, bg = "white",
             col = colorRampPalette(c("blue","white","red"))(100) )

circle.corr( pp, order = TRUE, bg = "white",
             col = colorRampPalette(c("blue","white","red"))(100) )

circle.corr( aa, order = TRUE, bg = "white",
             col = colorRampPalette(c("blue","white","red"))(100) )

###########################################################################################
#Random forest - for variable importance
head(ann)

#LATER, need to filter NDVI/EVI for quality - just looking quickly at this now.
#first remove columns with bad preds (for some reason, rugosity5, slopeu, and slovev produced all NAs)
preds <- ann[,-c(4,7,12:14)]
preds <- na.omit(preds)
# preds$month <- as.numeric(format(preds$timestamp, "%m"))
# preds$year <- as.numeric(format(preds$timestamp, "%Y"))
# preds$doy <- as.numeric(format(preds$timestamp, "%j"))
preds <- preds[,-1]
response <- preds$present
#preds <- preds[,-c(preds$present)]


preds.rf <- preds[,-3]
rf <- randomForest(preds.rf, response, na.rm=T, importance=T)
#decided to use variable imp plot based on mean decrease in GINI bc Calle & Urrea (2010) show its rankings are more stable.
#Also, got explanation from here: http://www.stat.berkeley.edu/~breiman/RandomForests/cc%5Fhome.htm#varimp
varImp <- importance(rf,type = 2)

#FIX TITLES TO INCLUDE SPP AUTOMATICALLY 
rf.graph <- paste0(fig.dir, spp, "/", spp, "_varImp_allvars.pdf")
pdf(file = rf.graph, width=8, height=8)
par(mar=c(8,6,6,4))
varImpPlot(rf, type=2,pch=16, col="blue",main=paste("Black-chinned Hummingbird \nRandom Forest Variable Importance", sep=' '))
#abline(v=abs(min(varImp)), col='red',lty='longdash', lwd=2)
dev.off()


#now do RF without location and time to see which env vars show up
preds.rf2 <- preds[,-c(1:3,15:17)] 
rf2 <- randomForest(preds.rf2, response, na.rm=T, importance=T)
varImp2 <- importance(rf2,type = 2)

#FIX TITLES TO INCLUDE SPP AUTOMATICALLY
rf2.graph <- paste0(fig.dir, spp, "/", spp, "_varImp_envVars.pdf")
pdf(file = rf2.graph, width=8, height=8 )
par(mar=c(8,6,6,4))
varImpPlot(rf2, type=2,pch=16, col="blue",main=paste("Black-chinned Hummingbird \nRandom Forest Variable Importance", sep=' '))
dev.off()
###########################################################################################
#Now try it with the party package based on this blog http://alandgraf.blogspot.com/2012/07/random-forest-variable-importance.html and the Strobl et al. paper)
#WOW, this ran for almost 24 hrs and crashed my machine. Need to rethink.
# data.controls <- cforest_unbiased(ntree=1000, mtry=round(sqrt(ncol(preds))))
# cf1 <- cforest(present~.,data=preds,control=data.controls)
# varimp(cf1)
# varimp(cf1,conditional=TRUE)


#Variables can be considered informative and important if their variable importance value is above the absolute value of the lowest negative-scoring variable. 
#“The rationale for this rule of thumb is that the importance of irrelevant variables varies randomly around zero” (Strobl et al. 2009b: 342). 
# barplot(sort(varImp), horiz=TRUE, xlab="Variable Importance \n(predictors to right of red dashed are significant)", ))
# abline(v=abs(min(varImp)), col='red',lty='longdash', lwd=2)
# 
# #Variables can be considered informative and important if their variable importance value is above the absolute value of the lowest negative-scoring variable. 
# #“The rationale for this rule of thumb is that the importance of irrelevant variables varies randomly around zero” (Strobl et al. 2009b: 342). 
# dotplot(sort(data.cforest.varimp), xlab=”Variable Importance
#         in DATA\n(predictors to right of dashed vertical line are
#                   significant)”, panel = function(x,y){
#                     panel.dotplot(x, y, col=’darkblue’, pch=16, cex=1.1)
#                     panel.abline(v=abs(min(data.cforest.varimp)), col=’red’,
#                                  lty=’longdash’, lwd=2)
###########################################################################################

# Compare environmental data using Kolmogorov-Smirnoff (KS) tests
# KS test has less power to detect a shift in the median but more power to detect changes in 
# the shape of the distributions. Null hypothesis is that both groups were sampled from populations 
# with identical distributions. It tests for any violation of that null hypothesis -- 
# different medians, different variances, or different distributions.

critical_D <- function(n1, n2){
  # http://www.soest.hawaii.edu/wessel/courses/gg313/Critical_KS.pdf
  # assuming alpha = 0.05
  # values > Da are significant?
  Da <- 1.36 * sqrt((n1+n2)/(n1*n2))
  return(Da)
}

seasons <- c("spring", "fall")
vars <- c("SRTM_elev", "EVI", "lwrf", "swrf", "t10m") #add other vars as necessary

# Compare presence vs absence points in seasons and years for each variable
ks_pa <- data.frame("year"=1, "season"=NA, "var"=NA, "Dstat"=1, "Pvalue"=1)
for (y in unique(years)){
  for (s in unique(seasons)){
    pres <- sf[sf$present == 1 & sf$year == y & sf$season == s,]
    abs <- sf[sf$present == 0 & sf$year == y & sf$season == s,]
    for (var in unique(vars)){
      compare <- ggplot(pres, aes(x=get(var))) + geom_density(alpha=0.3, fill="red") + 
        geom_density(data=abs, alpha=0.3, fill="blue") + theme_classic() + xlab(var) + 
        ggtitle(paste(y,s,var, sep=" - "))
      print(compare)
      ks <- ks.test(pres[,colnames(pres) %in% var], abs[,colnames(abs) %in% var])
      ks_pa = rbind(ks_pa, c(y,s,var,round(as.numeric(ks$statistic),4),round(ks$p.value,4)))
    }  
  }
}
ks_pa <- ks_pa[-1,] #delete first row of dummy data


# compare years. Are there any years that really differ?
ks_yrs <- data.frame("year1"=1, "year2"=1, "var"=NA, "Dstat"=1, "pvalue"=1, "signif"=NA)
for (var in unique(vars)){
  for (y in unique(years)){
    if(y==2013)
      next
    for(ynext in unique(years)){
      if (y > ynext)
        next
      if(y==ynext)
        next
    pres <- sf[sf$present == 1 & sf$year == y,]
    abs <- sf[sf$present == 0 & sf$year == y,]
    pres2 <- sf[sf$present == 1 & sf$year == ynext,]
    abs2 <- sf[sf$present == 0 & sf$year == ynext,]
    compare <- ggplot(pres, aes(x=get(var))) + geom_density(alpha=0.3, fill="red") + 
      geom_density(data=pres2, alpha=0.3, fill="blue") + 
      geom_density(data=abs, alpha=0.3, fill="grey60",col="red") + 
      geom_density(data=abs2, alpha=0.3, fill="grey60",col="blue") + 
      theme_classic() + xlab(var) + 
      ggtitle(paste(y,ynext,var, sep=" - "))
    print(compare)
    ks <- ks.test(pres[,colnames(pres) %in% var], pres2[,colnames(pres2) %in% var])
    Da <- critical_D(nrow(pres), nrow(pres2))
    if(Da > ks$statistic){ sig <- "Y" }
    else{ sig <- "N" }  
    ks_yrs <- rbind(ks_yrs, c(y, ynext, var, round(as.numeric(ks$statistic),4), round(ks$p.value,4), sig))
  }}
}
ks_yrs <- ks_yrs[-1,]


# compare seasons Are hb selecting diff in diff seasons? How influenced by background diffs is this?
ks_season <- data.frame("year"=1, "var"=NA, "Dstat"=1, "pvalue"=1, "signif"=NA)
for (var in unique(vars)){
  for (y in unique(years)){
      spr_pres <- sf[sf$present == 1 & sf$year == y & sf$season == "spring",]
      fal_pres <- sf[sf$present == 1 & sf$year == y & sf$season == "fall",]
      compare <- ggplot(spr_pres, aes(x=get(var))) + geom_density(alpha=0.3, fill="cadetblue") + 
        geom_density(data=fal_pres, alpha=0.3, fill="orange") + 
        theme_classic() + xlab(var) + ggtitle(paste(y, var, sep=" - "))
      print(compare)
      ks <- ks.test(spr_pres[,colnames(spr_pres) %in% var], fal_pres[,colnames(fal_pres) %in% var])
      Da <- critical_D(nrow(pres), nrow(pres2))
      if(Da > ks$statistic){ sig <- "N" }
      else{ sig <- "Y" }  
      ks_season <- rbind(ks_season, c(y, var, round(as.numeric(ks$statistic),4), round(ks$p.value,4), sig))
    }
}
ks_season <- ks_season[-1,]



#---------------------------------------------------------------------------------
#         GLM test comparing distribution of envr. data with presence
#---------------------------------------------------------------------------------

for (y in unique(years)){
  for (s in unique(seasons)){
    data <- sf[sf$year == y & sf$season == s,]
    
    fit <- glm(present ~ swrf + lwrf + SRTM_elev + EVI + t10m, family="binomial", data=data)
    
    #check the residual deviance and degrees of freedom in the model (should not be significant)
    check <- 1 - pchisq(fit$deviance, fit$df.residual)
    if (check > 0.05) {
      print (summary(fit))
    }
    else { print("WE NEED A BETTER MODEL") }
  }
}




