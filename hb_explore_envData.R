library(randomForest)
library(party)
library(sm)
library(ggplot2)

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
fig.dir <- "/Users/tcormier/Documents/820_Hummingbirds/migration_study/figures/"


###########################################################################################

#ann file
agan.file <- paste0(agan.dir, spp, "/combined/", spp, "_lag",lag,"_allYears.csv")

ann <- read.csv(agan.file, as.is=T)
ann$timestamp <- as.Date(ann$timestamp, format='%Y-%m-%d %H:%M:%S.000')
ann$present <- as.factor(ann$present)
ann$GlobCover <- as.factor(ann$GlobCover)
ann$month <- as.numeric(format(ann$timestamp, "%m"))
ann$year <- as.numeric(format(ann$timestamp, "%Y"))
ann$doy <- as.numeric(format(ann$timestamp, "%j"))

#glob cover labels - legend for map comes from http://dup.esrin.esa.it/files/p68/GLOBCOVER2009_Product_Description_Manual_1.0.pdf
#glob.label <- unique(ann$GlobCover)
#oops, missing some here - go back and figure out which classes are missing.
#glob.names <- c("closed Broadleaved","Closed to Open Grassland","Mosaic Vegetation Cropland","Closed to Open Broadleaved Evergreen and/or Semi-Deciduous Forest",
                "Closed Needle-Leaved Evergreen Forest","Urban","Water","Closed to Open Mixed Broadleaved and Needleleaved Forest","Mosaic Forest or Shrubland and Grassland",
                "Closed to Open Shrubland","Mosaic Cropland/Vegetation", "Closed Broadleaved Forest Regularly Flooded, Fresh Water",
                "Closed Broadleaved Semi-Deciduous and/or Evergreen Forest Regularly Flooded, Saline Water", "Sparse Vegetation",
                "Close to Open Grassland or Shurbland or Woody Vegetation on Regularly Flooded or Waterlogged soil, fresh brakish, or saline water",
                "Open Broadleaved Deciduous Forest/Woodland","Bare","Permanent Snow and Ice","No Data")

#separate pres from abs
pres <- ann[ann$present == 1,]
abs <- ann[ann$present == 0,]

#Separate by season
#first read in the migration timing table and format
migtime.file <- paste0(migtime.dir, "west_migration_", spp, ".txt")
migtime <- read.table(migtime.file, header=T, sep=" ")
years <- unique(ann$year)

#set up dfs that will hold data during spring and fall migrations
sf <- as.data.frame(matrix(data=NA, nrow=0, ncol=ncol(ann)+1))
names(sf) <- c(names(ann), "season")

for (yr in years) {
  seasonal <- seasonalSub(ann, migtime, yr)
  #append this year's data to sf table
  sf <- rbind(sf, seasonal)
}#end yr loop

###########################################################################################
# habitat utilization graphs (i.e. variable histograms)
var <- "t10m"

#subset by presence, absence
pres <- sf[sf$present==1,]
abs <- sf[sf$present==0,]

#subset by spring, fall
spr <- sf[sf$season == "spring",]
fall <- sf[sf$season == "fall",]

#END here on 7/9: refine this into a nice plot!
# Kernel density plots for spring and fall migrations
p <- ggplot(spr, aes(x=get(var), fill=present)) + geom_density(alpha=.3) + scale_fill_brewer()
#p + theme(legend.position="top")

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
circle.corr( yy, order = TRUE, bg = "gray65",
             col = colorRampPalette(c("blue","white","red"))(100) )

circle.corr( pp, order = TRUE, bg = "gray65",
             col = colorRampPalette(c("blue","white","red"))(100) )

circle.corr( aa, order = TRUE, bg = "gray65",
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