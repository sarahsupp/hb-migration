###################################################################################################
# Create matrix plot of correlation matrix
# Circle size represents strength of correlation
# Circle color represents direction of correlation
# From http://gallery.r-enthusiasts.com/graph/Correlation_matrix_circles_152
circle.corr <- function(corr, col=c("black","white"), bg = "white",
                        cex = 1, order = FALSE, title = "", ...){
  
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE),
                                   6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix!")
  n <- nrow(corr)
  m <- ncol(corr)
  
#   ## reorder the variables using principal component analysis
#   if (order) {
#     if(!n==m){
#       stop("The matrix must be squre if order is TRUE!")
#     }
#     x.eigen <- eigen(corr)$vectors[, 1:2]
#     e1 <- x.eigen[, 1]
#     e2 <- x.eigen[, 2]
#     alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
#     corr <- corr[order(alpha), order(alpha)]
#   }
  
  ## set up variable names
  rname <- rownames(corr)
  cname <- colnames(corr)
  if (is.null(rname))
    rname <- 1:n
  if (is.null(cname))
    cname <- 1:m
  rname <- as.character(rname)
  cname <- as.character(cname)
  
  ## calculate label-text width approximately
  par(mar = c(0, 0, 2, 0), bg = "white")
  plot.new()
  plot.window(c(0, m), c(0, n), asp = 1)
  xlabwidth <- max(strwidth(rname, cex = cex))
  ylabwidth <- max(strwidth(cname, cex = cex))
  
  ## set up an empty plot with the appropriate dimensions
  plot.window(c(-xlabwidth + 0.5, m + 0.5), c(0, n + 1 + ylabwidth),
              asp = 1, xlab="", ylab="")
  rect(0.5, 0.5, m + 0.5, n + 0.5, col = bg)  ##background color
  
  ## add variable names and title
  text(rep(-xlabwidth/2, n), n:1, rname, col = "red", cex = cex)
  text(1:m, rep(n + 1 + ylabwidth/2, m), cname, srt = 90, col = "red",
       cex = cex)
  title(title)
  
  ## add grid
  segments(rep(0.5, n + 1), 0.5 + 0:n, rep(m + 0.5, n + 1),
           0.5 + 0:n, col = "gray")
  segments(0.5 + 0:m, rep(0.5, m + 1), 0.5 + 0:m, rep(n + 0.5,
                                                      m), col = "gray")
  
  ## assign circles' fill color
  nc <- length(col)
  if(nc==1)
    bg <- rep(col, n*m)
  else{
    ff <- seq(-1,1, length=nc+1)
    bg2 = rep(0, n * m)
    for (i in 1:(n * m)){
      bg2[i] <- rank(c(ff[2:nc], as.vector(corr)[i]),
                     ties.method = "random")[nc]
    }
    bg <- (col[nc:1])[bg2]
  }
  
  ## plot n*m circles using vector language, suggested by Yihui Xie
  ## the area of circles denotes the absolute value of coefficient
  symbols(rep(1:m, each = n), rep(n:1, m), add = TRUE, inches = F,
          circles = as.vector(sqrt(abs(corr))/2), bg = as.vector(bg))
  
  # Create binary vector. 1 where corr stat is greater than threshold, 0 where not
  zz <- as.vector(ifelse(corr >= 0.75 | corr <= -0.75, 1, 0))
  
  points(rep(1:m, each = n)[zz==1], rep(n:1, m)[zz==1], pch=16, col="yellow")
}
#------------------------

###################################################################################################
# Subset observations by seasonal migration dates based on Sarah's text files of dates per year.
# ann is the annotated observation data. migtime is the table of spring and fall migration dates.
seasonalSub <- function(ann, migtime, yr) {
  # start by subsetting for year, since spr and fall migrations start and end at slightly different 
  # time each year
  ann.yr <- ann[ann$year == yr,]
  #now subset by spring start and end date for the given year (yr)
  spr.begin <- migtime$spr_begin[migtime$year == yr]
  spr.end <- migtime$spr_end[migtime$year == yr]
  ann.spr <- ann.yr[ann.yr$doy >= spr.begin & ann.yr$doy <= spr.end,]
  
  #now the same for fall
  fall.begin <- migtime$fal_begin[migtime$year == yr]
  fall.end <- migtime$fal_end[migtime$year == yr]
  ann.fall <- ann.yr[ann.yr$doy >= fall.begin & ann.yr$doy <= fall.end,]
  
  #add "season" field and put tables together
  ann.spr$season <- "spring"
  ann.fall$season <- "fall"
  
  sf <- rbind(ann.spr, ann.fall)
  sf$season <- as.factor(sf$season)
  
  return(sf)
  
}#end seasonalSub function

###################################################################################################
#For plotting transparency - this function will take a color and set transparency (on 0-255 scale)
#from: http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

###################################################################################################
#Another function for transparency - see if this is better
#from one of the responses here: http://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

###################################################################################################
#Function that submits request to movebank and returns the filepath of the resulting xml.
#requires: url= movebank url, un=username, pw=password, xy=file containing xypoints to be annotated - must be in specific movebank format
#xml=xml file that generates specifics of the request (see hb_movebank_createXML.py), req.outdir=directory where details of request are stored for
#record-keeping
#
submitMB <- function(url, un, pw, xy, xml,req.outdir){
  now <- Sys.time()
  #for file naming
  text.now <- format(now, format="%Y%m%d_%H%M%S")
  ret.xml <- paste(req.outdir, "/", unlist(strsplit(basename(xy), "\\."))[1], "_", unlist(strsplit(basename(xml), "\\."))[1], "_", text.now, ".xml", sep="")
  
  #system call to curl (was able to get this to work somehow when I couldn't get Rcurl to work! - just need to plug ahead for now.)
  curl.cmd <- paste("curl -o ", ret.xml, " -F \"request=@", xml, "\" -F \"tracks=@", xy, "\" -F \"login=", un, "\" -F \"password=", pw, "\" ", url, sep="")    
  system(curl.cmd)
  
  return(list(ret.xml,now))
}#end submitMB

###################################################################################################
#write text files for resubmitting failed mb requests
#requires same psql db set up that I have - this is not a super general function, but 
#could be adapted to use text logs instead of the db. 

#ALSO, requests should not be resubmitted unless you know why it failed and have resolved the issue.
# 
# writeFails <- function()
#   #
#   
#   
#   
#   
###################################################################################################
#clip eBird points by a polygon and returns clipped points.
#requires sp package
#pts is a dataframe of eBird points. sa is a SpatialPolygonsDataFrame representing the study area
#LATfield and LONfield refer to the names of the latitude and longitude fields in the point layer
clipPoints <- function(pts, sa, LONfield, LATfield) {
  #Convert to spatial points df
  coordinates(pts) <- c(LONfield, LATfield)
  
  #clip adata.sp to sa boundaries - reduce the number of points in analysis
  pts.na <- over(pts, sa)
  pts.na2 <- cbind(pts, pts.na)
  pts.clip <- pts.na2[!is.na(pts.na2$ras),]
  
  return(pts.clip)
}#end clipPoints
  

###################################################################################################
# Enhance pairs plots with simple correlations or histograms

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "black", ...)
}

###################################################################################################
#Filter annotated data - we can add to this as necessary. For now, super simple.

fil_ebird <- function(raw_ann, outfile) {
  #filter file
  #filter out all NA values
  fil1 <- raw_ann[complete.cases(raw_ann),]
  #Filter out rows where modis veg pixel reliability is anything but 0 (good data). That means
  #that all pixels used in the generation of this value were high quality (nothing missing, no clouds, no ice).
  fil2 <- fil1[fil1$VegIndex_pix_reliability == 0,]
  
  #write the output file
  write.csv(fil2, outfile, quote=F, row.names=F)
  
  #return fil2 in case we need to do anything else with it
  return(fil2)
  
}#end fil_ebird

###################################################################################################
# Clip points with a polygon
PointsPolyIntersect <- function(points, polygon) {
  data.clip <- over(points, polygon)
  ptid <- na.omit(data.clip) 
  pt.poly <- points[as.numeric(as.character(row.names(ptid))),]  
  return(pt.poly)
}
###################################################################################################
#formatting the eBird absent data date in various ways
formatAbsent <- function(adata) {
  adata$DAY <- format(adata$DAY, format="%j")
  adata$YEAR <- format(adata$YEAR, format="%Y")
  dt <- paste(adata$DAY, adata$YEAR, sep="_")
  d2 <- as.Date(dt, "%j_%Y")
  adata$DATE <- d2
  adata$MONTH <- format(adata$DATE, "%m")
  adata$DOM <- format(adata$DATE, "%d")
  return(adata)
}

#######################################################################################################
#identify season and time windows
ID_windows = function(yeardat, spring, peak, fall, timewindow){
  # This function separates each year into a given number of time windows, where window 2 (w2) begins 
  # on the date for the onset of spring migration and the last window (wn) ends in the time frame that 
  # follows the time frame containing the date for the end of fall migration (wn-1)
  # yeardat: a data.frame subsetted for a single year
  # spring: date for onset of spring migration
  # peak: date for peak latitude for the population
  # fall: date for end of fall migration
  # timewindow: numeric value indicating the time window used to aggregate observation data (e.g. 7 would define a week)
  # returns data for that year labeled and clipped to the spring-breeding-fall migration period (e.g. excludes "winter")
  
  yeardat$increment = 0
  yeardat$window = 0
  yeardat$season="winter"
  
  if(fall < spring | fall < peak) {
    print(paste0("ERROR: migration dates do not make sense, check data ", yeardat$SCI_NAME[1], yeardat$YEAR[1]))
    return(yeardat)
  }
  
  start <-  spring - timewindow
  n.windows <- floor(((fall - start)/timewindow) + 2)
  mid.window <- ceiling((peak - start)/timewindow)
  end <- start + n.windows * timewindow
  
  yeardat[yeardat$yday >= start & yeardat$yday < end,]$increment <- yeardat[yeardat$yday >= start & yeardat$yday < end,]$yday-start + 1
  yeardat[yeardat$yday >= start & yeardat$yday < end,]$window <- ceiling(yeardat[yeardat$yday >= start & yeardat$yday < end,]$increment/timewindow)
  yeardat[yeardat$yday >= start & yeardat$yday < end & yeardat$window < mid.window,]$season <- "spring"
  yeardat[yeardat$yday >= start & yeardat$yday < end & yeardat$window > mid.window,]$season <- "fall"
  yeardat[yeardat$yday >= start & yeardat$yday < end & yeardat$window == mid.window,]$season <- "breeding"
  
  #assign the start ydate for each window (for later aggregation and plotting)
  n=start
  yday1=c(NA,n)
  while(n <= end){
    n=n+5
    yday1=append(yday1, n)
  }
  yday1.df = data.frame("window"=0:length(yday1[-1]), "yday1"=yday1)  
  yeardat = merge(yeardat, yday1.df, by = intersect("window", "window"))
  
  # #assign id var for comparing time frames (id compares rows/locations)
  # if(yeardat$pres[1] %in% c(0,1)){ yeardat$compare.win = yeardat$window }
  # else if(yeardat$pres[1] == -5) { yeardat$compare.win = yeardat$window +1 }
  # else if(yeardat$pres[1] == -10) { yeardat$compare.win = yeardat$window +2 }
  # else if(yeardat$pres[1] == -15) { yeardat$compare.win = yeardat$window +3 }
  # else if(yeardat$pres[1] == 5) { yeardat$compare.win = yeardat$window -1 }
  # else if(yeardat$pres[1] == 10) { yeardat$compare.win = yeardat$window -2 }
  # else if(yeardat$pres[1] == 15) { yeardat$compare.win = yeardat$window -3 }
  
  return(yeardat)
}


#######################################################################################################
#assign time windows using the ID_windows function and label by season using the migration dates

assign_window_season = function(migdates, data, windowlength){
  #migdates is a dataframe with columns: spr_begin, peak_lat, fal_end, spr_spd, fal_spd, species, year
  #data is the main dataframe containing columns for the environmental data associated with hummingbird location points (or absence)
  #windowlength is the number of days in a time window (we started with three-day window frames, but are now using 5-day window frames)
  years = unique(migdates$year)
  for (y in years){
    yrdat = data[data$year == y,]
    yrmig = migdates[migdates$year == y,]
    #assign time frames from the alpha hulls (the last value is the desired time frame)
    yrdat_win = ID_windows(yrdat, spring=yrmig[1,1], peak=yrmig[1,2], fall=yrmig[1,3], timewindow=windowlength)
    
    if (y == years[1]){ assigned_data = yrdat_win }
   else{ assigned_data = rbind(assigned_data, yrdat_win) }
  }
  return(assigned_data)
}


importANDformat = function(path, prescode, migdates, alpha_window){
  #path is a pathname to load the RData file for eBird occurrence + remote sensing data
  #prescode is the code to indicate presence (-1:-3 = previous windows, 1 = present, 0 = absent, 2:4 = next windows)
  #migdates is a dataframe with migration dates for assigning season
  #alpha_window is a number indicating the number of days used to assign alpha hulls (we use 5)
  require(lubridate)
  
  dat = get(load(path))
  dat$pres = as.factor(prescode)
  dat$id = as.factor(row.names(dat))
  #lubridate to pull month and year from filename
  dat$month = as.factor(month(as.Date(dat$timestamp)))
  dat$year = as.factor(year(as.Date(dat$timestamp)))
  #dat$year =  ordered(as.factor(year(as.Date(dat$timestamp))), levels=c(2008:2014))
  dat$yday = as.numeric(yday(as.Date(dat$timestamp)))
  #assign window id codes for later comparisons
  dat = assign_window_season(migdates, dat, alpha_window)
  dat$season = as.factor(dat$season)
  dat$window = as.factor(dat$window)
  dat <- dat[,c("id", "location.long", "location.lat", "window", "year", "month", "yday", "yday1", "season", "EVI", "SRTM_elev", "t10m", "pres")]
  return(dat)
}

cleanColNames = function(migdates){
  #strip extra "" and "\\" from fields
  names(migdates) = gsub('.{1}$', '', substring(names(migdates),3))   #make the column names look nicer
  migdates$spr_begin = as.numeric(gsub("\"", "", migdates$spr_begin, fixed=TRUE))
  migdates$peak_lat = as.numeric(gsub("\"", "", migdates$peak_lat, fixed=TRUE))
  migdates$fal_end = as.numeric(gsub("\"", "", migdates$fal_end, fixed=TRUE))
  migdates$spr_spd = as.numeric(gsub("\"", "", migdates$spr_spd, fixed=TRUE))
  migdates$fal_spd = as.numeric(gsub("\"", "", migdates$fal_spd, fixed=TRUE))
  migdates$species = gsub("\"", "", migdates$species, fixed=TRUE)
  migdates$year = as.factor(gsub("\"", "", migdates$year, fixed=TRUE))
  return(migdates)
}


## From http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_%28ggplot2%29/#Helper%20functions 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=0.95, .drop=TRUE) {
  ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
  ##   data: a data frame.
  ##   measurevar: the name of a column that contains the variable to be summariezed
  ##   groupvars: a vector containing names of columns that contain grouping variables
  ##   na.rm: a boolean that indicates whether to ignore NA's
  ##   conf.interval: the percent range of the confidence interval (default is 95%)
  
  require(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
               )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Multiple plot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


RS_means = function(dat){
  # calculates and returns means, sd, se, and 95% CI for the main RS variables 
  dat.lat = summarySE(dat, measurevar="location.lat", groupvars=c("pres", "yday1", "year", "window"), na.rm=TRUE, conf.interval=0.95)
  dat.lon = summarySE(dat, measurevar="location.long", groupvars=c("pres", "yday1", "year", "window"), na.rm=TRUE, conf.interval=0.95)
  dat.evi = summarySE(dat, measurevar="EVI", groupvars=c("pres", "yday1", "year", "window"), na.rm=TRUE, conf.interval=0.95)
  dat.t10m = summarySE(dat, measurevar="t10m", groupvars=c("pres", "yday1", "year", "window"), na.rm=TRUE, conf.interval=0.95)
  dat.elev = summarySE(dat, measurevar="SRTM_elev", groupvars=c("pres", "yday1", "year", "window"), na.rm=TRUE, conf.interval=0.95)
  
  #merge data.frames, and name columns appropriately (check col numbers)
  dat.sum=cbind(dat.lat, dat.lon[,6:9], dat.evi[,6:9], dat.t10m[,6:9], dat.elev[,6:9])
  names(dat.sum) = c("pres", "yday1", "year", "window", "N", "mean.lat", "sd.lat", "se.lat", "ci.lat",
                     "mean.lon", "sd.lon", "se.lon", "ci.lon", "mean.EVI", "sd.EVI", "se.EVI", "ci.EVI",
                     "mean.t10m", "sd.t10m", "se.t10m", "ci.t10m", "mean.elev", "sd.elev", "se.elev", "ci.elev")
  
  return(dat.sum)
}

#calculate Akaike weights sensu Hobbs & Hilborn 2006
#http://www.planta.cn/forum/files_planta/2006_nobbs_ecology_153.pdf
Akaike.weight <- function(AIC.vec){
  AIC.vec <- AIC.vec-min(AIC.vec,na.rm=T)
  Akai.weight.denom <- sum(exp(-2*AIC.vec))
  Akai.weight <- (exp(-2*AIC.vec)) / Akai.weight.denom
  return(Akai.weight)
}

