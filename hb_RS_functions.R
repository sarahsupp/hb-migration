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
  
  