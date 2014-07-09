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
