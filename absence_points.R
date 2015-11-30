require(dplyr)
require(gdata)
require(tidyr)

window.size <- 5

dat.path <- paste0("/home/lorra/Dropbox/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/t", window.size, "/")

getAbsenceWindows <- function(test) {
  test$window <- as.numeric(test$window)
  nonzerowindow <- sort(unique(test$window))[-1]
  
  # abs is a window by window matrix showing where absences are taken from for each window
  abswindows <- matrix(FALSE, nrow=length(nonzerowindow), ncol=length(nonzerowindow))
  lowerTriangle(abswindows) <- TRUE
  
  #calculate the window at which peak migration starts (counted as beginning of fall)
  peakwin <- filter(test, DAY==peak) %>%
    select(window) %>%
    unique(.) %>%
    as.numeric(.)
  
  # get the index/rownumber for the peakwindow
  i.peakwin <- which(nonzerowindow==peakwin)
  
  # get rid of all absences earlier than peak for the fall presences. 
  abswindows[(i.peakwin + 1):nrow(abswindows), 1:(i.peakwin - 1)] <- FALSE
  abswindows <- data.frame(window=nonzerowindow, abswindows)
  names(abswindows) <- c("preswindow", nonzerowindow)
  # gather and remove "FALSE" rows to get the full set of absence data
  m.abswindows <- gather(abswindows, "abswindow", "value", -preswindow) %>%
    filter(value==TRUE) %>%
    select(-value)
  
  # merge to the presences file so can get all the info needed
  prsabs <- merge(test, m.abswindows, by.x="window", by.y="abswindow") %>%
    mutate(window=preswindow) %>%
    select(-preswindow)
  
  return(prsabs)
}

for(sp in c("bchu", "bthu", "cahu", "rthu", "ruhu")) {
  if(sp=="rthu") {
    prs <- read.table(paste0(dat.path, "t", window.size, "_", sp, "_humdat_east.txt"), sep=",",  header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  } else {
    prs <- read.table(paste0(dat.path, "t", window.size, "_", sp, "_humdat_west.txt"), sep=",",  header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  }
  prs <- as.data.frame(sapply(prs, function(x) gsub("\\\\", "", x)), stringsAsFactors = FALSE)
  l.prs <- split(prs, prs$YEAR)
  l.abs <- lapply(l.prs, getAbsenceWindows)
  abs <- do.call("rbind", l.abs)
  
  # get window/year start days
  winstart <- group_by(prs, YEAR, window) %>%
    summarise(absday = min(DAY))
  
  # get the start day of the presence window that it corresponds to
  abs <- merge(abs, winstart)
  
  save(abs, file=paste0(sp, "_abs.rda"))
}

# get summaries
require(dplyr)
require(tidyr)
load("bchu_abs.rda")

cahu.summary <- group_by(abs, YEAR, window) %>%
  summarise(count = n()) %>%
  spread(key=YEAR, value=count)

# need to do this for presences too and work out how to display.

for(sp in c("bchu", "bthu", "cahu", "rthu", "ruhu")) {
  if(sp=="rthu") {
    prs <- read.table(paste0(dat.path, "t", window.size, "_", sp, "_humdat_east.txt"), sep=",",  header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  } else {
    prs <- read.table(paste0(dat.path, "t", window.size, "_", sp, "_humdat_west.txt"), sep=",",  header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  }
  prs <- as.data.frame(sapply(prs, function(x) gsub("\\\\", "", x)), stringsAsFactors = FALSE)
  prs$window <- as.numeric(prs$window)
  prs.summary <- group_by(prs, YEAR, window) %>%
    summarise(count = n()) %>%
    spread(key=YEAR, value=count)
  save(prs.summary, file=paste0(sp, "_prs_summary.rda"))
}