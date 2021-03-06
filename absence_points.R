require(dplyr)
require(gdata)
require(tidyr)
require(ggplot2)
require(cowplot)

window.size <- 5

dat.path <- paste0("/home/lorra/Dropbox/hb_migration_data/ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/t", window.size, "/")
out.path <- paste0("/home/lorra/Dropbox/hb_migration_data/pseudoabsences/t", window.size, "/")
db.path <- "/home/lorra/Dropbox/hb_migration_data"

getAbsenceWindows <- function(test, direction) {
  test$window <- as.numeric(test$window)
  nonzerowindow <- sort(unique(test$window))[-1]

  #calculate the window at which peak migration starts (counted as beginning of fall)
  peakwin <- filter(test, DAY==peak) %>%
    select(window) %>%
    unique(.) %>%
    as.numeric(.)
  
  # get the index/rownumber for the peakwindow
  i.peakwin <- which(nonzerowindow==peakwin)
  
  # abs is a window by window matrix showing where absences are taken from for each window
  abswindows <- matrix(FALSE, nrow=length(nonzerowindow), ncol=length(nonzerowindow))
  
  if(direction=="comingfrom") {
    lowerTriangle(abswindows) <- TRUE
    # get rid of all absences earlier than peak for the fall presences. 
    abswindows[(i.peakwin + 1):nrow(abswindows), 1:(i.peakwin - 1)] <- FALSE
  } else if (direction=="goingto") {
    # get rid of all absences later than peak for the spring presences. 
    upperTriangle(abswindows) <- TRUE
    abswindows[1:(i.peakwin - 1), (i.peakwin + 1):nrow(abswindows)] <- FALSE
  }
  
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
  
  # This section subsamples so that the max number of absences is 3*presences ---------------------------------
  # find how many presences there are by window
  prsXwindow <- group_by(test, window) %>%
    summarise(no.pres=n())
  
  absXwindow <- group_by(prsabs, window) %>%
    summarise(no.abs=n())
  
  prsabsXwindow <- merge(prsXwindow, absXwindow)
  
  prsabsXwindow$sample <- prsabsXwindow$no.abs >= 3*prsabsXwindow$no.pres
  
  l.prsabs <- split(prsabs, prsabs$window)
  
  # sampling function
  sample.abs <- function(wdw) {
    if(prsabsXwindow[prsabsXwindow$window==wdw, "sample"]) {
      out <- l.prsabs[[wdw]][sample(nrow(l.prsabs[[wdw]]), 3*prsabsXwindow[prsabsXwindow$window==wdw, "no.pres"]),]
    } else {
      out <- l.prsabs[[wdw]]
    }
    return(out)
  }
  
  l.prsabs <- lapply(as.character(prsabsXwindow$window), function(x) sample.abs(x))
  prsabs <- do.call("rbind", l.prsabs)
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
  l.abs.comingfrom <- lapply(l.prs, function(x) getAbsenceWindows(x, "comingfrom"))
  abs.comingfrom <- do.call("rbind", l.abs.comingfrom)
  
  l.abs.goingto <- lapply(l.prs, function(x) getAbsenceWindows(x, "goingto"))
  abs.goingto <- do.call("rbind", l.abs.goingto)
  
  # get window/year start days
  winstart <- group_by(prs, YEAR, window) %>%
    summarise(absday = min(DAY))
  
  # get the start day of the presence window that it corresponds to
  abs.comingfrom <- merge(abs.comingfrom, winstart)
  abs.goingto <- merge(abs.goingto, winstart)
  
  write.csv(abs.comingfrom, file=paste0(out.path, "/", sp, "_abs_comingfrom.txt"))
  write.csv(abs.goingto, file=paste0(out.path, "/", sp, "_abs_goingto.txt"))
}

# create unique point/data files for sending to WH
files <- list.files(out.path, full.names = TRUE)
input <- lapply(files, read.csv)
input.l <- do.call("rbind", input)
input.u <- unique(input.l[,c("LATITUDE", "LONGITUDE", "absday", "TIME", "YEAR")])
write.csv(input.u, file=paste0(out.path, "/unique.absence.points.csv"))

test <- unique(input.u[,c("LATITUDE", "LONGITUDE", "absday", "YEAR")])

# code to run after files back from WH - annotates individual species files with environmental variables
env_dat <- read.csv(paste0(db.path, "/gesp_Quality_EVI_SRTM_t10m_250.csv"), stringsAsFactors = FALSE)

# SRTM_elev and t10m were not on the final file so have included previous version of these
other <- read.csv(paste0(db.path, "/gesp_fcl_abs_combined_250.csv"))
other <- select(other, timestamp, location.long, location.lat, SRTM_elev, t10m)
env_dat <- merge(env_dat, other)
env_dat$timestamp <- strptime(env_dat$timestamp, "%Y-%m-%d %H:%M:%S")
if(!file.exists(paste0(out.path, "/annotated"))) dir.create(paste0(out.path, "/annotated"))

files <- list.files(out.path, pattern = "_abs_")

chck_rows <- lapply(files, function(f) {
  sp_dat <- read.csv(paste0(out.path, f), stringsAsFactors = FALSE)
  sp_dat$timestamp <- strptime(paste(as.Date(sp_dat$absday - 1, origin = paste0(sp_dat$YEAR, "-01-01"), value = 1), sp_dat$TIME),
                               "%Y-%m-%d %H:%M:%S")
  
  # annotate species absence points with environmental data
  sp_dat_ann <- merge(env_dat, sp_dat, 
                      by.x = c("location.lat", "location.long", "timestamp"), 
                      by.y = c("LATITUDE", "LONGITUDE", "timestamp")) %>%
    mutate(id=rownames(.), pres=0, year = YEAR, month = MONTH, yday = absday, yday1 = absday) %>%
    dplyr::select(id, location.long, location.lat, window, year, month, yday, yday1, season, EVI, SRTM_elev, t10m, pres)
  
  sp_dat_ann$window <- as.factor(sp_dat_ann$window)
  sp_dat_ann$year <- as.factor(sp_dat_ann$year)
  sp_dat_ann$month <- as.factor(sp_dat_ann$month)
  sp_dat_ann$season <- as.factor(sp_dat_ann$season)
  sp_dat_ann$pres <- as.factor(sp_dat_ann$pres)
  
  save(sp_dat_ann, file=paste0(out.path, "/annotated/", gsub(".txt", ".rda", f)))
  chck <- c(nrow(sp_dat), nrow(sp_dat_ann))
})

# create and plot summaries of presences
prs.summary <- list()
for(sp in c("bchu", "bthu", "cahu", "rthu", "ruhu")) {
  if(sp=="rthu") {
    prs <- read.table(paste0(dat.path, "t", window.size, "_", sp, "_humdat_east.txt"), sep=",",  header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  } else {
    prs <- read.table(paste0(dat.path, "t", window.size, "_", sp, "_humdat_west.txt"), sep=",",  header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  }
  prs <- as.data.frame(sapply(prs, function(x) gsub("\\\\", "", x)), stringsAsFactors = FALSE)
  prs$window <- as.numeric(prs$window)
  prs.summary[[sp]] <- group_by(prs, YEAR, window) %>%
    summarise(count = n()) %>%
    mutate(species=sp)
}

prs.summary.full <- do.call("rbind", prs.summary)

prs.plot <- ggplot(data=prs.summary.full, aes(x=window, y=count, colour=YEAR)) +
  geom_line() +
  facet_wrap(~species, scales="free")

save_plot("~/Dropbox/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/presences.pdf", prs.plot,
          base_aspect_ratio = 1:3)

# create and plot summaries of absences
abs.summary <- list()
for(sp in c("bchu", "bthu", "cahu", "rthu", "ruhu")) {
  
  load(paste0(out.path, "/", sp, "_abs.rda"))
  abs.summary[[sp]] <- group_by(abs, YEAR, window) %>%
    summarise(count = n()) %>%
    mutate(species=sp)
}

abs.summary.full <- do.call("rbind", abs.summary)

abs.plot <- ggplot(data=abs.summary.full, aes(x=window, y=count, colour=YEAR)) +
  geom_line() +
  facet_wrap(~species, scales="free")

save_plot("~/Dropbox/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/figures/absences.pdf", abs.plot,
          base_aspect_ratio = 1:3)
