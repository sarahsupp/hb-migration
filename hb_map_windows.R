dat.dir <- "~/Dropbox/hb_migration_data/"

alpha_window <- 5

source("hb_RS_functions.R")
pres_spring <- list()
pres_fall <- list()

# for each species, we are selecting records from 2010. For each season, we select the
# first and last window, as well as those that are 1/4, 1/2 and 3/4 way along the route.
for(sp in c("bchu", "bthu", "cahu", "ruhu", "rthu")) {
  spfiles = list.files(path = paste0(dat.dir, "ebird_annotated_raw/combined/"), pattern = glob2rx(paste0(sp,"*.RData")), recursive=FALSE, full.names=TRUE)
  mfiles = list.files(path = paste0(dat.dir, "ebird_raw/eBird_checklists_2008-2014/aggregate_by_species/"), pattern = glob2rx(paste0(sp,"*migration*")), recursive=FALSE, full.names=TRUE)
  
  #import migration dates from previous analysis
  if(sp == "rthu"){ migdates = read.table(mfiles[1], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")}
  if(sp %in% c("bchu", "bthu", "cahu", "ruhu")){ migdates = read.table(mfiles[2], header=TRUE, sep=",", quote="", fill=TRUE, as.is=TRUE, comment.char="")}
  migdates=cleanColNames(migdates)
  
  #load all the files into a list and remove rows with incomplete predictors
  pres = importANDformat(spfiles[2], 1, migdates, alpha_window)
  pres <- pres[which(complete.cases(subset(pres, select = -yday1))),]
  
  pres$sp <- sp
  
  # windows selection
  p_spring <- subset(pres, season == "spring" & year == 2010)
  p_fall <- subset(pres, season == "fall" & year == 2010)
  
  
  # get 5 evenly spread windows for each sp/ssn combination
  win_list <- as.numeric(as.character(unique(p_fall$window)))
  win_fall <- data.frame(window = c(win_list[ceiling(length(win_list)*0.1)], 
                                    win_list[ceiling(length(win_list)*0.3)], 
                                    win_list[ceiling(length(win_list)*0.5)], 
                                    win_list[ceiling(length(win_list)*0.7)],
                                    win_list[ceiling(length(win_list)*0.9)]), 
                         win_id = as.factor(1:5))
  
  win_list <- as.numeric(as.character(unique(p_spring$window)))
  win_spring <- data.frame(window = c(win_list[ceiling(length(win_list)*0.1)], 
                                      win_list[ceiling(length(win_list)*0.3)], 
                                      win_list[ceiling(length(win_list)*0.5)], 
                                      win_list[ceiling(length(win_list)*0.7)],
                                      win_list[ceiling(length(win_list)*0.9)]), 
                           win_id = as.factor(1:5))
  
  pres_spring[[sp]] <- merge(p_spring, win_spring)
  pres_fall[[sp]] <- merge(p_fall, win_fall)
}



library(ggplot2)
library(plyr)

pres_spring <- ldply(pres_spring)
pres_spring$season <- "Spring"

pres_fall <- ldply(pres_fall)
pres_fall$season <- "Autumn"

pres <- rbind(pres_spring, pres_fall)

# get species scientific names
pres[pres == "bchu"] <- "A. alexandri"
pres[pres == "bthu"] <- "S. platycercus"
pres[pres == "cahu"] <- "S. calliope"
pres[pres == "ruhu"] <- "S. rufus"
pres[pres == "rthu"] <- "A. colubris"

pres$season <- factor(pres$season, levels = c("Spring", "Autumn"))

find_hull <- function(df) df[chull(df$location.long, df$location.lat), ]

hulls <- ddply(pres, c("sp", "season", "win_id"), find_hull)

ggplot() + geom_polygon(data=map_data("world"), aes(x=long, y=lat, group=group), fill="gray") + 
  coord_fixed(ratio=1.3, xlim=c(-130,-50), ylim=c(20,60)) + 
  geom_polygon(data = hulls, aes(x = location.long, y = location.lat, fill = win_id), alpha = 0.5) +
  facet_grid(sp ~ season) + 
  scale_fill_viridis_d(name = "Window") + 
  theme_classic() + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())

ggsave("~/Dropbox/NASA_Hummingbirds/P10_eBird_Migration_multiple topics/2-Mechanisms/manuscript_RSOS/Revisions/hb_window_map.png")

library(dplyr)
win_summary <- group_by(pres, sp, season, win_id) %>% summarise(count = n())
