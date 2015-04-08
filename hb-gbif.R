# Query GBIF data for North American migratory hummingbird species in the winter range:
# Archilochus colubris, A. alexandri, Selasphorus platycercus, S. calliope, and S. rufus
# tutorial for rgbif can be found @ROpenSci: http://ropensci.org/tutorials/rgbif_tutorial.html
# Sarah R. Supp 2014

devtools::install_github("ropensci/rgbif")
library(rgbif)
library(plyr)
library(XML)
library(httr)
library(maps)
library(ggplot2)


species = c("Archilochus colubris", "Archilochus alexandri", "Selasphorus platycercus", 
           "Stellula calliope", "Selasphorus rufus") #note: must use synonym for "Selasphorus calliope"
centralamerica = c("MX", "GT", "SV", "HN", "CR", "PA", "BZ", "NI")
centralamerica_names = c("Mexico","Honduras","Guatemala","El Salvador","Panama","Belize","Costa Rica","Nicaragua")


#get data for each hb species in winter range (don't include usa), and plot on map of central america
for (s in unique(species)){
  key <- name_backbone(name=s, kingdom='animal')$speciesKey
  locs = data.frame('name'='name', 'key'=1, 'decimalLatitude'=1.0, 'decimalLongitude'=1.0, 'elevation'=1.0)
  
  for (country in unique(centralamerica)){
  dat <- occ_search(taxonKey=key, return='data', limit=3000, country=country, hasCoordinate=TRUE,
                  fields=c('name', 'key', 'decimalLatitude', 'decimalLongitude', 'elevation'))
  
  if(head(dat) != "no data found, try a different search"){ #data doesn't exist for all countries
    if(ncol(dat) < 5){
      dat$elevation = NA
    }
  
    locs = rbind(locs, dat)
    }}
  
  locs=locs[-1,] #delete the first row of dummy data
  locs$name <- s # so that we just get one symbol type
  locs = locs[locs$decimalLatitude < 35 & locs$decimalLatitude > 7, ] #so it only encompasses the main winter range of species - central america
  locs = locs[locs$decimalLongitude < -78 & locs$decimalLongitude > -120, ]
  
  gbifmap(input=locs, mapdatabase="world", region=centralamerica_names)
}
