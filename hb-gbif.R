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
           "Selasphorus calliope", "Selasphorus rufus")
centralamerica = c("MX", "GT", "SV", "HN", "CR", "PA", "BZ", "NI")

#count all GBIF observation records
occ_count(basisOfRecord='OBSERVATION')

#count all GBIF records with lat long 
occ_count(georeferenced=TRUE)

occ_count(country=mexico_code, georeferenced=T, type="year")

#get data for each hb species
for (s in unique(species)){
  key <- name_backbone(name=s, kingdom='animal')$speciesKey
  dat <- occ_search(taxonKey=key, return='data', limit=3000, country=mexico_code, hasCoordinate=TRUE, decimalLatitude='33,10',
                  fields=c('name', 'key', 'decimalLatitude', 'decimalLongitude', 'elevation'))
  dat$name <- s # so that we just get one symbol type
  dat = dat[dat$decimalLatitude < 35, ] #so it only encompasses the main winter range of species
  gbifmap(input=dat, mapdatabase="world", 
          region=c("Mexico","Honduras","Guatemala","El Salvador","Panama","Belize","Costa Rica","Nicaragua"))
}
