# Query GBIF data for North American migratory hummingbird species in the winter range:
# Archilochus colubris, A. alexandri, Selasphorus platycercus, S. calliope, and S. rufus
# Sarah R. Supp 2014

library(rgbif)
library(plyr)
library(XML)
library(httr)
library(maps)
library(ggplot2)

species = ("Archilochus colubris", "Archilochus alexandri", "Selasphorus platycercus", 
           "Selasphorus calliope", "Selasphorus rufus")

#count all GBIF observation records
occ_count(basisOfRecord='OBSERVATION')

#count all GBIF records with lat long 
occ_count(georeferenced=TRUE)

#count all records from Mexico
occ_count(country='MEXICO')

head(name_lookup(query='Archilochus alexandri', rank="species", return="data"))

key <- name_backbone(name='Archilochus alexandri', kingdom='animal')$speciesKey
dat <- occ_search(taxonKey=key, return='data', limit=300)
dat$name <- "Puma concolor" # so that we just get one symbol type
gbifmap(input=dat, mapdatabase="world", region="Mexico")

