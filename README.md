hb-migration
============

This repository is to collect the code for projects to describe and predict the migration of North American hummingbird species for 2004-2013 and into the future. Species include: anna's, allen's, black-chinned, broad-billed, broad-tailed, calliope, costa's, rufous, and ruby-throated hummingbird. Contains code for cleaning eBird data, describing and analyzing species migration paths and timing, exploring environmental and climate correlations with migration, hummingbird physiology, and predicting response to climate change. Generates all results and figures for the paper(s) related to the project. We are working with movebank to automate data annotation with movebank for "non-track" data points. Collaborators on this and related projects include S. R. Supp, F. A. La Sorte, T.A. Cormier, M. C.W. Lim, G. Bohrer, D. Powers, S. Wethington, S. Goetz, P.A. Beck, and C. H. Graham.
Code was written by SRS, TC, and modified from FAL and PAB.

The code and data in this repository allow for the analyses and figures to be fully replicated using eBird data, BBL data,  GBIF data, icosahedron map of North America, and remote sensing products.

The project and code in this repository are still under development.

**Requirements:**
R 3.x, R packages `chron`, `fields`, `knitr`, `gamm4` , `geosphere`, `ggplot2`, `ggmap`, `maps`, `maptools`, `mapdata`, `mgcv`, `party`, `plyr`, `randomForest`, `raster`, `reshape2`, `rgdal`, `Rmisc`, `RPostgreSQL`, `SDMTools`,  `sirad`, `sm`, `sp`, `spaa`, `timeDate`, `XML`, Python 2.7 or greater, Python modules `argparse`, `os`, `xml.etree.ElementTRee`, `xml.dom`, and files containing functions specific to this code.

The analyses can be replicated by changing the working directory at the top of the file .R to the location on your computer where you have stored the .R, .txt, and .csv files, and by modifying the pathnames for output to where you would like to store the results. Questions about the code should be directed to SRS (sarah@weecology.org).


**Data use:**
Data is currently privately housed, but may be requested for the purposes of replication. If you wish to use the data for additional research, you should contact the Graham lab (catherine.graham@stonybrook.edu) and Cornell Lab of Ornithology to request eBird data (http://ebird.org/ebird/data/download), or the corresponding author (sarah@weecology.org).

*Datasets needed*
* eBird observation localities with species and date
* eBird list of complete checklists
* eBird effort data by hexagon cell
* icosahedron cell data
* mapdata
* environmental rasters and data from movebank, interpolated to daily values
* Bird Banding Laboratory encounter and recovery dataset
* GBIF data, can be queried using R


**Included Files:**

*Project 1 - Supp et al. 201X. Ecosphere*
* `hb-migration.r` script -- cleans up the trait and community datasets, runs descriptive analyses, and outputs figures.
* `migration-fxns.r` script -- stores functions for subsetting, analyzing, and plotting the data.
* `BBLmapcode.r` script -- maps foreign recaptures from Bird Banding Lab data. M. Lim

*Project 2 - eBird + remote sensing*
* `hb-gbif.r` script -- loops through the species to identify observed wintering locations in Central America from GBIF data online.
* `hb-physiological-demand.r` -- Calculates standard operative temperature (Te), cost of thermoregulation (Ct) and cost of transport/flight (Cf) from the environmental and geographic variables.
* `Standard_Operative_Temperature_from_meteorology.R` -- GET FROM EXTREME_LIMITS: Functions to calculate standard operative temperature from radiation, temperature and geographic variables for individuals.
* `Diffuse_fraction_of_solar_radiation.r` -- GET FROM EXTREME_LIMITS: Functions to calculate diffuse fraction of solar radiation, solar zenith, and other variables needed for incorporating radiation into standard operative temperature measures.
* `aggregate_hb_checklists.r` script -- separates the original queried file into separate files for each species, writes to file
* `hb_movebank_createXML.py` script -- creates xml files for submission to movebank.org for annotation with environmental variables.
* `hb-submit_request_movebank.r` script -- submits annotation request to movebank using curl
* `hb-submit_movebank_request.sh` script -- generates and submits movebank annotation request using point locations and a request xml generated by `hb_movebank_createXML.py`. 
* `hb_movebank_downloadAnnotation_csv.sh` script -- checks status of movebank requests and downloads annotated files when they are available.
* `hb_prep_mb-ebd_tracks.R` script -- prepares the eBird data for annotation with movebank environmental data
* `hb_RS_functions.R` script -- 
* `hb_analysis_envData.R` script --
* `hb_explore_envData.R` script --

*Extra scripts from exploratory or past analyses*
* `eBird-effort.r` script -- makes maps of total citizen science checklist effort per hex cell for each year, using the effort file downloaded from the eBird portal. Contains daily data. (not needed)
* `eBird-effort-new.r` script -- processes the eBird effort data (not needed)
* `hb-centroids.r` script -- current a "junk" script. Is not necessary to run analyses. (not needed)
* `present_figs.r` script -- generated figures for a presentation in May 2014 (not needed)

**License:** This code is available under a BSD 2-Clause License.
