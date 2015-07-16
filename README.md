hb-migration
============

This repository is to collect the code for projects to describe and predict the migration of North American hummingbird species for 2004-2013 and into the future. Species include: anna's, allen's, black-chinned, broad-billed, broad-tailed, calliope, costa's, rufous, and ruby-throated hummingbird. Contains code for cleaning eBird data, describing and analyzing species migration paths and timing, exploring environmental and climate correlations with migration, hummingbird physiology, and predicting response to climate change. Generates all results and figures for the paper(s) related to the project. We are working with movebank to automate data annotation with movebank for "non-track" data points. Collaborators on this and related projects include S. R. Supp, F. A. La Sorte, T.A. Cormier, Kevin C. Guay, M. C.W. Lim, G. Bohrer, D. Powers, S. Wethington, S. Goetz, P.A. Beck, and C. H. Graham.
Code was written by SRS, TC, KG, and modified from FAL and PAB.

The code and data in this repository allow for the analyses and figures to be fully replicated using eBird data, BBL data,  GBIF data, icosahedron map of North America, and remote sensing products.

The project and code in this repository are still under development.

#### Contents
* Overview
  * Data use
  * Datasets needed
  * License
* Scripts
  * Requirements
  * Project 1 - Supp et al. 201X. Ecosphere
  * Project 2 - eBird + remote sensing
  * Extra scripts from exploratory or past analysis
* Query MoveBank (MB)
  * PostgreSQL Database Setup
  * File Structure
  * Submit a Query

#### Overview
##### Data use
Data is currently privately housed, but may be requested for the purposes of replication. If you wish to use the data for additional research, you should contact the Graham lab (catherine.graham@stonybrook.edu) and Cornell Lab of Ornithology to request eBird data (http://ebird.org/ebird/data/download), or the corresponding author (sarah@weecology.org).

##### Datasets needed
* eBird observation localities with species and date
* eBird list of complete checklists
* eBird effort data by hexagon cell
* icosahedron cell data
* mapdata
* environmental rasters and data from movebank, interpolated to daily values
* Bird Banding Laboratory encounter and recovery dataset
* GBIF data; can be queried using R

##### License:
This code is available under a BSD 2-Clause License.

#### Scripts

##### Requirements
R 3.x, R packages `chron`, `fields`, `knitr`, `gamm4` , `geosphere`, `ggplot2`, `ggmap`, `lme4`, `ltm`, `maps`, `maptools`, `mapdata`, `mgcv`, `party`, `plyr`, `randomForest`, `Rarity`, `raster`, `reshape2`, `rgdal`, `Rmisc`, `RPostgreSQL`, `SDMTools`,  `sirad`, `sm`, `sp`, `spaa`, `timeDate`, `XML`, Python 2.7 or greater, Python modules `argparse`, `os`, `xml.etree.ElementTRee`, `xml.dom`, and files containing functions specific to this code.

The analyses can be replicated by changing the working directory at the top of the file .R to the location on your computer where you have stored the .R, .txt, and .csv files, and by modifying the pathnames for output to where you would like to store the results. Questions about the code should be directed to SRS (sarah@weecology.org).

##### Project 1 - Supp et al. 2015. Ecosphere. http://dx.doi.org/10.1890/ES14-00290.1
* `hb-migration.r` script -- cleans up the trait and community datasets, runs descriptive analyses, and outputs figures.
* `migration-fxns.r` script -- stores functions for subsetting, analyzing, and plotting the data.
* `BBL_Appendix.r` script -- replicates the figures in Appendix A for the manuscript.
* `BBLmapcode.r` script -- original script that maps foreign recaptures from Bird Banding Lab data. M. Lim
* `eBird-effort.r` script -- combines eBird files with icosahedron map to calc effort. Use this one.
* `eBird-effort-new.r` script -- combines eBird files with icosahedron map
* `aggregate_hb_checklists.r` -- used to initially combine eBird observations, which came as separate month-year files

##### Project 2 - eBird + remote sensing
* `agg_hb_checklists_paper2.r` script -- aggregates the monthly checklists FAL sent from eBird, and separates them based on species
* `hb-gbif.r` script -- loops through the species to identify observed wintering locations in Central America from GBIF data online.
* `hb-physiological-demand.r` -- Calculates standard operative temperature (Te), cost of thermoregulation (Ct) and cost of transport/flight (Cf) from the environmental and geographic variables.
* `Standard_Operative_Temperature_from_meteorology.R` -- GET FROM EXTREME_LIMITS: Functions to calculate standard operative temperature from radiation, temperature and geographic variables for individuals.
* `Diffuse_fraction_of_solar_radiation.r` -- GET FROM EXTREME_LIMITS: Functions to calculate diffuse fraction of solar radiation, solar zenith, and other variables needed for incorporating radiation into standard operative temperature measures.
* `aggregate_hb_checklists.r` script -- [tc] separates the original queried file into separate files for each species, writes to file. \*similar script (newer) by kg: `hb_mb_combine_datasets.R`
* `hb_prep_mb-ebd_tracks.R` script -- prepares the eBird data for annotation with movebank environmental data
* `hb_RS_functions.R` script -- [tc]
* `hb_analysis_envData.R` script --
* `hb_explore_envData.R` script --

###### Alpha Hulls
* `hb_alpha_hulls.R` script -- [kg] calculates alpha hulls for each species, year and week for the day of observation, preceeding week and following week. Generates csv file to be submitted to movebank.
* `hb_alpha_functions.R` script -- [kg] functions that are useful for the alpha hull calculation script.

###### MoveBank Requests
* `hb_movebank_createXML.py` script -- [tc] creates xml files for submission to movebank.org for annotation with environmental variables.


* `hb-submit_request_movebank.r` script -- [tc] submits annotation request to movebank using curl
* `hb-submit_movebank_request.sh` script [LEGACY] -- [tc] generates and submits movebank annotation request using point locations and a request xml generated by `hb_movebank_createXML.py`.
* `hb_mb_checkStatus_retrieveAnnotation.R` script [tc] checks the status of movebank queries and downloads data when ready.
* `hb_movebank_downloadAnnotation_csv.sh` script -- [tc] checks status of movebank requests and downloads annotated files when they are available.
* `hb_mb_combine_datasets.R` -- [kg] combine annotated data from movebank into a single dataframe/ CSV (with a column for each variable). Output: one CSV and one RData file per species, abs or prs, and time window (focal, minus, plus).

##### Extra scripts from exploratory or past analyses
* `eBird-effort.r` script -- makes maps of total citizen science checklist effort per hex cell for each year, using the effort file downloaded from the eBird portal. Contains daily data. (not needed)
* `eBird-effort-new.r` script -- processes the eBird effort data (not needed)
* `hb-centroids.r` script -- current a "junk" script. Is not necessary to run analyses. (not needed)
* `present_figs.r` script -- generated figures for a presentation in May 2014 (not needed)


#### Query MoveBank (MB)

##### PostgreSQL Database Setup

The `hb-submit_movebank_request.R` and `hb_mb_checkStatus_retrieveAnnotation.R` scripts use a PostgreSQL database to keep track of submitted movebank queries. The following SQL scripts set up the requsite tables.
* `create-table-access_key_list.sql`
* `create-table-access_key_lut.sql`

You will need to change the following line in `hb-submit_movebank_request.R` and `hb_mb_checkStatus_retrieveAnnotation.R` to reflect your database information:

      con <- dbConnect(drv="PostgreSQL", port="5433",host="fusion.whrc.ad", dbname="kguay")

##### File Structure

* `download_annotations` -- MB data is downloaded here
* `env_request_xmls` -- Metadata, used by MB, to identify environmental variables
  * `submit_lists` -- List of `env_request_xmls` to submit
* `status_xmls` -- Status files downloaded from MB. This is how we know when a job is done or failed
* `submit_csv` -- CSV files for each query (e.g. `bchu_abs_movebank_submit_fcl.csv`)
  * `meta` -- Metadata about the requests, including (for each species / year / time window) alpha value, number of samples in hull, number of outliers, total samples
  * `submit_lists` -- Contains csv files describing which files in `submit_csv` to submit to MB
* `submitted_requests` -- XML files used by MB to keep track of submitted requests.

##### Submit a Query

1. Run `hb_alpha_hulls.R`

1. Change the following lines in `hb-submit_movebank_request.R`:

        # The following paths need to be changed (path/to/)

        [17] source("path/to/hb_RS_functions.R")
        [39] xy.csv <- paste0(path/to/movebank/submit_csv/submit_lists/rthu.csv")
        [43] xml.csv <- "path/to/movebank/env_request_xmls/submit_lists/submit_06_24_2015.csv"

        [46] un <- "tcormier"
        [49] pw <- "<password>"

        [53] req.outdir <- paste0("path/to/hummingbird/movebank/submitted_requests/")
        [56] req.outdir <- paste0("path/to/hummingbird/movebank/downloaded_annotations/")

        [62] con <- dbConnect(drv="PostgreSQL", port="5433",host="fusion.whrc.ad", dbname="kguay")

1. Choose which variables you want to query: `movebank/env_request_xmls/submit_lists/list.csv`

        e.g. submit_06_24_2015.csv:
        path/to/movebank/env_request_xmls/Quality_EVI.xml
        path/to/movebank/env_request_xmls/lwrf_swrf_t10m_rh20m.xml
        path/to/movebank/env_request_xmls/MODIS_NPP.xml
        path/to/movebank/env_request_xmls/SRTM_elev.xml

1. Choose which species you want to submit: `movebank/submit_csv/submit_lists/`

        e.g. rthu only: rthu.csv:
        # focal
        path/to/movebank/submit_csv/rthu_abs_movebank_submit_fcl.csv
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_fcl.csv
        # +/- 5 days
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_min_5.csv
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_pls_5.csv
        # +/- 10 days
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_min_10.csv
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_pls_10.csv
        # +/- 15 days
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_min_15.csv
        path/to/movebank/submit_csv/rthu_prs_movebank_submit_pls_15.csv

1. Submit MB queries: run `hb-submit_movebank_request.R`

1. Check the status of the queries: run `hb_mb_checkStatus_retrieveAnnotation.R`

1. Combine individual data products (e.g. NDVI, temp, etc.) into a single dataframe: run `hb_mb_combine_datasets.R`
