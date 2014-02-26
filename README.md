hb-migration
============

This repository is to collect the code for a project to describe and predict the migration of North American hummingbird species for 2004-2013 and into the future. Species include: anna's, allen's, black-chinned, broad-billed, broad-tailed, calliope, costa's, rufous, and ruby-throated hummingbird. Contains code for cleaning eBird data, describing and analyzing species migration paths and timing, exploring environmental and climate correlations with migration, hummingbird physiology, and predicting response to climate change. Collaborators on this project include S. R. Supp, F. A. La Sorte, T. Cormier, M. Lim, G. Bohrer, D. Powers, S. Wethington, and C. H. Graham,.
Code was written by SRS, TC, and FAL.

The code and data in this repository allow for the analyses and figures to be fully replicated using eBird data, icosahedron map of North America, and remote sensing products.

The project and code in this repository are still under development.

**Requirements:**
R 3.x, R packages `knitr`, `geosphere`, `ggplot2`, `ggmap`, `maps`, `mgcv`, `raster`, `reshape`, `spaa`, and files containing functions specific to this code.

The analyses can be replicated by changing the working directory at the top of the file .R to the location on your computer where you have stored the .R, .txt, and .csv files, and by modifying the pathnames for output to where you would like to store the results. Questions about the code should be directed to SRS (sarah@weecology.org).


**Data use:**
Data is currently privately housed, but may be requested for the purposes of replication. If you wish to use the data for additional research, you should contact the Graham lab (catherine.graham@stonybrook.edu) and Cornell Lab of Ornithology to request eBird data (http://ebird.org/ebird/data/download).


**Included Files:**
* `hb-migration.r` script -- cleans up the trait and community datasets, runs descriptive analyses, and outputs figures.
* `migration-fxns.r` script -- stores functions for subsetting, analyzing, and plotting the data.
* `eBird-effort.r` script -- makes maps of total citizen science checklist effort per hex cell for each year. Contains daily data.
* `hb-centroids.r` script -- current a "junk" script. Is not necessary to run analyses.
* `BBLmapcode.r` script -- maps foreign recaptures from Bird Banding Lab data.

**License:** This code is available under a BSD 2-Clause License.
