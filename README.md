hb-migration
============

This repository is to collect the code for a project to describe and predict the migration of North American hummingbird species for 2004-2013 and into the future. Species include: anna's, allen's, black-chinned, broad-billed, broad-tailed, calliope, costa's, rufous, and ruby-throated hummingbird. Contains code for cleaning eBird data, analyzing species migration paths and timing, and predicting response to climate change. Collaborators on this project include S. R. Supp, C. H. Graham, F. A. La Sorte, T. Cormier, D. Powers, S. Wethington.
Code was written by S. R. Supp.

The code and data in this repository allow for the analyses and figures to be fully replicated using eBird data and remote sensing products.

The project and code in this repository are still under development.

**Requirements:**
R 3.x, R packages `knitr`, `ggplot2`, `ggmap`, `reshape`, `vegan`, and files containing functions specific to this code.

The analyses can be replicated by changing the working directory at the top of the file .R to the location on your computer where you have stored the .R, .txt, and .csv files, and by modifying the pathnames for output to where you would like to store the results. Questions about the code should be directed to SRS (sarah@weecology.org).


**Data use:**
Data is currently privately housed, but may be requested for the purposes of replication. If you wish to use the data for additional research, you should contact the Graham lab (catherine.graham@stonybrook.edu) and Cornell Lab of Ornithology to request eBird data (http://ebird.org/ebird/data/download).


**Included Files:**
`hb-migration.r` script -- cleans up the trait and community datasets, runs descriptive analyses, and outputs figures.
`migration-fxns.r` script -- stores functions for subsetting, analyzing, and plotting the data.

**License:** This code is available under a BSD 2-Clause License.
