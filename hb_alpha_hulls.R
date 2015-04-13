# Compute weekly alpha hulls from ebird data
# The main code will be within a function s.t. it can be adapted for each species

# Kevin Guay (kguay@whrc.org)
# on 13 Apr 2015
# mo 13 Apr 2015

# load required packages
require(plyr)
require(rworldmap)

# set base path depending on OS (Windows/ Unix)
if(Sys.info()['sysname'] == 'Windows'){
  path.prefix <- 'C:/'
}else{
  path.prefix <- '/mnt/arctic/c/'
}

