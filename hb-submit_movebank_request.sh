#!/bin/bash

# Purpose: This shell script uses Rcurl to submit an annotation request to movebank.
# The arguments required include: url=url to which to submit request; xy=a csv file containing a list
# of xy coordinates for annotation (see "/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/test_calliope.csv")
# for an example of how it must be structured); xml=an xml file containing the request parameters (i.e. sensor, product, 
# interpolation etc. See /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/NDVI_EVI.xml for
# example); un=movebank username; pw=movebank password; outdir=output directory for returned xml document with access key.
#
# Author: Tina Cormier
# Date: March 12, 2014
# Status: Runs


if [ "$#" -ne 6 ]; then
    echo "usage: hb-submit_movebank_request.sh <url> <xy> <xml> <un> <pw> <outdir>"
    echo
    echo "example: hb-submit_movebank_request.sh http://www.bioinfo.mpg.de/orn-gateway/request-annotation-xml.jsp /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/test_calliope.csv /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/NDVI_EVI.xml tcormier nohg3ITh /Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests/"
    exit
fi



url=${1}
xy=${2}
xml=${3}
un=${4}
pw=${5}
outdir=${6}
###############################
#hardcoded for testing:
#
#address to submit request
url=http://www.bioinfo.mpg.de/orn-gateway/request-annotation-xml.jsp
#
#Text file of xy locations of points to be annotated. Must have the following columns
#with no spaces: Timestamp, location-long, location-lat, height-above-ellipsoid (optional)
xy=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/test_calliope.csv

#XML containing movebank request details (see hb_movebank_createXML.py)
xml=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/NDVI_EVI.xml

#user name
un=tcormier

#password
pw=nohg3ITh

#outdir
outdir=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests/
outfile=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests/test.xml
###############################

#output file naming
now=$(date +%Y%m%d_%H%M%S)
outfile=$(echo ${outdir}/$(basename ${xy} .csv)_$(basename ${xml} .xml)_request_${now}.xml)

#access key list - hardcoded and is appended to, not overwritten.
ak=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests_accessKeys.txt

#access key request lookup table
ak_lut=/Users/tcormier/Documents/820_Hummingbirds/prelim_analyses/movebank/auto_submit/submitted_requests_accessKeys_lut.csv

#curl request - should return xml with access key
curl -o ${outfile}  -F "request=@${xml}" -F "tracks=@${xy}" -F "login=${un}" -F "password=${pw}" ${url}

#add access key to list - may just be able to parse the lut
#accessKey=$(cat ${outfile} | awk -F "=" '{print $4}' | awk -F " " '{print $1}' | sed 's/"//g')
#echo ${accessKey} >> ${ak}

#append info to access key lut
if [ ! -f ${ak_lut} ]; then
	echo AccessKey,date_time,tracks,xml,status > ${ak_lut}
	echo ${accessKey},${now},${xy},${xml},requested >> ${ak_lut}
else
	echo ${accessKey},${now},${xy},${xml},requested >> ${ak_lut}
fi





